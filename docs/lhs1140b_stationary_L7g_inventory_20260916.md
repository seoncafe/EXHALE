# L7g data inventory: where the association and recombination energy goes

Item L7g of `docs/PLAN_20260913_lhs_stationary.md`, the data inventory that
`docs/PLAN_20260916_rev3.md` section 7 requires before any model is written.
No source file was edited by this item and nothing was built. The predecessor
is `docs/lhs1140b_stationary_L7g_20260915.md`, the ledger audit, which closed
the double-count question and left exactly one open: whether
`h2_vibrational_heat_fraction` should hold back more of the 4.478 eV that R15
releases.

Every number below is marked MEASURED (computed here, from the tables and
expressions the tree holds, or from arithmetic on values it holds) or READ
(taken from a source document or from a source file).

---

## 0. Verdict

**The data support the bounded two-parameter approximation, not a
level-resolved cascade.** Three of the four ingredients a cascade needs are
missing, and one of the three is missing inside the tree itself:

1. **REVISED 2026-09-16: the nascent state distribution IS now read, and it
   settles the split.** Orel (1987) for M = H, Esposito & Capitelli (2009) for
   M = H and Paolini, Ohlinger & Forrey (2011) for M = He all resolve the
   product (v, J) and all put the nascent molecule in the highest bound states
   of the ladder; Schwenke (1988) for M = H2 gives the total rate only.
   MEASURED on the code's own ladder, the levels they name leave
   D0 - E_int = 0.0000 to 0.0180 eV, **so f_trans is below 0.01 of D0 and the
   code's f_trans = 0 is data-supported rather than merely conservative**
   (section 1.6). Hollenbach & McKee (1979), the source of the expression the
   code uses, still gives a partition only for GRAIN formation and for the H-
   route, and none for the gas-phase three-body channel. **What the four do
   NOT give is a tabulated distribution**: Orel's state-resolved numbers are
   "available on request", Esposito & Capitelli's are in figures, and both are
   needed as an initial condition for a cascade.
2. **For M = H a calculation of the right scope now exists and has been read,
   but the rate coefficients themselves are not distributed with it; for
   M = H2 and M = He nothing of the right scope is in hand at all.**
   Bossion, Scribano, Lique & Parlant (2018), MNRAS 480, 3718, supplied by the
   user 2026-09-16 as `references/Bossion_2018_MNRAS_480_3718.pdf` and read in
   the published version, covers **exactly 260 of this code's own 302 ladder
   levels** (MEASURED: their listed envelope of highest states per v reproduces
   260 of the 302 levels of `molecular_infrared_data.f90` to the level, and
   they state "A total of 260 initial states were taken into account in the
   calculations"), reaching 49931.7 K = 4.3028 eV, over translational
   temperatures 100 to 15000 K, **with the three-body collisional dissociation
   computed from each initial (v, j)**. Its inverse is the nascent
   distribution. **But the paper distributes no table**: it carries eight
   figures and one room-temperature comparison table, no supporting
   information, no repository and no data-availability statement (READ,
   verified over the whole text). The next-best obtained set, Wrathmall,
   Gusdorf & Flower (2007), covers the lowest 108 of 318 bound levels, which on
   this code's ladder reaches 32711.4 K (2.819 eV, MEASURED), while three-body
   association forms H2 within a few kT of 4.478 eV. Lique (2015) and
   Balakrishnan, Forrey & Dalgarno (1999), both supplied by the user
   2026-09-16 and read in the published version, DO distribute their numbers
   but cover far less of the ladder: **MEASURED on the distributed file
   itself, Lique covers 54 of the code's 302 levels to 1.8881 eV (v = 0 to 3),
   Jozwiak et al. (2024) for helium 53 of them to 1.8451 eV (v = 0 to 3), and
   Balakrishnan et al. (1999), which Jozwiak supersedes, 20 of them to
   1.105 eV. No set that distributes its numbers reaches past v = 3, and the
   nascent molecule is born at v = 10 to 14** (sections 2, 2.1, 2.2a).
3. **The radiative network the tree holds is incomplete exactly where a
   cascade would start.** MEASURED on `molecular_infrared_data.f90`: the line
   list is a connected network over the 302-level X ladder (every one of the
   1833 lines has BOTH its upper and its lower level in the ladder, largest
   residual 0.099 K), but 15 levels near dissociation carry no downward
   transition at all, and they are all of v = 13 and v = 14 plus (10,0),
   (11,0), (11,2), (12,0) and (12,2). A cascade solver built on this table
   alone would trap population in precisely the levels a nascent molecule is
   born in (section 3). Bossion et al. (2018) reach only five of those fifteen
   ((10,0), (11,0), (11,2), (12,0), (12,2), MEASURED), so **the ten levels of
   v = 13 and v = 14 have neither a collisional rate nor a downward radiative
   transition in anything now available**.

**What the inventory does settle, and it settles the open question of the
predecessor memo in this layer's favor.** The code's present validity claim,
that the collider density is far above n_cr of EVERY level of the cascade, was
supported at the code site by the v = 1 critical density alone. It can now be
supported by an ALL-LEVEL bound measured from the code's own A values:
the largest total spontaneous decay rate over all 302 bound levels is
**A_tot,max = 5.594e-06 s^-1** (MEASURED), 6.6 times the v = 1, J = 0 value of
8.532e-07 s^-1 that the fraction uses. That maximum is confirmed from a second,
independent source: Cloudy's `transprob_X.dat`, the Wolniewicz, Simbotin &
Dalgarno (1998) set, gives **5.5945e-06 s^-1 at (1,28)**, agreeing to four
digits (MEASURED, section 3).

**REVISED 2026-09-16, and the ARGUMENT changes although the answer does not.**
The bound was first carried on the H2 channel alone: A_tot,max through the
code's own Hollenbach & McKee fit gives n_cr <= 5.5e+07 cm^-3 at 808 K, against
the layer's n(H2) = 4e9 to 5e13 cm^-3 (READ), and the heat fraction would then
be 1 to within 1.4e-02 at the thinnest point. **Measured H2-H2 rates retire
that form of the bound.** Le Bourlot, Pineau des Forets & Flower (1999), read
from the Cloudy tables, put the thermal v = 1 to v' = 0 rate 21 to 107 times
BELOW the code's `gamma_10_H2` over this layer, so the H2-only n_cr is really
**5.9e+09 cm^-3 at 808 K and 1.8e+08 at 1527 K**, and a pure-H2 collider at the
thin end (n(H2) = 4e9, 1527 K) would give **f_quench = 0.958, not 1**
(MEASURED, sections 2.4 and 2.5).

**f_quench = 1 survives on the FULL collider mix, and atomic hydrogen carries
it.** MEASURED at cell 1 with the best rate for each collider (Lique 2015 for
H, Le Bourlot et al. 1999 for H2, Lee or Jozwiak for He), the collider sum is
**99.46 per cent atomic H, 0.53 per cent He and 0.006 per cent H2**, giving
q = 8.0e+02 s^-1 against A_tot,max and **1 - f_quench = 7.0e-09**; at 1527 K,
f_quench > 1 - 1e-3 needs only n(HI) > 3.0e+08 cm^-3. The remaining
approximation still runs the safe way: the collisional coefficient used is the
v = 1 one, while de-excitation of the closely spaced high levels is faster.

**What Bossion et al. (2018) changes, now that it has been read in the
published version: the reason the cascade is not buildable, not the verdict.**
Before, the obstruction was that no calculation of the right scope existed in
obtainable form. It does exist, its scope is adequate for M = H, and its
scope is now known exactly (260 of 302 levels, 100-15000 K, state-resolved
dissociation included). What obstructs the cascade now is (a) the numbers are
not published with the paper, (b) M = H2 and M = He are untouched by it,
(c) their treatment of the quasi-bound molecules assigns to the DISSOCIATIVE
channel exactly the near-dissociation population a detailed-balance inversion
would need as the nascent distribution (READ: "quasi-bound H2 molecules are
counted as pertaining to the three-body dissociative channel. Around 5000 K
they contribute up to 50 per cent of the dissociative channel, and their
contribution decreases with temperature to reach 30 per cent at 15 000 K"),
and (d) at this layer's 200 to 1600 K the set sits at the edge of its own
accuracy: their trajectories start at a collision energy of 0.1 eV = 1160 K
(READ), their stated mean deviation against the quantum calculation above a
few hundred K is "of the order of a factor 3" (READ), and their
room-temperature v = 1 to v = 0 rate is 1.28e-13 against the measured
3.0 +/- 1.5e-13 cm^3 s^-1 of Heidner & Kasper (1972), their Table 1 (READ).
The recommendation of section 6 is therefore unchanged.

**So the split between the prompt translational share and the internal share
does not move this layer's answer at all.** Both shares end as heat where
f_quench = 1. The split matters only for a shallower base or the top of a
thinner layer, which is what the code comment already says and what this
inventory now supports with data rather than with an assertion.

**The one thing the inventory finds that is NOT degenerate is elsewhere in the
network: R6.** READ from Larsson, McCall & Orel (2008), the published review
in `references/Larsson_2008CPL_462_145.pdf`, p. 149: "the peak in the
vibrational distribution of H2(v) was found to occur at v = 5-6 [45], which
agrees very well with the TSR experiment [56]." On this code's own ladder
(MEASURED) v = 5 lies at 2.2927 eV and v = 6 at 2.6664 eV above (0,0), so
**25 to 29 per cent of R6's 9.250 eV is born as H2 internal energy, not as
kinetic energy of the fragments**, and R6 carries the proton loop, 49 per cent
of the H2 loss at cell 1 (READ from the predecessor memo). In this layer it
makes no difference, for the same reason: f_quench = 1. Where f_quench < 1 it
is the largest unbooked internal share in the network, larger than R15's.

---

## 1. Energy partition of `H + H + M -> H2 + M` (R15, 4.478 eV)

### 1.1 What the code does now

READ, `molecular_reaction_heat.f90` lines 469 to 498: the whole
q(R15) = 4.478 eV is multiplied by `h2_vibrational_heat_fraction(T, n_HI,
n_H2)`, which is Hollenbach & McKee's `(1 + n_cr/n)^-1` built from the v = 1
decay rate A_10 = 8.3e-07 s^-1 and the v = 1 collisional coefficients of their
eq. (6.29), with helium deliberately absent
(`h2_vibrational_relaxation.f90`, module header).

### 1.2 What Hollenbach & McKee (1979) actually gives

READ from `references/Hollenbach_1979_ApJS_41_555.pdf`, p. 586, section VI c
(scanned text, so the arithmetic below is quoted as figures rather than
letter-for-letter): grain formation distributes the 4.48 eV of binding energy
as about 4.2 eV into vibrational and rotational excitation of the molecule,
about 0.2 eV into thermal energy of the molecule leaving the grain, and a
negligible amount to the grain. Their eq. (6.43) carries that split
explicitly, `0.2 + 4.2 [1 + n_cr(H2)/n]^-1` eV per formation: **only the
internal 4.2 eV is subject to the branching, the 0.2 eV is deposited
promptly**. For the H- route, eq. (6.44), they state that the entire 3.53 eV
is assumed to be excitation energy because that route matters as a heat source
only where collisional de-excitation delivers it as heat anyway.

**They give nothing for gas-phase `H + H + M`.** Their section VI c names
grain formation and the H- route (their eq. 6.42) as the two formation
mechanisms, and the three-body channel does not appear. Row R16 of rev3
section 0 is confirmed.

Two further details of that page, READ, that bear on the code:
- their n_cr for the heat fraction, eq. (6.45), is built from `gamma_10^H2`
  and **`gamma_20^H`** (their own words: "we used the rates y10H2 and y20H
  given in equation (6.29) and assumed that the molecule generally radiates in
  steps dv = 2 with A = 10^-6 s^-1"), whereas the code uses `gamma_10^H` for
  the atomic collider. The two atomic fits differ in shape:
  `gamma_10^H = 1.0e-12 T^1/2 exp[-(1000/T)]` against
  `gamma_20^H = 1.6e-12 T^1/2 exp[-(400/T)^2]` (eq. 6.29, READ).
  The code's choice is the slower of the two at these temperatures and so is
  the conservative one.
- their A for the heat fraction is 1e-06 s^-1 for a dv = 2 step, not the
  8.3e-07 s^-1 of the v = 1 total decay.

### 1.3 What the other published sources in hand give

| source | what it gives for R15 | product states? |
|---|---|---|
| Cohen & Westberg (1983) JPCRD 12, 531, p. 559, `references/Cohen_1983JPCRD_12_531.pdf` | the total rate the code uses: k1(H2) = 2.8e-31 T^-0.6 cm^6 s^-1 over 50-5000 K, and separately k1(H) = 8.8e-33 cm^6 s^-1, temperature independent, over the same range (READ) | none |
| Flower & Harris (2007) MNRAS 377, 705, `references/Flower_2007_MNRAS_377_705.pdf` | the total rate only, k3 = 1.44e-26 T^-1.54 cm^6 s^-1, obtained by detailed balance from the collisional dissociation rate through the Saha relation (READ, their section 3.2 and eq. 12). They DO solve level-resolved rate equations, 49 bound levels to about 20000 K, but the three-body channel enters as a total and no nascent distribution over levels is stated | none |
| Glover & Jappsen (2007) ApJ 666, 1, p. 17 sec. 3.3, `references/Glover_2007_ApJ_666_1.pdf` | no three-body formation heating at all; for the H-, H2+ and grain routes, verbatim (READ): "We assume that essentially all of this energy goes into rotational and vibrational excitation of the resulting H2 molecule and, hence, is radiated away at low gas densities and is converted by collisional de-excitation into heat at high gas densities." | none |
| Glover & Abel (2008) MNRAS 388, 1627, `references/Glover_2008_MNRAS_388_1627.pdf` | three-body formation is present as a rate (their reactions 30 and 31, with k31 = k30/8 following Palla et al. 1983) and the heating treatment is referred to Glover & Jappsen (2007). Their section 3.3.3 varies k30 between the Abel et al. (2002) and the Flower & Harris (2007) values and reports the effect on the temperature (READ) | none |
| Palla, Salpeter & Stahler (1983) ApJ 271, 632, `references/Palla_1983_ApJ_271_632.pdf` | the rate and the third-body ratio k(H2) = k(H)/8 that the later work carries | none |
| Lepp & Shull (1983) ApJ 270, 578, `references/Lepp_1983_ApJ_270_578.pdf` | level-resolved collisional DISSOCIATION kinetics, from which the collisional dissociation rate coefficient the detailed-balance route uses is built | none for recombination |
| Yelle (2004) Icarus 170, 167, p. 170, `references/Yelle_2004Icarus_170_167.pdf` | the state of practice this code follows: "Chemical heating rates for the reactions in Table 1 are calculated from heats of formation" (READ), with no product-state partition anywhere | none |

**Superseded 2026-09-16 by section 1.6.** Until the four gas-phase
calculations arrived, the only published partitions were for OTHER formation
routes and bracketed the prompt translational share between 0 and 0.045 of D0
(Glover & Jappsen's "essentially all" into internal, and Hollenbach & McKee's
0.2 of 4.48 eV for grain formation). That transferred bound is no longer
needed: the gas-phase three-body product states are now read directly.

### 1.4 Does the third-body identity matter

**For the rate, yes, and the code is currently applying the M = H2
coefficient to the total density.** READ, `mol_rates.f90` R15:
`rk_R15_3body_H2(T, n) = k3b_H_H_to_H2(T) * n` with n the TOTAL particle
density and `k3b = 2.8e-31 T^-0.6`, which is Cohen & Westberg's k1(H2).
MEASURED from their p. 559 table:
- at 1000 K, k1(H2) = 4.44e-33 and k1(H) = 8.8e-33 cm^6 s^-1, so atomic H is
  **2.0 times** the efficiency of H2, not 1/8 of it as the Palla et al. ratio
  used in the primordial literature would give at this temperature;
- for a monatomic inert third body their second data sheet gives
  k1(Ar) = 7.0e11 T^-1 L^2 mol^-2 s^-1 over 77-5000 K, which is
  1.93e-30 T^-1 cm^6 s^-1 (MEASURED conversion), i.e. 1.93e-33 at 1000 K,
  **0.43 of k1(H2)**. Helium, the dominant third body in this configuration,
  has no entry in either sheet; argon is its nearest tabulated analogue and
  sits below H2.

So applying k1(H2) to the total density overstates the helium contribution by
roughly a factor 2 and understates the atomic-hydrogen one by roughly a factor
2, and the two partly cancel. **No published source in hand gives a helium
third-body coefficient for this reaction.** It is on the request list.

**For the product states, unknown.** The resonance-complex mechanism the
state-resolved calculations use makes the stabilizing collision with M the
step that fixes the nascent level, so the identity of M is expected to matter,
and the two calculations that would say so (Orel 1987 for M = H, Schwenke 1988
for M = H2) could not be obtained.

### 1.5 Temperature range

The layer runs 200 to 1600 K (READ: ghost 226.02 K, cell 1 808.3 K, cell 25
1527.3 K, from the predecessor memo section 6). Cohen & Westberg's rate is
recommended over 50-5000 K and covers it outright. Hollenbach & McKee's
collisional fits are stated good to 20 per cent above about 500 K and 50 per
cent above about 300 K (READ, quoted at the code site), so the ghost row at
226 K sits below both statements; it carries almost no H2 destruction and
reformation, but the fact is recorded.

---

### 1.6 THE NASCENT DISTRIBUTION, READ FROM THE FOUR GAS-PHASE CALCULATIONS

All four were supplied by the user 2026-09-16 and read in the published
version. **Together they settle question 1 qualitatively and give a
data-supported bound on f_trans for the first time.**

| source | third body | resolves product states? | temperature range | what it gives |
|---|---|---|---|---|
| **Orel (1987) JCP 87, 314**, `references/Orel_1987JCP_87_314.pdf` | **H** | **YES**, individual (v, J) | cross sections at E/k = 50, 100, 300, 1000, 2000 K; the state-resolved rate coefficients are plotted to 350 K only (their Figs. 5 to 7) | resonance complex theory with quasiclassical trajectories from six orbiting resonances, (14,5), (14,4), (13,8), (12,11), (12,12), (11,13). The state-resolved numbers are NOT tabulated: "The full results are available on request" |
| **Esposito & Capitelli (2009) JPCA 113, 15307**, `references/Esposito_2009_JPCA_113_15307.pdf` | **H** | **YES**, final (w, k) | **300 to 10000 K** | quasiclassical three-body dynamics plus orbiting resonance theory, with detailed balance applied to direct dissociation from bound states. Distributions over final vibrational states at 300, 1000 and 10000 K and over final rotational states at 300 and 3000 K, all in figures, none tabulated |
| **Schwenke (1988) JCP 89, 2076**, `references/Schwenke_1988JCP_89_2076.pdf` | **H2**, the third body the code's own coefficient belongs to | **NO**, total rate only | 100 to 5000 K | a new ab initio H4 surface and full four-body quasiclassical trajectories within resonance complex theory; Tables I to VIII are surface parameters and Table IX is the total rate constant. No product-state table anywhere |
| **Paolini, Ohlinger & Forrey (2011) PRA 83, 042713**, `references/Paolini_2011PRA_83_042713.pdf` | **He** and Ar | **YES**, per (v, j) | figures to 1000 K; compared with experiment at 77 and 300 K | quantum calculation on a Sturmian basis, resonant and nonresonant contributions on equal footing. **The only published helium third-body calculation in hand** |

**What they say about where the molecule is born, verbatim:**

- Orel (1987), section III: "the distribution is peaked near this line, that
  is, in the highest bound states of the H2 manifold. This distribution
  broadens and extends to lower vibrational states as the collision energy
  increases", and "v = 14, v = 13, v = 12, and v = 11 make up the largest
  contribution to the rate coefficient", with "v = 14 dominating the rate
  constant for low temperature, but becoming overshadowed by v = 13 at higher
  T". Within each v the top rotational levels carry it: for v = 13 "the
  highest four states making up approximately 80% of the total rate
  coefficients". His conclusion: "We have shown that this distribution is
  peaked in the highest bound states of the hydrogen molecule. Due to the
  highly excited nature of the products, the nascent distribution will be
  greatly changed by subsequent collisions."
- Esposito & Capitelli (2009): "It is clear that at low temperature only high
  lying w states are significant, but this changes gradually when the
  temperature increases, with much more similar weights among rates relative
  to all w for T = 10 000 K", and "For T = 300 K, there is a clear prevalence
  of high w, which is less pronounced for T = 1000 K. At T = 10 000 K, the
  distribution is roughly an inverted parabola with a maximum at w = 7."
  Rotation is decisive: "with variation of 2 or 3 orders of magnitude among
  k = 0 and the maximum rotational quantum number compatible with the selected
  final vibration."
- Paolini et al. (2011), on helium: "These results include recombination to
  either of the two most weakly bound vibrational levels for each rotational
  level j <= 20. Recombinations to more strongly bound levels were found to
  make a negligible contribution."

**THE f_trans BOUND, MEASURED on the code's own ladder.** The prompt
translational share of one association event is `D0 - E_int(nascent)`. Taking
the levels these calculations name:

| nascent level | E_int [eV] | D0 - E_int [eV] | share of D0 |
|---|---|---|---|
| highest bound level of v = 11, (11,12) | 4.4601 | 0.0180 | 0.40 % |
| highest bound level of v = 12, (12,10) | 4.4764 | 0.0017 | 0.04 % |
| highest bound level of v = 13, (13,7) | 4.4717 | 0.0064 | 0.14 % |
| highest bound level of v = 14, (14,4) | 4.4781 | 0.0000 | 0.00 % |
| the two most weakly bound levels overall, Paolini's helium channel | 4.4781, 4.4764 | 0.0000, 0.0017 | 0.00 %, 0.04 % |
| the J = 0 origin of v = 11, an extreme lower bound on E_int | 4.0775 | 0.4006 | **8.95 %** |

**So f_trans <= 0.09 with certainty, and the physically indicated value is
below 0.01**, since both Orel and Esposito & Capitelli put the population at
the top J of each v and Paolini's helium channel is the two most weakly bound
levels outright. **The code's present term, which puts all of D0 into the
internal share (f_trans = 0), is therefore DATA-SUPPORTED and not merely the
conservative choice.**

Two consequences that follow and are recorded rather than acted on:

1. **The two-parameter form collapses to one parameter.** With f_trans ~ 0 the
   whole 4.478 eV is internal and f_quench governs all of it, which is exactly
   what `molecular_reaction_heat.f90` computes. What remains unsupported is
   the level resolution of f_quench, not the split.
2. **The nascent molecule is born within 0.02 eV of the dissociation limit, so
   it is also the most easily re-dissociated.** Orel says so ("the nascent
   distribution will be greatly changed by subsequent collisions"), Schwenke
   names the same competition ("the primary bottleneck occurs in the
   relaxation of bound vibrational states close to the dissociation limit"),
   and Bossion et al. book the quasi-bound population into the dissociative
   channel for the same reason. **The effective association rate and the
   effective deposited heat are therefore coupled through the same cascade**,
   and a model that resolves one has to resolve the other.

**Cross-checks between the four and against the code's rate, MEASURED or
READ:**
- Esposito & Capitelli reproduce Orel at 300 K: 1.6e-32 cm^6 s^-1 against
  Orel's 1.19e-32, "obtained both with the six QB 'RBC states'" (READ).
- Their direct-dissociation-plus-detailed-balance contribution runs from
  1.19e-32 cm^6 s^-1 at 300 K to 8.53e-33 at 10000 K (READ), against the
  code's Cohen & Westberg k1(H) of 8.8e-33, temperature independent: a factor
  1.35 at 300 K (MEASURED).
- Schwenke's own verdict on his M = H2 result: "the theoretical results about
  a factor of 2 too small over the temperature range 300-5000 K" against
  experiment, and at 100 K "too small by a factor of about 3" (READ). His
  reference for the experimental recommendation is Cohen & Westberg (1983),
  the same evaluation the code carries, so the code's rate is the measured
  one and Schwenke's is the theory that underestimates it.

## 2. State-resolved collisional de-excitation of H2(v, J) by H, H2 and He

Nothing state-resolved is in the tree. MEASURED by search: no LAMDA-format or
BASECOL-format collisional file exists anywhere under
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/`; the only H2 level file with
(v, J) labels is `src/utils/h2_shielding_lbl/PDR71/data/Levels/level_h2.dat`,
302 levels, which carries energies and degeneracies and no rates. The whole
collisional input of the code for this purpose is the two v = 1 fits of
Hollenbach & McKee eq. (6.29).

| source | collider | levels covered | temperature range | obtained |
|---|---|---|---|---|
| Wrathmall, Gusdorf & Flower (2007) MNRAS 382, 133 | H | "the energetically lowest 54 levels of ortho- and para-H2, that is, 108 of the total of 318 bound levels" (READ); on this code's ladder the 108th level lies at 32711.4 K = 2.819 eV (MEASURED) | to 6000 K (READ) | YES, `references/Wrathmall_2007_MNRAS_382_133.pdf` |
| **Jozwiak, Thibault, Viel, Wcislo & Lique (2024) A&A 685, A113** | **He** | 1059 transitions between levels with internal energies below 15000 cm^-1 (READ). **MEASURED against this code's ladder: 53 of the 302 levels, to 21411.2 K = 1.8451 eV, v = 0 to 3.** Quantum close coupling on the Thibault et al. (2017) surface, revising Flower et al. (1998) | **20 to 8000 K** (READ) | **YES**, `references/Jozwiak_2024A&A_685_A113.pdf`, supplied by the user 2026-09-16. **The numbers ARE distributed**, verbatim: "The ortho- and para-H2 data are available at the CDS via anonymous ftp to cdsarc.cds.unistra.fr (130.79.128.5) or via https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/685/A113", and "The complete dataset with all state-to-state coefficients will be available online from the BASECOL (Dubernet et al. 2013) website" |
| **Lique (2015) MNRAS 453, 810** | H | "We computed cross-sections for all the ro-vibrational levels with internal energies eps_vj < 22 000 K", with the basis excluding levels above 45000 K and levels with j > 24 (READ, their section 2). **MEASURED against this code's ladder: 55 of the 302 levels, reaching 21941.9 K = 1.891 eV, v = 0 to 4 only.** Nearly exact time-independent quantum, hydrogen exchange included rigorously; NO three-body dissociation | 100 to 5000 K (READ) | **YES, paper AND data.** `references/Lique_2015_MNRAS_453_810.pdf` and the supplementary file `references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat`, both supplied by the user 2026-09-16. See section 2.1 for what the file contains |
| **Bossion, Scribano, Lique & Parlant (2018) MNRAS 480, 3718** | H | 260 initial (v, j) states with internal energy at or below 50000 K, the highest per v being (0,30), (1,29), (2,26), (3,25), (4,24), (5,22), (6,20), (7,18), (8,16), (9,14), (10,11), (11,8), (12,4) (READ, their section 2). **MEASURED against this code's ladder: those 260 states ARE 260 of the code's own 302 levels, reaching 49931.7 K = 4.3028 eV; the 42 left out are 1 to 6 high-J levels of each v plus ALL of v = 13 (8 levels) and ALL of v = 14 (5 levels).** Inelastic, reactive (with ortho-para conversion) and **state-resolved three-body collisional dissociation** from each initial (v, j) | 100 to 15000 K translational; 20 collision energies from 0.1 eV (1160 K) to 2.7 eV (31330 K) (READ) | YES, the PAPER, `references/Bossion_2018_MNRAS_480_3718.pdf`, supplied by the user 2026-09-16. **The RATE COEFFICIENTS are NOT distributed with it**: eight figures and one room-temperature comparison table, no supporting information, no repository, no data-availability statement (READ, whole text) |
| **Balakrishnan, Forrey & Dalgarno (1999) ApJ 514, 520** | **He** | their Table 1 lists the levels used, v = 0 with j = 0 to 9 and v = 1 with j = 0 to 9, "selected ro-vibrational levels below the v = 2 excitation threshold" (READ). **MEASURED against this code's ladder: 20 of the 302 levels, reaching 12817.3 K = 1.105 eV.** Close-coupled quantum with full ro-vibrational coupling. Tabulated: pure rotational excitation within v = 0 (their Table 2), pure rotational de-excitation within v = 0 (Table 3, with a fifth-order fit in Table 4), and v = 0 to v = 1 ro-vibrational excitation for eight channels (Table 5, fit in Table 6) | 10 to 5000 K for the pure rotational de-excitation; the fits are stated valid 100-5000 K (rotational) and 100-3000 K (ro-vibrational) (READ) | **YES**, `references/Balakrishnan_1999_ApJ_514_520.pdf`, supplied by the user 2026-09-16. The numbers are in the paper's own tables |
| Flower & Roueff (1999) JPhB 32, 3399, and the companion Letters | H2 | rovibrational transitions induced by ground-state H2 | not stated in the abstract | NO, request list |

**Verdict for question 2, revised 2026-09-16 on the published Bossion et al.
(2018): for M = H the calculation exists and its scope is adequate, but the
numbers are not published; for M = H2 and M = He nothing of the right scope is
in hand.** Wrathmall et al. stop at 2.82 eV while three-body association forms
H2 within a few kT of 4.478 eV. Bossion et al. reach 4.3028 eV over 260 of
this code's own 302 levels, and because their three-body collisional
dissociation is resolved from each initial (v, j), **its detailed-balance
inverse IS the state-resolved three-body association rate into each (v, j) for
M = H**, which answers question 1 as well as question 2. What stands in the
way is not the physics but the distribution of the numbers.

Four limits of that set, all READ from the published paper and all relevant
here:

- **The rate coefficients are not in the paper.** Section 2 describes the
  method, section 3 shows the behavior in eight figures, and the only tabulated
  numbers are the five room-temperature v = 1 to v = 0 thermal rates of their
  Table 1. There is no supporting information, no repository and no
  data-availability statement anywhere in the text. The cooling function built
  on these rates is deferred to "Coppola et al., in preparation".
- **The quasi-bound molecules are booked into the dissociative channel**, and
  they are exactly the near-dissociation population the inversion needs:
  verbatim, "In the present calculations, quasi-bound H2 molecules are counted
  as pertaining to the three-body dissociative channel. Around 5000 K they
  contribute up to 50 per cent of the dissociative channel, and their
  contribution decreases with temperature to reach 30 per cent at 15 000 K."
  An inversion would therefore put into the continuum a large share of what a
  dense layer would see as bound, vibrationally hot H2. The paper states the
  condition for the assumption: "The above assumption remains valid as long as
  the quasi-bound molecules have a lifetime lower than the typical collision
  time. In high-density media, collisions are more frequent and a part of the
  quasi-bound molecules may interact before dissociating." **This layer is a
  high-density medium by that measure**, n(H2) = 4e9 to 5e13 cm^-3.
- **v = 13 and v = 14 are not initial states**, so the inversion gives no
  formation into the 13 levels of the code's ladder that belong to them.
- **At 200 to 1600 K the set is at the edge of its accuracy.** Trajectories
  were run at 20 collision energies from 0.1 eV (1160 K) upward, while eq. (6)
  integrates a Maxwellian from zero; the stated mean deviation against the
  quantum calculation above a few hundred K is "of the order of a factor 3";
  and their own Table 1 gives 1.28e-13 cm^3 s^-1 for the room-temperature
  v = 1 to v = 0 relaxation against the measured 3.0 +/- 1.5e-13 of Heidner &
  Kasper (1972), with the paper stating that "The temperature of the
  experiments (300 K) is out of the accuracy domain of the QCT calculations."
  This layer runs 200 to 1600 K.

### 2.1 What Lique (2015) changes: the collisional side of f_quench, measured

Lique's calculation is the accurate one for H2-H at these temperatures and it
supersedes Wrathmall et al. for the levels both cover, but it covers FEWER of
them: 55 of the code's 302 levels against Wrathmall's 108, because Lique cuts
at an internal energy of 22000 K while Wrathmall counts the 108 lowest levels
whatever their energy (MEASURED).

**The size of the revision against Wrathmall et al. (2007), READ from the
published paper:**
- pure rotational within v = 0, their Table 1: "The maximum deviations are
  less than 20-30 per cent";
- ro-vibrational, their Fig. 4: "The differences are much more significant
  than for the pure rotational excitation. In addition, unlike pure rotational
  excitation, the maximum deviation is at low temperatures. The difference can
  be as high as an order of magnitude", with Wrathmall et al. UNDERESTIMATING
  the ro-vibrational relaxation at low temperature because they omit the
  reactive channels, and generally OVERESTIMATING it at high temperature.

**And it gives a direct correction factor on the collisional coefficient the
code actually uses.** Their Table 2 gives the thermal v = 1 to v' = 0
vibrational relaxation by H at room temperature: 1.8e-13 cm^3 s^-1 (this
work), against the measured 3.0 +/- 1.5e-13 of Heidner & Kasper (1972),
0.7e-13 for Wrathmall et al., 0.8e-13 for Martin & Mandy (1995) and 1.1e-13
for Garcia & Lagana (1986) (READ).

**THE SUPPLEMENTARY DATA FILE IS NOW IN THE TREE, so the correction is
measured across the layer's own temperature range rather than at one point.**
`references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat`, supplied by the user
2026-09-16. Its format, READ from its own header: a title block naming the
paper, then one line per transition, `v  j  v' j'` followed by fifty rate
coefficients, "k(cm3 s-1) (T) , T= 100 to 5000K by steps of 100K".
MEASURED on the file:

- **1431 state-to-state transitions.**
- **Level envelope: v = 0 with j <= 17, v = 1 with j <= 14, v = 2 with j <= 11,
  v = 3 with j <= 8**, initial and final alike. On this code's ladder that is
  **54 of the 302 levels, reaching 21911.0 K = 1.8881 eV**. The distributed
  file therefore stops one manifold short of the paper's stated criterion
  (internal energies below 22000 K would admit v = 4, j = 0 at 21941.9 K);
  v = 4 is absent from the file.
- **Temperature grid: 100 to 5000 K in steps of 100 K, 50 points.**

Reconstructing Lique's own eq. (5) from the file, a Boltzmann average over
j = 0 to 8 of v = 1 summed into j' = 0 to 12 of v' = 0, gives
**1.7809e-13 cm^3 s^-1 at 300 K**, which reproduces the 1.8e-13 of his
Table 2 (MEASURED, and an independent check that the summation is the one he
performed).

**The correction factor on the code's `gamma_10_H`, MEASURED across the
layer:**

| T [K] | Lique, thermal v = 1 -> 0 [cm^3 s^-1] | code `gamma_10_H` | code / Lique |
|---|---|---|---|
| 200 | 8.39e-14 | 9.53e-14 | 1.14 |
| 300 | 1.78e-13 | 6.18e-13 | **3.47** |
| 400 | 4.80e-13 | 1.64e-12 | 3.42 |
| 600 | 1.89e-12 | 4.63e-12 | 2.45 |
| 800 | 4.29e-12 | 8.10e-12 | 1.89 |
| 1000 | 7.49e-12 | 1.16e-11 | 1.55 |
| 1200 | 1.13e-11 | 1.51e-11 | 1.33 |
| 1400 | 1.57e-11 | 1.83e-11 | 1.17 |
| 1600 | 2.03e-11 | 2.14e-11 | 1.05 |
| 2000 | 3.04e-11 | 2.71e-11 | 0.89 |
| 5000 | 1.02e-10 | 5.79e-11 | 0.57 |

**Over this layer's 200 to 1600 K the code's atomic-hydrogen coefficient is
1.05 to 3.47 times the accurate quantum value, largest near 300 K.** The ratio
crosses 1 between 1700 and 1800 K, so above that the same fit UNDERSTATES the
rate; at 100 K it understates it by a factor 100, far outside anything this
layer reaches. Both readings push f_quench the same way as before and neither
threatens the bound: with the measured correction, n_cr(all-level, against
atomic H) moves from 6.78e+05 to **1.27e+06** cm^-3 at 808 K and from 2.76e+05
to **3.01e+05** at 1527 K (MEASURED), three to four orders below the layer's
densities. The conservative bound of section 0 is built on the H2 channel
(n_cr <= 5.5e+07 at 808 K) and is untouched, since Lique gives no H2-H2 rate.

### 2.2 What Balakrishnan et al. (1999) changes: the helium collider, measured

The code leaves helium out of the collider sum by design, and the module
header states that this "understates the de-excitation rate, so it errs
towards radiating the energy away rather than towards claiming it as heat".
**That is directionally right and quantitatively negligible for the
vibrational channel, which is the channel f_quench is about.**

MEASURED from their Table 5 by detailed balance (their Table 1 gives the level
energies; (1,0) lies 4161 cm^-1 = 5986.8 K above (0,0), which reproduces this
code's own ladder value of 5987.0 K to 0.2 K, an independent check on the
ladder):

| T [K] | He, k(1,0 -> 0,0) and the summed partner channel | code's `gamma_10_H2` | code's `gamma_10_H` | He/H2 |
|---|---|---|---|---|
| 700 | 6.8e-17 | 6.70e-14 | 6.34e-12 | 1.0e-03 |
| 1000 | 3.1e-16, summed >= 5.5e-16 | 1.89e-13 | 1.16e-11 | 2.9e-03 |
| 1500 | 1.6e-15, summed >= 3.4e-15 | 6.37e-13 | 1.99e-11 | 5.3e-03 |
| 3000 | 2.1e-14, summed >= 5.3e-14 | 4.40e-12 | 3.92e-11 | 1.2e-02 |

Even at n(He) = 2.6 times n(HI) + n(H2), helium adds 7.5e-03 of the H2
channel's contribution to the collider sum at 1000 K on this estimate.

**CORRECTED 2026-09-16 on the distributed Jozwiak data (section 2.2b): that
number was built from ONE state-to-state channel and understates the thermal
v = 1 to v = 0 total by about an order of magnitude. The conclusion that
helium is negligible survives, and is stronger for a different reason: the
collider sum at the base is 99.7 per cent ATOMIC HYDROGEN, not H2.**

**For PURE ROTATIONAL transitions the picture is the opposite**, and the paper
says so of its own Table 2 (READ): the He rates "are typically several times
larger than those for collisions with H". That matters for a rotational
cascade and not for the vibrational branching f_quench represents.

**What it does not do: reach the levels that matter.** Twenty levels, v = 0
and v = 1 only, to 1.105 eV. The nascent molecule is born at v = 10 to 14.

### 2.2a Jozwiak et al. (2024) supersedes Balakrishnan et al. for the helium
collider, and distributes its numbers

Supplied by the user 2026-09-16, `references/Jozwiak_2024A&A_685_A113.pdf`,
read in the published version. It is the current H2-He set and it replaces
both Flower et al. (1998) and Balakrishnan et al. (1999) for this purpose:

- **1059 transitions between levels with internal energies below 15000 cm^-1**
  (READ). MEASURED against this code's ladder that is **53 of the 302 levels,
  to 21411.2 K = 1.8451 eV, v = 0 to 3**, against Balakrishnan's 20 levels to
  1.105 eV.
- **Temperatures 20 to 8000 K**, against Balakrishnan's 10 to 5000 K.
- **The numbers ARE distributed**, verbatim: "The ortho- and para-H2 data are
  available at the CDS via anonymous ftp to cdsarc.cds.unistra.fr
  (130.79.128.5) or via
  https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/685/A113", and "The complete
  dataset with all state-to-state coefficients will be available online from
  the BASECOL (Dubernet et al. 2013) website". The CDS table is not in this
  workspace and was not fetched by this item.
- It supersedes the earlier work where the two differ, in their own words
  (READ): the agreement with previous calculations is good for pure rotational
  transitions between low-lying levels, but "we do find significant
  discrepancies for rovibrational processes involving highly-excited
  rotational and vibrational states", attributed to the broader range of
  intramolecular distances and the accuracy of the surface.

### 2.2b The Jozwiak CDS tables, read and measured

The CDS tables themselves were supplied by the user 2026-09-16 and are at
`references/Jozwiak_2024J_A+A_685_A113/`. Format, READ from the `ReadMe`
(CDS catalogue J/A+A/685/A113):

- `oh2-lev.dat` and `ph2-lev.dat`, the level lists, three fields per record:
  level label, rotational quantum number J, vibrational quantum number v.
  **26 ortho levels and 27 para levels, 53 in all**, being every H2 level with
  internal energy below 15000 cm^-1.
- `oh2-rat.dat` and `ph2-rat.dat`, one record per transition: transition label,
  initial and final level labels, then the rate coefficients in cm^3 s^-1.
  **520 ortho and 539 para records, 1059 in all**, out of the 676 and 729
  possible; the rest did not converge to the 10 per cent criterion. Both
  directions and the elastic entries are present.
- **Temperature grid: 43 points, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120,
  140, 160, 180, 200, 250, 300, ... , 1000, 1100, ... , 1500, 1750, 2000,
  3000, 4000, 5000, 6000, 7000, 8000 K.**

MEASURED against the code's own ladder:

- **All 53 levels are levels of this code's 302-level ladder**, no exceptions.
- **Envelope: v = 0 with J <= 17, v = 1 with J <= 14, v = 2 with J <= 11,
  v = 3 with J <= 7.** Highest internal energy **21411.2 K = 1.8451 eV**.
- **v = 11 to 14 are ABSENT.** The set does not touch the levels a nascent
  molecule occupies, and it stops essentially where Lique's H set stops.
- Internal consistency check on the file: for the para (1,0) and (0,0) pair the
  ratio of the upward to the downward rate reproduces exp(-dE/T) to four digits
  at 200, 500, 1000 and 2000 K, with dE = 5987.0 K taken from the code's own
  ladder (MEASURED). The tables obey detailed balance and the parsing is right.

**The helium de-excitation rate, measured from the distributed tables.** The
thermal v = 1 to v' = 0 rate, Boltzmann-averaged over the J of v = 1 with the
ortho and para weights and summed into every J' of v = 0, is the quantity to
compare with the code's `gamma_10_H2`:

| T [K] | 200 | 400 | 808 | 1000 | 1500 |
|---|---|---|---|---|---|
| Jozwiak He, thermal v = 1 -> 0 [cm^3 s^-1] | 3.53e-18 | 6.87e-17 | 2.64e-15 | 7.02e-15 | 4.87e-14 |
| n(He) k_He / n(H2) `gamma_10_H2` at He/H = 2.13 | 2.0e-03 | 9.5e-03 | 5.6e-02 | 7.9e-02 | 1.6e-01 |

This supersedes the single-channel estimate of section 2.2, which was about an
order of magnitude low. It also revises Balakrishnan et al. (1999) downward:
the single channel (1,0) to (0,0) is **3.74e-17 cm^3 s^-1 at 1000 K here
against 3.11e-16 obtained by detailed balance from Balakrishnan's Table 5**, a
factor 8.3 (MEASURED), which is the kind of difference the paper announces for
"rovibrational processes involving highly-excited rotational and vibrational
states".

**And against the collider sum the code actually forms, helium is 0.1 per
cent.** MEASURED at cell 1 of the LHS 1140 b molecular base (T = 808.3 K,
x2 = 0.3553, He/H = 2.13, all READ from the predecessor memo and the case
input), per hydrogen nucleus:

| term | value |
|---|---|
| n(HI) `gamma_10_H` | 5.32e-12, **99.66 per cent of the coded sum** |
| n(H2) `gamma_10_H2` | 1.80e-14, 0.34 per cent |
| what helium would add, n(He) k_He | 5.63e-15, **0.106 per cent of the coded sum** |

**So the module header's statement that omitting helium "understates the
de-excitation rate" is right in sign and worth one part in a thousand at the
base**, because the collider sum there is dominated by atomic hydrogen and not
by H2 at all. Helium reaches a sixth of the H2 CHANNEL by 1500 K, but the H2
channel is itself a third of a per cent of the total where the base sits.

**What it does not do: reach the levels that matter.** It stops at v = 3,
1.8451 eV, so the helium half of a level-resolved cascade is no more buildable
than the hydrogen half, and there is still no helium THIRD-BODY rate
coefficient for the association itself (section 1.4).

### 2.3 The Cloudy c25.00 H2 data directory, read and measured

`~/CLOUDY/c25.00/data/h2/` (read-only; nothing in it was edited). Every
collision file has the same shape: a magic number, a `#>>refer` line naming the
published source, one line of temperatures, then one record per transition,
`v_upper J_upper v_lower J_lower` followed by the DE-EXCITATION rate
coefficients in cm^3 s^-1. MEASURED against this code's 302-level ladder:

| file | collider, source | transitions | levels, of 302 | highest level | max E [eV] | v = 11-14? | T grid |
|---|---|---|---|---|---|---|---|
| `coll_rates_H_07.dat` | H, Wrathmall, Gusdorf & Flower (2007) | 1485 | 78 | (0,23) | 2.9330 | no, v <= 3 | 60 pts, 100-6000 K |
| `coll_rates_H_15.dat` | H, Lique (2015) | 1431 | 54 | (3,8) | 1.8881 | no, v <= 3 | 50 pts, 100-5000 K |
| `coll_rates_H_99.dat` | H, Le Bourlot, Pineau des Forets & Flower (1999) | 627 | 51 | (3,8) | 1.8881 | no, v <= 3 | 9 pts, 100-6000 K |
| `coll_rates_He_LeBourlot.dat` | **He**, Le Bourlot et al. (1999) | 627 | 51 | (3,8) | 1.8881 | no, v <= 3 | 9 pts, 100-6000 K |
| **`coll_rates_He_ORNL.dat`** | **He, "Lee, T. G. et al. 2007, ApJ, in preparation"** | **3251** | **300** | (12,10) | **4.4764** | **YES, to v = 14** | 41 pts, 1-10000 K |
| `coll_rates_H2ortho_LeBourlot.dat` | **H2 (ortho perturber)**, Le Bourlot et al. (1999) | 361 | 39 | (1,13) | 1.6355 | no, v <= 2 | 9 pts, 100-6000 K |
| `coll_rates_H2para_LeBourlot.dat` | **H2 (para perturber)**, Le Bourlot et al. (1999) | 627 | 51 | (3,8) | 1.8881 | no, v <= 3 | 9 pts, 100-6000 K |
| `coll_rates_H2ortho_ORNL.dat` | H2 (ortho), Wan et al. (2018) | 240 | 32 | (0,31) | 4.3643 | no, **v = 0 only** | 64 pts, 1-10000 K |
| `coll_rates_H2para_ORNL.dat` | H2 (para), Wan et al. (2018) | 240 | 32 | (0,31) | 4.3643 | no, **v = 0 only** | 64 pts, 1-10000 K |
| `coll_rates_Hp.dat` | H+, Gerlich (1990) | 44 | 10 | (0,9) | 0.6202 | no, v = 0 only | 6 pts, 3-6000 K |

Every level in every file is a level of this code's ladder; none introduces a
level the ladder does not have.

**The find: `coll_rates_He_ORNL.dat` covers the whole ladder.** 300 of the 302
levels, every one except (0,31) and (14,4); levels per v run 31, 31, 29, 28,
26, 24, 23, 21, 19, 17, 15, 13, 11, 8, 4 for v = 0 to 14; **all 36 levels with
v >= 11 carry at least one downward transition**, and fourteen of the fifteen
levels the tree's line list leaves radiatively unconnected have downward
collisional rows here. **Collisionally, a level-resolved cascade for the HELIUM
collider is buildable from this file.** Its stated source, however, is "Lee,
T. G. et al. 2007, ApJ, in preparation": no such paper is in ADS, so **no
number is taken from it in this memo** and it is request item 2.

**And the H2-H2 collider is NOT missing after all, which corrects section 6.3.**
Le Bourlot et al. (1999) supply a state-resolved H2-H2 set, separately for the
ortho and the para perturber, to (3,8). Wan et al. (2018), the ORNL H2-H2
files, are pure rotational within v = 0 and carry no vibrational transition at
all, which matches their own abstract ("Rate coefficients for pure rotational
quenching in H2(v1 = 0, j1) + H2(v2 = 0, j2) collisions from initial levels of
j1 = 2-31"). Both papers are now in `references/`
(`LeBourlot_1999_MNRAS_305_802.pdf`, `Wan_2018_ApJ_862_132.pdf`).

### 2.4 The H2-H2 rate, and it settles the code's own open question

The code's `gamma_10_H2` carries Hollenbach & McKee's eq. (6.29) with 12000 in
the exponent, and `h2_vibrational_relaxation.f90` records that Burton et al.
(1990) reprint the same coefficient with 18100 there and that "Which is right
has NOT been checked against the original rate source". **The Le Bourlot data
settle it in favour of the Burton reading.** MEASURED, the thermal v = 1 to
v' = 0 rate from `coll_rates_H2*_LeBourlot.dat` (Boltzmann average over the J
of v = 1, summed into every J' of v = 0):

| T [K] | ortho-H2 perturber | para-H2 perturber | 3:1 mix | code fit, 12000 | Burton reading, 18100 |
|---|---|---|---|---|---|
| 808 | 1.03e-15 | 6.93e-16 | 9.47e-16 | 1.011e-13, **x107** | 4.85e-15, x5.1 |
| 1000 | 3.80e-15 | 2.63e-15 | 3.51e-15 | 1.893e-13, **x54** | 1.18e-14, x3.4 |
| 1527 | 3.42e-14 | 2.50e-14 | 3.19e-14 | 6.714e-13, **x21** | 7.17e-14, x2.3 |

**So the coefficient the code uses is 21 to 107 times too large over this
layer, and the reading it did not take is 2.3 to 5.1 times too large.** This is
reported and not applied: it is a source edit, outside this item.

### 2.5 What that does to the f_quench bound, and it changes the ARGUMENT

**The conservative bound of section 0 was built on the H2 channel alone and
that bound does not survive.** MEASURED, n_cr(all-level) against an H2 collider
alone moves from 5.53e+07 to **5.91e+09** cm^-3 at 808 K and from 8.33e+06 to
**1.75e+08** at 1527 K. At the thin end of the layer, n(H2) = 4e9 cm^-3 at
1527 K with **no other collider**, that gives f_quench = **0.958**, not 1.

**The conclusion survives, on the full collider mix, because atomic hydrogen
carries it.** MEASURED at cell 1 (T = 808.3 K, x2 = 0.3553, He/H = 2.13,
n(H2) = 5e13 so n_H = 2.81e14), with the best rate for each collider (Lique
2015 for H, Le Bourlot et al. 1999 for H2, and the Jozwiak or Lee value for
He):

| channel | share of the collider sum |
|---|---|
| atomic H | **0.9946** |
| He | 5.3e-03 |
| H2 | 5.9e-05 |

giving q = 8.0e+02 s^-1 against A_tot,max = 5.594e-06, i.e. **1 - f_quench =
7.0e-09**. And f_quench > 1 - 1e-3 at 1527 K needs only n(HI) > 3.0e+08 cm^-3
(MEASURED). **So f_quench = 1 in this layer is not in doubt; what changes is
that the claim must rest on the atomic-hydrogen channel, and the H2-only
statement of section 0 is retired.**

Note on the counting: Wrathmall et al. say 318 bound levels where this code's
ladder has 302. MEASURED: the code's ladder and the Meudon PDR 302-level file
agree to 0.053 K level by level, so the difference is a different bound-level
criterion near dissociation and not a different molecule.

---

### 2.6 Chasing the ORNL H2-He file to a published source: it does not reach one

Three leads were supplied by the user 2026-09-16 and read in the published
version. **Verdict: `coll_rates_He_ORNL.dat` is NOT attributable to a published
source, and where a published set overlaps it, it is superseded.**

| lead | what it actually is | does it publish the file's contents? |
|---|---|---|
| **Lee, Balakrishnan, Forrey, Stancil, Shaw, Schultz & Ferland (2008) ApJ 689, 1105**, `references/Lee_2008_ApJ_689_1105.pdf` | H2-**H2** rotational quenching, 2 to 10000 K, coupled-channel on a published (H2)2 surface. Its own title says so | **No.** It cites the H2-He set only in the text, as "T.-G. Lee et al. (2008, in preparation)" (READ, verbatim, its section 4). **That reference does not appear in its reference list at all**, so it cannot be chased, and ADS has no later paper matching it |
| **Mack, Clark, Forrey, Balakrishnan, Lee & Stancil (2006) PRA 74, 052718**, `references/Mack_2006PA_74_052718.pdf` | cold He + H2 near dissociation, numerically exact close coupling on the **Muchnick & Russek (MR) surface** (J. Chem. Phys. 100, 4336, 1994) | **No.** Cross sections at translational energies below 1000 cm^-1 and zero-temperature rate coefficients, all in figures (their Figs. 4 to 16). No thermal rate table, no temperature grid |
| **Ohlinger, Forrey, Lee & Stancil (2007) PRA 76, 042712**, `references/Ohlinger_2007PRA_76_042712.pdf` | H2 dissociation by He, quantum coupled-states with an L^2 Sturmian basis, same **MR surface** | **No.** Its tables are resonance parameters (energies and widths, Tables I and II); everything else is dissociation cross sections in figures. **It DOES reach the near-dissociation levels** ("the highest bound vibrational level for a given initial rotational level ranging from 0 to 20", READ, and its text names v = 14, j = 0-3 and v = 13, j = 1-5), but as cross sections over translational energy, not as the thermal rates a cascade needs |

So the file's lineage is identified, and it is the Forrey/Stancil group on the
Muchnick & Russek surface with Mack (2006) and Ohlinger (2007) as the published
fragments; **the thermal state-to-state table itself was never published.**

**AND IT IS SUPERSEDED WHERE IT CAN BE CHECKED.** MEASURED, over the 387
transitions `coll_rates_He_ORNL.dat` shares with the published, distributed
Jozwiak et al. (2024) CDS tables:

| quantity | at 300 K | at 1000 K |
|---|---|---|
| median ORNL / Jozwiak | **3.73** | **4.24** |
| 16th to 84th percentile | 1.20 to 11.2 | 1.47 to 13.9 |
| full range | 0.61 to 65 | 0.84 to 68 |

and channel by channel: (0,2) to (0,0) agrees to 1 and 9 per cent, (0,10) to
(0,8) to 30 per cent, while (1,0) to (0,0) is 8.5 and 8.2 times high, (1,2) to
(0,0) 9.4 times, (3,0) to (2,0) 5.8 and 5.2 times and (2,0) to (0,0) 45 and 27
times. **That is exactly the pattern Jozwiak et al. announce** (READ): good
agreement "for pure rotational transitions between low-lying rotational levels"
and "significant discrepancies for rovibrational processes involving
highly-excited rotational and vibrational states", which they attribute to the
broader range of intramolecular distances and the accuracy of their surface
against the older ones.

The file itself is internally consistent: MEASURED, every one of its 3251 rows
runs downward in energy, so it carries de-excitation only, as its format states.

**Consequence for the cascade.** The helium collisional half is no longer
"available pending a citation": the only set that reaches v = 11 to 14 rests on
a superseded surface, has no published thermal table, and is high by a median
factor 4 on exactly the rovibrational channels a cascade would use. **What is
needed is a modern H2-He set extended above v = 3**, which does not exist;
Jozwiak et al. (2024) is the modern set and it stops at v = 3.

## 3. Radiative decay: the A values, and what the tree already holds

**The published source is complete.** Wolniewicz, Simbotin & Dalgarno (1998),
ApJS 115, 293, abstract (READ, verbatim): "Accurate calculations of the
quadrupole moment of H2, carried out with a 494 term variational
representation of the electronic eigenfunction, are reported, and the
quadrupole transition probabilities connecting all the bound rovibrational
levels of H2 are presented." Obtained, `references/Wolniewicz_1998_ApJS_115_293.pdf`.

**What the tree holds is the later equivalent and it is a connected network,
with one hole.** MEASURED on `molecular_infrared_data.f90` (the electric
quadrupole plus magnetic dipole list of Roueff et al. 2019, A&A 630, A58,
table 2, as the generated file states):

- 302 levels spanning 0 to 51965.8 K (4.4781 eV), i.e. the whole bound X
  ladder to dissociation; 1833 lines.
- Every line's upper AND lower level is in the ladder: matching each line's
  `Tu - dE` against the level list leaves a largest residual of 0.099 K and no
  unmatched line. The list is therefore usable as a transition network without
  any external level set.
- 285 of the 302 levels have at least one downward transition.
- The 17 that do not are (0,0) and (0,1), which cannot decay, **and 15 levels
  near dissociation: all of v = 13 (J = 0, 1, 2, 4, 6) and all of v = 14
  (J = 0 to 4), plus (10,0), (11,0), (11,2), (12,0) and (12,2)** (v, J labels
  MEASURED by matching the code's ladder to the PDR level file).
- Largest total spontaneous decay rate over all levels, A_tot,max =
  5.594e-06 s^-1; the v = 1, J = 0 total is 8.532e-07 s^-1, which reproduces
  Hollenbach & McKee's adopted 8.3e-07 s^-1 to 3 per cent. The largest values
  are high-J rotational levels of low v, not the near-dissociation levels: the
  the maxima by vibrational level run 5.51e-06 (v = 0), 5.59e-06 (v = 1), 4.90e-06 (v = 10),
  1.01e-06 (v = 13) s^-1.
- Cross-referenced against Bossion et al. (2018), MEASURED: five of those
  fifteen ((10,0), (11,0), (11,2), (12,0), (12,2)) do have collisional rates in
  that calculation, and the other ten, every level of v = 13 and v = 14, have
  neither a collisional rate nor a downward radiative transition in anything
  now available.
- Cross-referenced against the two collisional sets whose numbers are now in
  hand, MEASURED: Lique (2015) covers 54 of the 302 levels to 1.8881 eV and
  Jozwiak et al. (2024) 53 of them to 1.8451 eV, so **neither reaches any of
  the fifteen**, which all lie above 3.84 eV. The radiative holes and the
  collisional holes are the same region of the ladder.
- Independent confirmation that Wolniewicz, Simbotin & Dalgarno (1998) is the
  right repair for the fifteen: Nesterenok et al. (2019) build their own H2
  network on exactly that source (READ, verbatim): "We take into account 298
  rotational levels of the ground electronic state of H2 molecule for which the
  Einstein coefficients are given by Wolniewicz, Simbotin & Dalgarno (1998).
  The level energies of H2 are taken from Dabrowski (1984)." Their 298 against
  this code's 302 is a difference of bound-level criterion near dissociation.
- **AND THE REPAIR IS NOW MEASURED, not only argued.** Cloudy c25.00 carries
  the Wolniewicz set as `~/CLOUDY/c25.00/data/h2/transprob_X.dat`, which its
  own header attributes to "Wolniewicz, L., Simbotin, I., and Dalgarno, A.,
  1998, ApJS, 115, 293-313" (READ) and which gives `EU VU JU EL VL JL A`.
  MEASURED on it: **4661 X-to-X lines over 301 levels, 299 of them with a
  downward transition** ((0,0) and (0,1) cannot decay). **Fourteen of the
  tree's fifteen unconnected levels are connected here**, with total downward
  rates A_tot = 5.00e-06 (10,0), 4.58e-06 (11,0), 4.60e-06 (11,2), 3.85e-06
  (12,0), 3.83e-06 (12,2), 2.76e-06 (13,0), 2.74e-06 (13,1), 2.69e-06 (13,2),
  2.46e-06 (13,4), 1.94e-06 (13,6), 1.28e-06 (14,0), 1.23e-06 (14,1),
  1.13e-06 (14,2) and 9.59e-07 s^-1 (14,3). **The one that stays unconnected is
  (14,4)**, the topmost level of this code's ladder, which Wolniewicz does not
  treat at all; that is also the level absent from `coll_rates_He_ORNL.dat` and
  it explains Nesterenok et al.'s 298 against this code's 302.
- **The all-level maximum is confirmed from outside the tree.** MEASURED on
  `transprob_X.dat`, max A_tot = **5.5945e-06 s^-1 at (1,28)**, against the
  5.594e-06 measured from the tree's own Roueff-based line list: the two
  independent sets agree to four digits, so the bound the memo carries is not
  an artefact of one line list.
- The (v, J = 0) origins, MEASURED from the ladder, for use anywhere a
  vibrational energy is needed: v = 1 0.5159, v = 2 1.0027, v = 3 1.4608,
  v = 4 1.8908, v = 5 2.2927, v = 6 2.6664, v = 7 3.0113, v = 8 3.3266,
  v = 9 3.6109, v = 10 3.8622, v = 11 4.0775, v = 12 4.2529, v = 13 4.3830,
  v = 14 4.4601 eV.

**There is no cascade solver.** `molecular_infrared_cooling.f90` uses these
tables in LTE: it forms `f_u(T) = g_u exp(-T_u/T)/Q(T)` with the partition
function of `caloric_eos`, sums the optically thin line emission and the
absorption of the diluted base field, and tabulates the result against
temperature. The levels are populated by a Boltzmann factor, not by a rate
network, so nothing in the code can currently carry a nascent distribution.

---

## 4. What thermalizes the prompt translational share

A nascent H2 leaves the association with some kinetic energy, and the code
deposits that share promptly. That is legitimate only if it thermalizes faster
than transport and faster than the reactions that consume the molecule.

MEASURED, with the base state READ from the predecessor memo section 6
(cell 1 at 808.3 K, base at 1 microbar and 3e13 cm^-3, n(H2) = 4e9 to
5e13 cm^-3 across the layer) and the planet parameters READ from
`LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH2.13/input.inp`
(R_p = 1.1274e+09 cm, M_p = 3.3449e+28 g, so g = 1756 cm s^-2):

| quantity | value |
|---|---|
| speed of an H2 carrying 0.086 eV (kT at 1000 K) | 2.87e+05 cm s^-1 |
| speed of an H2 carrying 0.5 eV | 6.92e+05 cm s^-1 |
| elastic collision time at the base, n = 3e13 cm^-3, sigma = 1e-15 cm^2 | 1.2e-04 s (0.086 eV) to 4.8e-05 s (0.5 eV) |
| the same with a deliberately small sigma = 1e-16 cm^2 | 1.2e-03 to 4.8e-04 s |
| elastic collision time at the top of the layer, n = 4e9 cm^-3 | 0.9 s to 0.4 s (sigma = 1e-15), 9 s to 4 s (sigma = 1e-16) |
| pressure scale height at the base, mu = 2 to 4 | 2.7e+06 to 5.3e+06 cm at 226 K, 9.5e+06 to 1.9e+07 cm at 808 K |
| transport time over a scale height at a speed of 1e4 cm s^-1 (forty times the largest base velocity the P44 note records) | 2.5e+02 to 5.0e+02 s |
| chemical turnover time of H2 at cell 1, n(H2)/R15 rate, with the rate 6.92e+04 cm^-3 s^-1 READ from the predecessor memo | 2.9e+08 to 7.2e+08 s |

**The margin is at least five orders at the base and at least two orders at
the top of the layer, against the shortest competing time, and twelve orders
against the chemistry.** Prompt deposition of the translational share is safe
everywhere in this layer, by any reading of the elastic cross section within
an order of magnitude of gas kinetic. Row R40's second clause is therefore
answered: the assumption holds here, measured, and the statement belongs at
the code site with these numbers.

---

## 5. R5, R6, R7 and R16: who receives the energy

The ledger currently deposits all of q for each of these four (READ,
`molecular_reaction_heat.f90` lines 486 to 500; only R15 carries a fraction
and only R23 has its photon removed).

| # | reaction | q [eV] (READ, predecessor memo) | what the published record says about the products | status |
|---|---|---|---|---|
| R5 | H2+ + e -> H + H | +10.948 | **SETTLED 2026-09-16, and it does NOT all become heat.** Giusti-Suzor, Bardsley & Derkits (1983) PRA 28, 682 treat recombination through the lowest doubly excited state of H2, the 1Sigma_g+ state dominated by the 1sigma_u^2 configuration, whose potential curve is constructed so that "the final segment is used to ensure dissociation to the limit H(1s)+H(2s)" (READ, verbatim, their section IV A). Takagi (2002) Phys. Scr. T96, 52 states it as the general result (READ, verbatim, his section 3): "One of the hydrogen atoms of the dissociative product is in the ground electronic state, but another is in the excited state, whose principle quantum number is n = 2: H2+(v) + e -> H(1s) + H(n = 2). The product of n = 2 is dominantly produced only for the low vibrational molecular ion and at low collision energies. Except for this condition, the halves of product atoms are distributed to the highly excited states." **So at least 10.1988 eV of the 10.948 eV is electronic excitation of one hydrogen atom, leaving 0.749 eV as kinetic energy: 93.2 per cent of q(R5) is misdirected in the present ledger** (MEASURED from the two energies) | **SETTLED in direction; the branching to n >= 3 is not, and the code's H2+ is outside Takagi's low-v condition** (the code's own site note records that the H2+ lifetime against R5/R8/R9, ~1e-3 s, is far shorter than its radiative vibrational relaxation, ~1 s). Giusti-Suzor et al. state their own range: electron energies below 0.5 eV, H2+ in the lowest three vibrational states, the lowest doubly excited state only, and "For electrons with energy greater than 2 eV, or for ions in more highly excited vibrational states, the contributions of other doubly excited states become more important" |
| R6 | H3+ + e -> H2 + H | +9.250 | **the H2 fragment is born vibrationally hot, now read in BOTH primary papers.** Kokoouline, Greene & Esry (2001) Nature 412, 891 (READ, verbatim): "this direct-pathways calculation results in a peak vibrational distribution for the channels H2(v, j) that peaks at v ~ 5-6, which is in accord with recent experimental results", alongside their three-body branching "sigma(H + H + H)/sigma_DR = 0.70 +/- 0.07" which "overlaps with the experimental result of 0.75 +/- 0.08". Strasser et al. (2001) PRL 86, 779 measured it in the TSR (READ, verbatim): "The resulting vibrational distribution [Fig. 1(b)] is found to be wide with a peak around v = 5", and "The vibrational state distribution obtained here provides the first quantitative evidence about product excitation in two-body DR of H3+." MEASURED on the code's ladder, v = 5 is 2.2927 eV and v = 6 is 2.6664 eV, so **the PEAK of the distribution is 25 to 29 per cent of q(R6)**, internal and subject to the same radiate-or-thermalize branching as R15 | **the PEAK is settled, the MEAN is not.** Both papers state the distribution is broad and neither tabulates it; Strasser et al. give the spread as "Primarily, Ek varies by ~4.5 eV depending on the vibrational state of the H2 fragment", the full 0 to 4.478 eV range of the ladder, and quote "an uncertainty of about one level" on the fitted distribution |
| R7 | H3+ + e -> 3H | +4.771 | three atoms, all kinetic. The branching between R6 and R7 is confirmed from outside the code: the same review gives the three-body branching ratio as 0.70 +/- 0.07, and the code's rates 5.04e-08 and 2.16e-08 (300/Te)^0.65 (READ, `mol_rates.f90`) give exactly 0.70 and 0.30 (MEASURED) | SETTLED, no internal share |
| R16 | HeH+ + e -> He + H | +11.753 | **SETTLED 2026-09-16 on the published paper.** Guberman (1994) PRA 49, R4277, abstract: "The dominant dissociative route is the C 2Sigma+ state leading to H(n = 2) atoms", and in the text (READ, verbatim): "Calculated cross sections along the X and A states are three to four orders of magnitude smaller than C state cross sections. Therefore the dissociation products will nearly always include an excited n = 2 H atom, in agreement with the recent TARN II storage ring results." **So 10.1988 eV of the 11.753 eV is electronic excitation, leaving 1.554 eV as kinetic energy: 86.8 per cent of q(R16) is misdirected in the present ledger** (MEASURED). The code already carries an H(n = 2) population and a Lyman-alpha channel it could be delivered to | **SETTLED in direction.** Ranges READ at the source: the calculation is for 3HeH+ and "the calculated cross sections are sensitive to the isotopomer under study", and it covers the ground-state ion at electron energies 0.001 to 0.33 eV |

R15's own entry is section 1. R23 was corrected by the predecessor item.

**Note on priority.** R16 is the largest share in one event but HeH+ is a trace
carrier here; R6 is a smaller share of a channel that carries 49 per cent of
the H2 loss at cell 1. Neither moves this layer, because f_quench = 1 makes
internal and prompt energy the same thing here. Both matter for a shallower
base.

---

## 6. The model the data support, and its bounds

### 6.1 The form

```
  Q_R15 = k3b(T) n_M n_HI^2  *  D0 * [ f_trans + (1 - f_trans) * f_quench ]
```

with D0 = 4.478 eV, and the same two-parameter treatment applied to the
internal share of R6 (and, if the request-list papers confirm it, to R16's
electronic share, which is a different recipient and not this form).

### 6.2 The bounds, and where each comes from

| parameter | bound | source of the bound |
|---|---|---|
| `f_trans`, the prompt translational share of D0 | **0 <= f_trans <= 0.09, with the indicated value below 0.01** (revised 2026-09-16 from gas-phase three-body calculations; the earlier 0 to 0.05, transferred from grain formation, is retired) | Orel (1987) for M = H: the nascent distribution "is peaked near this line, that is, in the highest bound states of the H2 manifold", with v = 11 to 14 dominant and the top J of each v carrying it; Paolini, Ohlinger & Forrey (2011) for M = He: "recombination to either of the two most weakly bound vibrational levels for each rotational level j <= 20"; Esposito & Capitelli (2009) for M = H: "at low temperature only high lying w states are significant". MEASURED on the code's ladder, the highest bound level of v = 11 to 14 leaves D0 - E_int = 0.0000 to 0.0180 eV, 0.00 to 0.40 per cent of D0; the 0.09 is the extreme bound obtained by putting the whole population at the J = 0 origin of v = 11, which no source supports. **The code's f_trans = 0 is data-supported, not merely conservative** (section 1.6) |
| `f_quench`, the collisionally thermalized share of the internal energy | in THIS layer, 1 - 1.4e-02 <= f_quench <= 1 - 1.1e-06 | MEASURED: `(1 + n_cr/n)^-1` with n_cr = A_tot,max/gamma and A_tot,max = 5.594e-06 s^-1 over all 302 levels, gamma the Hollenbach & McKee v = 1 H2 coefficient, n(H2) = 4e9 to 5e13 cm^-3 |
| direction of the remaining approximations in f_quench | net effect small, and now measured rather than argued | (a) the v = 1 collisional coefficient understates de-excitation of the closely spaced high levels, which RAISES f_quench; (b) helium is left out of the collider sum, which LOWERS f_quench, but by 7.5e-03 of the H2 channel at 1000 K even at n(He) = 2.6 times n(HI) + n(H2) (MEASURED from Balakrishnan et al. 1999 Table 5, section 2.2); (c) the code's `gamma_10_H` is 1.05 to 3.47 times Lique's (2015) accurate value over this layer's 200 to 1600 K, MEASURED from the supplementary data file across the whole range (section 2.1), which LOWERS f_quench when corrected, moving n_cr against atomic H from 6.78e+05 to 1.27e+06 cm^-3 at 808 K and from 2.76e+05 to 3.01e+05 at 1527 K, still three to four orders below the layer |

**Neither parameter is to be adjusted to recover a base temperature.** With
f_quench = 1 in this layer the two are degenerate and the base temperature is
insensitive to both: the largest movement the whole family can produce here is
f_trans times the 1.4e-02 that f_quench falls short at the thinnest point,
i.e. below 1e-03 of the R15 heat, against a base heating that exceeds the base
cooling by a factor 48 (READ, predecessor memo section 4). **The 808 K base is
not explained by the association heat fraction.**

### 6.3 What would change the recommendation

**Revised 2026-09-16 on the published Bossion, Scribano, Lique & Parlant
(2018).** Reading it moves the obstruction without removing it: the
calculation of the right scope exists for M = H, and what is missing is the
numbers, two other colliders and a corner of the physics.

**Does a level-resolved cascade become buildable for M = H, now that Orel
(1987) is in hand? Not yet, and the gap is smaller and sharper.** The
qualitative initial condition is settled (the molecule is born in the top
levels, section 1.6) and that is what fixes f_trans. A cascade needs the
initial condition as NUMBERS, and none of the three state-resolving papers
tabulates them: Orel writes "The full results are available on request",
Esposito & Capitelli present theirs in figures, and Paolini et al. likewise.
Bossion et al. (2018) would give them by detailed balance over 260 of this
code's own 302 levels, and distributes nothing either. What still blocks it,
in order:

1. **The numbers are not published, and a fourth paper that USES them does
   not distribute them either.** Bossion et al. (2018) carries figures and one
   room-temperature table with no supporting information, repository or
   data-availability statement. **Nesterenok, Bossion, Scribano & Lique (2019)
   MNRAS 489, 4520 was checked for exactly this** (supplied by the user
   2026-09-16, `references/Nesterenok_2019MNRAS_489_4520.pdf`, read in the
   published version) and it carries **no data-availability statement, no
   supporting information, no repository link, no fit formula and no table of
   rate coefficients**: its Table 1 lists only which reference each collider
   is taken from, and its Appendix A describes the shock calculation. **So
   request item 1 does not close.** What that paper DOES establish is that the
   needed extension EXISTS, in its own words (READ, verbatim, section 2.1):
   "For this study, Bossion et al. (2018) calculations have been extended to
   all the bound states of H2 molecule (internal energy up to dissociation
   limit ~55 100 K) and up to 20 000 K of collisional energy for the rate
   constant ... For this high-temperature regime, we extended the QCT
   cross-sections up to 100 000 K of collisional energy in order to ensure
   convergence on the rate constants. We used an energy step of 2000 K up to
   that limit." **That set covers the whole bound ladder, v = 13 and v = 14
   included, with the three-body collisional dissociation in it**, which is
   precisely the initial condition and the collisional half a cascade needs
   for M = H. It exists, it is held by those four authors, and it is
   unpublished. Item 1 of the request list now names them. They carry the
   quasi-bound caveat unchanged (READ, verbatim): "we considered the
   quasi-bound states as pertaining to the three-body dissociation channel.
   This assertion remains valid as long as the average lifetime of the
   quasi-bound molecules is lower than the typical collision time; this is
   true for low- to moderate-density media" -- and this layer is not one.
   (MEASURED, no LAMDA- or BASECOL-format collisional file for H2-H exists
   anywhere under `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/`.)
2. **The inversion is distorted exactly where it matters.** Quasi-bound H2 is
   counted into the dissociative channel, up to 50 per cent of it at 5000 K
   (READ), and the paper itself states the assumption fails in high-density
   media. This layer is one. Using the inversion as a nascent distribution
   would therefore need the quasi-bound share put back, which the published
   figures do not resolve.
3. **M = H2: the PRODUCT states are untouched, but the COLLISIONAL rates are
   not, which corrects an earlier statement in this memo.** Schwenke (1988), the
   only three-body calculation for M = H2, gives the TOTAL rate and no
   product-state table (section 1.6). But a state-resolved H2-H2 COLLISIONAL
   set does exist and is in hand: Le Bourlot, Pineau des Forets & Flower (1999),
   carried by Cloudy as `coll_rates_H2ortho_LeBourlot.dat` and
   `coll_rates_H2para_LeBourlot.dat`, 51 levels to (3,8) (section 2.3). It
   stops at v = 3 like every other set whose numbers are distributed, and it
   shows the code's own `gamma_10_H2` to be 21 to 107 times too large over this
   layer (section 2.4).
4. **M = He is untouched by Bossion et al., and the best helium set, now in
   the tree as numbers, stops at v = 3.** Jozwiak et al. (2024) supersedes
   Balakrishnan et al. (1999), covers 53 of the code's 302 levels to 1.8451 eV
   over 20 to 8000 K, and its CDS tables are read and measured in section
   2.2b; v = 11 to 14 are absent there. What it does settle numerically is that
   helium is 0.1 per cent of the collider sum the code forms at the base
   (MEASURED), because that sum is 99.5 per cent atomic hydrogen.
   **REVISED TWICE 2026-09-16.** Cloudy's `coll_rates_He_ORNL.dat` covers 300
   of this code's 302 levels, reaches v = 14 and gives every one of the 36
   levels with v >= 11 at least one downward transition, over 1 to 10000 K
   (MEASURED, section 2.3), so at first sight the helium collisional half
   looked available. **Chasing its source closed that door** (section 2.6): the
   thermal table was never published, and against the published Jozwiak et al.
   (2024) set it is high by a median factor 3.7 to 4.2 on the rovibrational
   channels while agreeing to within 30 per cent on the pure rotational ones.
   **So for helium too, nothing usable reaches above v = 3.** For the helium THIRD BODY of the
   association itself, Paolini, Ohlinger & Forrey (2011) resolves the product
   states, putting the nascent molecule in "the two most weakly bound
   vibrational levels for each rotational level j <= 20" (section 1.6), but it
   tabulates no numbers either and stops at 1000 K.
5. **The radiative holes are CLOSED except one, measured 2026-09-16.**
   `~/CLOUDY/c25.00/data/h2/transprob_X.dat`, the Wolniewicz, Simbotin &
   Dalgarno (1998) set, connects fourteen of the fifteen; only (14,4), the
   topmost level, stays unconnected, and Wolniewicz does not treat it at all
   (section 3). Its all-level maximum A_tot, 5.5945e-06 s^-1 at (1,28), agrees
   with the tree's own list to four digits. The paragraph below is the state
   BEFORE that measurement and is kept for the record.
   *(superseded)* **The radiative holes remain, and ten of them are now holes
   twice over.**
   Of the 15 near-dissociation levels with no downward transition in the
   code's line list, Bossion et al. cover five collisionally; the other ten,
   all of v = 13 and v = 14, have neither a collisional rate nor a radiative
   decay in anything now available (MEASURED). The radiative side is fixable
   from Wolniewicz, Simbotin & Dalgarno (1998), which is now in `references/`
   and gives the quadrupole probabilities connecting ALL bound levels; the
   collisional side of those ten is not. **Nesterenok et al. (2019) confirm
   that this is the right repair**: they build their own network on exactly
   that source (READ, verbatim) -- "We take into account 298 rotational levels
   of the ground electronic state of H2 molecule for which the Einstein
   coefficients are given by Wolniewicz, Simbotin & Dalgarno (1998). The level
   energies of H2 are taken from Dabrowski (1984)" -- 298 levels against this
   code's 302, the difference being the bound-level criterion near
   dissociation.

6. **The two sets that DO distribute their numbers stop at v = 4.** Lique
   (2015) reaches 1.891 eV over 55 of the code's levels and Balakrishnan et
   al. (1999) 1.105 eV over 20 of them (MEASURED). They are the accurate data
   for the bottom of the ladder and they settle the two collisional
   corrections of section 6.2, but they say nothing about the levels a nascent
   molecule occupies.

7. **The nascent level is also the most easily re-dissociated, so the rate
   and the heat are coupled.** Orel: "Due to the highly excited nature of the
   products, the nascent distribution will be greatly changed by subsequent
   collisions"; Schwenke: "the primary bottleneck occurs in the relaxation of
   bound vibrational states close to the dissociation limit"; Bossion et al.
   book the quasi-bound population into the dissociative channel for the same
   reason. A cascade that resolves the deposited heat has to resolve the
   effective three-body rate at the same time, which is a larger change than
   the heat term alone.

**So the recommendation of section 6.1 and 6.2 stands unchanged in FORM, and
one of its two parameters is now pinned by data.** f_trans is below 0.01 from
three gas-phase calculations that resolve the product states (section 1.6), so
the two-parameter family collapses to the one-parameter form the code already
computes, and the remaining question is only the level resolution of f_quench.

**The cascade is still not buildable, and after this round the obstruction is
a single one.** Of the four ingredients:

| ingredient | state 2026-09-16 |
|---|---|
| the radiative network over the whole ladder | **SOLVED, measured 2026-09-16**: Cloudy's `transprob_X.dat` (Wolniewicz et al. 1998) connects 14 of the 15 holes; only (14,4), which Wolniewicz does not treat, stays out |
| H2-H collisional rates to v = 3 | **in hand as numbers**, `Rates_H_H2.dat`, 1431 transitions, 100-5000 K (also Cloudy's `coll_rates_H_07.dat`, 78 levels to 2.933 eV) |
| H2-He collisional rates to v = 3 | **in the tree as numbers**, the Jozwiak et al. (2024) CDS tables, 1059 transitions, 20-8000 K |
| H2-He collisional rates over the WHOLE ladder | **NOT USABLE, chased 2026-09-16.** Cloudy's `coll_rates_He_ORNL.dat` does cover 300 of 302 levels to v = 14, but its thermal table was never published, and where the published Jozwiak et al. (2024) set overlaps it, it is high by a median factor 3.7 to 4.2 on the rovibrational channels (MEASURED, section 2.6). **What is needed is a modern H2-He set above v = 3, which does not exist** |
| H2-H2 collisional rates to v = 3 | **in hand as numbers**, Cloudy's Le Bourlot files, 51 levels to (3,8) |
| H2-H collisional rates ABOVE v = 3, and the nascent distribution by detailed balance from the state-resolved three-body dissociation | **THE ONE THING STILL MISSING FOR M = H.** The set exists, covers all bound states to the dissociation limit, and is held unpublished by Nesterenok, Bossion, Scribano & Lique |

For M = H2, the third body the code's own rate coefficient belongs to, no
state-resolved product calculation exists at all (Schwenke 1988 gives the
total rate), and nothing changes that.

### 6.4 Uncertainty of the recommendation

- The f_quench = 1 conclusion for this layer is robust: the all-level bound is
  1.9 orders below the thinnest n(H2) in the layer and 5.9 orders below the
  base value, and the two measured corrections to the collider sum (helium
  omitted, `gamma_10_H` 1.05 to 3.47 times too large over 200 to 1600 K) move
  n_cr against atomic H by less than a factor 2 at the base (MEASURED,
  sections 2.1 and 2.2).
- **f_trans is no longer the weak part** (revised 2026-09-16). Three
  gas-phase three-body calculations that resolve the product states, two for
  M = H (Orel 1987, Esposito & Capitelli 2009) and one for M = He (Paolini,
  Ohlinger & Forrey 2011), all put the nascent molecule in the highest bound
  levels, where D0 - E_int is 0.0000 to 0.0180 eV on this code's own ladder
  (MEASURED). What none of them supplies is the tabulated distribution, so
  f_trans is pinned as a NUMBER below 0.01 while the initial condition of a
  cascade is still not in hand. For M = H2, the third body the code's own rate
  coefficient belongs to, Schwenke (1988) gives the total rate only, so the
  product states of that channel rest on the M = H and M = He results by
  analogy.
- The R6 internal share, 25 to 29 per cent of 9.250 eV, is now read in both
  primary papers rather than in a review, but it is the position of the PEAK
  and not a mean: both describe the distribution as broad, neither tabulates
  it, and Strasser et al. quote about one vibrational level of uncertainty.
  The mean could be lower or higher.
- **R5 and R16 are new and unbooked.** The published product states put 93.2
  per cent of q(R5) and 86.8 per cent of q(R16) into electronic excitation of
  a hydrogen atom rather than into fragment kinetic energy (MEASURED from
  10.1988 eV against 10.948 and 11.753 eV). Both are outside the present
  ledger, both are inside their sources' stated domains only in part (the
  code's H2+ is vibrationally hot, outside Takagi's low-v condition; Guberman
  computed 3HeH+ and calls the cross sections isotopomer sensitive), and the
  bound on what correcting them would do here is set by R15 carrying 95.1 per
  cent of `heat_mol_chem` at cell 1 (READ, 4.95e-07 of 5.205e-07 erg cm^-3
  s^-1), so **every other channel together is at most 4.9 per cent of the base
  heating** (MEASURED).
- Nothing here has been run. Every temperature and density is READ from the
  predecessor memo's measurements on `.L7f/g_newbc` and `.L7f/heatpp_noir`.

---

## 7. Source table

Obtained means the published version is in `references/` and was read there.

| source | bibcode | what it gives | validity | obtained |
|---|---|---|---|---|
| Hollenbach & McKee (1979) ApJS 41, 555 | 1979ApJS...41..555H | the `(1 + n_cr/n)^-1` form; the 4.2 + 0.2 eV split for GRAIN formation; the v = 1 and v = 2 collisional fits eq. (6.29) | fits stated good to 20% above ~500 K, 50% above ~300 K; NO gas-phase three-body product distribution | YES |
| Cohen & Westberg (1983) JPCRD 12, 531 | 1983JPCRD..12..531C | k1(H2), k1(H), k1(Ar) for H + H + M | 50-5000 K (H2, H), 77-5000 K (Ar); log k1(H2) uncertain to 0.2 at 300 K rising to 0.4 at 5000 K; no He, no product states | YES |
| Flower & Harris (2007) MNRAS 377, 705 | 2007MNRAS.377..705F | k3 = 1.44e-26 T^-1.54 by detailed balance; level-resolved populations over 49 levels to 20000 K | primordial collapse conditions; three-body enters as a total only | YES, added |
| Glover & Jappsen (2007) ApJ 666, 1 | 2007ApJ...666....1G | the "essentially all into internal excitation" assumption, verbatim | H-, H2+ and grain routes; three-body formation heating not included | YES |
| Glover & Abel (2008) MNRAS 388, 1627 | 2008MNRAS.388.1627G | the three-body rate uncertainty and its thermal effect; heating referred to Glover & Jappsen | primordial gas | YES, added |
| Palla, Salpeter & Stahler (1983) ApJ 271, 632 | 1983ApJ...271..632P | the three-body rate and the k(H2) = k(H)/8 ratio | primordial collapse | YES, added |
| Lepp & Shull (1983) ApJ 270, 578 | 1983ApJ...270..578L | level-resolved collisional dissociation kinetics | the dissociation side only | YES, added |
| Wolniewicz, Simbotin & Dalgarno (1998) ApJS 115, 293 | 1998ApJS..115..293W | quadrupole transition probabilities connecting ALL bound rovibrational levels of H2 | complete for the X state | YES, added |
| Wrathmall, Gusdorf & Flower (2007) MNRAS 382, 133 | 2007MNRAS.382..133W | H2-H rovibrational rates, lowest 108 of 318 bound levels | to 6000 K; stops at 2.819 eV of internal energy | YES, added |
| **Lique (2015) MNRAS 453, 810** | 2015MNRAS.453..810L | nearly exact quantum H2-H (de-)excitation rates with the hydrogen exchange channels included; the thermal room-temperature v = 1 to v' = 0 relaxation, 1.8e-13 cm^3 s^-1, against which the code's `gamma_10_H` is 3.4 times too large (MEASURED) | 100-5000 K; internal energies below 22000 K, MEASURED as 55 of this code's 302 levels to 1.891 eV, v = 0 to 4; no three-body dissociation | YES, supplied by the user 2026-09-16; the rate coefficients are distributed as supplementary data to the article (not in this workspace) |
| **Balakrishnan, Forrey & Dalgarno (1999) ApJ 514, 520** | 1999ApJ...514..520B | close-coupled quantum H2-He rotational and ro-vibrational rates, tabulated in the paper; the only helium collider data in hand | 10-5000 K (rotational de-excitation), fits valid 100-5000 K and 100-3000 K; v = 0 and v = 1 with j = 0 to 9 only, MEASURED as 20 of this code's 302 levels to 1.105 eV | YES, supplied by the user 2026-09-16 |
| **Jozwiak, Thibault, Viel, Wcislo & Lique (2024) A&A 685, A113** | 2024A&A...685A.113J | the current H2-He state-to-state rates, 1059 transitions, superseding Flower et al. (1998) and Balakrishnan et al. (1999); the numbers are distributed at the CDS (J/A+A/685/A113) and BASECOL | 20-8000 K on a 43-point grid; internal energies below 15000 cm^-1, MEASURED as 53 of this code's 302 levels to 1.8451 eV, v = 0 to 3, v = 11 to 14 absent; the tables obey detailed balance to four digits (MEASURED) | YES, **paper AND CDS tables**, both supplied by the user 2026-09-16: `references/Jozwiak_2024A&A_685_A113.pdf` and `references/Jozwiak_2024J_A+A_685_A113/` (`ReadMe`, `oh2-lev.dat`, `ph2-lev.dat`, `oh2-rat.dat`, `ph2-rat.dat`) |
| **Nesterenok, Bossion, Scribano & Lique (2019) MNRAS 489, 4520** | 2019MNRAS.489.4520N | the record that the Bossion et al. (2018) calculation was extended "to all the bound states of H2 molecule (internal energy up to dissociation limit ~55 100 K)" with the three-body collisional dissociation in it; and that its own network uses Wolniewicz et al. (1998) A values over 298 levels | C-type shock modelling; **it distributes NO rate coefficients**: no data-availability statement, no supporting information, no repository, no fit formula, no table of rates (READ, whole text). The quasi-bound caveat is carried unchanged and is stated to hold only "for low- to moderate-density media" | YES, supplied by the user 2026-09-16 |
| **Bossion, Scribano, Lique & Parlant (2018) MNRAS 480, 3718** | 2018MNRAS.480.3718B | QCT H2-H inelastic, reactive and **state-resolved three-body collisional dissociation** rate coefficients over 260 initial (v, j), which MEASURED are 260 of this code's own 302 ladder levels, to 49931.7 K = 4.3028 eV | 100-15000 K translational, on 20 collision energies from 0.1 eV upward; mean deviation against the quantum calculation "of the order of a factor 3" above a few hundred K; quasi-bound H2 booked into the dissociative channel (up to 50 per cent of it at 5000 K), an assumption the paper states fails in high-density media; v = 13 and v = 14 absent | **YES, the paper**, supplied by the user 2026-09-16, `references/Bossion_2018_MNRAS_480_3718.pdf`. **The rate coefficients themselves are NOT distributed with it** |
| Larsson, McCall & Orel (2008) CPL 462, 145 | 2008CPL...462..145L | the H3+ DR three-body branching 0.70 +/- 0.07 and the H2(v) peak at v = 5-6 | storage-ring energies, rotationally cold ions | YES |
| Yelle (2004) Icarus 170, 167 | 2004Icar..170..167Y | the state of practice: chemical heating from heats of formation, no product partition | hot Jupiter thermosphere | YES |
| **Orel (1987) JCP 87, 314** | 1987JChPh..87..314O | the nascent (v, J) distribution of `H + H + H -> H2(v,J) + H`: peaked in the highest bound states, v = 11 to 14 dominant, top J of each v carrying it | six orbiting resonances; cross sections at E/k = 50 to 2000 K; the state-resolved rate coefficients are plotted to 350 K and the numbers are "available on request", not tabulated | YES, supplied by the user 2026-09-16 |
| **Esposito & Capitelli (2009) JPCA 113, 15307** | 2009JPCA..11315307E | the same for M = H over a much wider temperature range, with the rotational distribution resolved; total rate 1.6e-32 cm^6 s^-1 at 300 K reproducing Orel's 1.19e-32 | 300 to 10000 K; distributions in figures, not tabulated | YES, supplied by the user 2026-09-16 |
| **Schwenke (1988) JCP 89, 2076** | 1988JChPh..89.2076S | the TOTAL three-body rate for M = H2, the third body the code's own coefficient belongs to, from an ab initio H4 surface and four-body trajectories | 100-5000 K; "the theoretical results about a factor of 2 too small over the temperature range 300-5000 K" against experiment; NO product-state table anywhere | YES, supplied by the user 2026-09-16 |
| **Paolini, Ohlinger & Forrey (2011) PRA 83, 042713** | 2011PhRvA..83d2713P | the only helium third-body calculation in hand, state-resolved: recombination goes to "either of the two most weakly bound vibrational levels for each rotational level j <= 20" | figures to 1000 K; compared with experiment at 77 and 300 K; He and Ar only | YES, supplied by the user 2026-09-16 |
| **Kokoouline, Greene & Esry (2001) Nature 412, 891** | 2001Natur.412..891K | the H3+ DR three-body branching 0.70 +/- 0.07 and the H2(v, j) product distribution peaking at v ~ 5-6 | direct 2p pathways only; the indirect np pathways are outside the calculation | YES, supplied by the user 2026-09-16 |
| **Strasser et al. (2001) PRL 86, 779** | 2001PhRvL..86..779S | the measured H2(v) distribution after two-body H3+ DR, "wide with a peak around v = 5"; the first quantitative product-excitation evidence for this channel | storage ring, vibrationally cold H3+ at a rotational temperature of 0.23(3) eV; about one vibrational level of uncertainty; the distribution is not tabulated | YES, supplied by the user 2026-09-16 |
| **Guberman (1994) PRA 49, R4277** | 1994PhRvA..49.4277G | the HeH+ DR route: the C 2Sigma+ state, "the dissociation products will nearly always include an excited n = 2 H atom" | 3HeH+, ground-state ion, electron energies 0.001 to 0.33 eV; cross sections stated to be isotopomer sensitive | YES, supplied by the user 2026-09-16 |
| **Giusti-Suzor, Bardsley & Derkits (1983) PRA 28, 682** | 1983PhRvA..28..682G | H2+ DR through the lowest doubly excited state, whose curve is built "to ensure dissociation to the limit H(1s)+H(2s)"; the rate for ground-state ions is smaller than for v = 1 and 2 by a factor of about 6 | electron energies below 0.5 eV, H2+ in the lowest three vibrational states, one doubly excited state only | YES, supplied by the user 2026-09-16 |
| **Takagi (2002) Phys. Scr. T96, 52** | 2002PhST...96...52T | the H2+ DR product electronic state by initial vibrational state: "One of the hydrogen atoms of the dissociative product is in the ground electronic state, but another is in the excited state, whose principle quantum number is n = 2", holding "only for the low vibrational molecular ion and at low collision energies" | MQDT extended to dissociative states; H2+, HD+ and D2+ over v = 0 to 18 | YES, supplied by the user 2026-09-16 |
| **Le Bourlot, Pineau des Forets & Flower (1999) MNRAS 305, 802** | 1999MNRAS.305..802L | the state-resolved H, He and H2 collisional sets Cloudy carries; **the only state-resolved H2-H2 vibrational set in hand**, and the one that shows the code's `gamma_10_H2` to be 21 to 107 times too large | "rovibrational transitions between all the energy levels up to approximately 20 000 K above the ground state" (READ); MEASURED as 51 levels to (3,8) = 1.8881 eV, 9 points 100-6000 K | YES, obtained 2026-09-16, `references/LeBourlot_1999_MNRAS_305_802.pdf` |
| **Wan, Yang, Stancil et al. (2018) ApJ 862, 132** | 2018ApJ...862..132W | the ORNL H2-H2 files; **pure rotational within v = 0 only**, so no vibrational transition and no bearing on f_quench | "pure rotational quenching in H2(v1 = 0, j1) + H2(v2 = 0, j2) collisions from initial levels of j1 = 2-31" (READ); MEASURED as 32 levels, v = 0 only, 64 points 1-10000 K | YES, obtained 2026-09-16, `references/Wan_2018_ApJ_862_132.pdf` |
| "Lee, T. G. et al. 2007/2008, ApJ, in preparation" | none in ADS | the source Cloudy names for `coll_rates_He_ORNL.dat`, the only collisional set in hand that covers the whole bound ladder (300 of 302 levels, to v = 14, 1-10000 K, MEASURED) | **NEVER PUBLISHED**, chased 2026-09-16 through the three leads of section 2.6; superseded where checkable, median factor 3.7 to 4.2 high against Jozwiak et al. (2024) on the rovibrational channels | **NO, and not worth chasing further.** No number is taken from the file |
| **Lee, Balakrishnan, Forrey, Stancil, Shaw, Schultz & Ferland (2008) ApJ 689, 1105** | 2008ApJ...689.1105L | H2-**H2** rotational quenching, 2-10000 K, coupled channel; the only place the H2-He set is mentioned, and only as "T.-G. Lee et al. (2008, in preparation)" in the text, absent from its reference list | rotational transitions only; not the H2-He set | YES, supplied by the user 2026-09-16 |
| **Mack, Clark, Forrey, Balakrishnan, Lee & Stancil (2006) PRA 74, 052718** | 2006PhRvA..74e2718M | cold He + H2 near dissociation, close coupling on the Muchnick & Russek (1994) surface; cross sections and zero-temperature rate coefficients, in figures | translational energies below 1000 cm^-1; **no thermal rate table** | YES, supplied by the user 2026-09-16 |
| **Ohlinger, Forrey, Lee & Stancil (2007) PRA 76, 042712** | 2007PhRvA..76d2712O | H2 dissociation by He, coupled states with an L^2 Sturmian basis, same Muchnick & Russek surface; **reaches the near-dissociation levels** (the highest bound v for each j = 0 to 20, including v = 14, j = 0-3 and v = 13, j = 1-5) | resonance parameters in tables, dissociation cross sections in figures; **no thermal rate table** | YES, supplied by the user 2026-09-16 |
| Roueff et al. (2019) A&A 630, A58 | 2019A&A...630A..58R | the H2 line list and ladder the code carries (302 levels, 1833 lines) | full X state; 15 near-dissociation levels carry no downward line in what the code holds | the DATA are in the tree; the paper itself is not in `references/` |

---

## 8. Request list

Papers that could not be obtained in published form, plus one set of numbers
that is not published at all. None of them is cited for a value anywhere
above. Retrieval was attempted through the ADS link gateway (PUB_PDF,
ADS_PDF, ARTICLE) and, for the OUP titles, directly; the publishers refuse
anonymous access.

**Updated 2026-09-16, second revision. Every paper of the original list has
now been supplied by the user and read in the published version**, and all of
them have moved to the source table of section 7: Bossion et al. (2018),
Orel (1987), Schwenke (1988), Esposito & Capitelli (2009), Paolini, Ohlinger &
Forrey (2011), Lique (2015), Balakrishnan, Forrey & Dalgarno (1999),
Kokoouline, Greene & Esry (2001), Strasser et al. (2001), Guberman (1994),
Giusti-Suzor, Bardsley & Derkits (1983) and Takagi (2002).

**An attribution correction carried through the whole memo:** Phys. Rev. A 28,
682 (1983) is Giusti-Suzor, Bardsley & Derkits, "Dissociative recombination in
low-energy e-H2+ collisions", not Guberman; the earlier attribution here was
wrong and is fixed everywhere.

**Third revision, 2026-09-16.** Three further items were supplied and read:
the Lique (2015) supplementary data file, which **closes what was item 2**;
Jozwiak et al. (2024) A&A 685, A113, which supersedes Balakrishnan et al.
(1999) for the helium collider and distributes its numbers at the CDS; and
Nesterenok, Bossion, Scribano & Lique (2019) MNRAS 489, 4520, checked for
whether it distributes the Bossion et al. rates. **It does not**, so item 1
stays open, but it establishes that the needed extension to all bound states
exists and names a fourth author who holds it.

**Fourth and fifth revisions, 2026-09-16.** The Cloudy c25.00 H2 data
directory was read and measured (section 2.3): it closes the radiative holes,
supplies the H2-H2 collisional set the memo had reported missing (Le Bourlot et
al. 1999, now in `references/`), and contains one file covering the whole bound
ladder for helium whose source is unpublished. **That source was then chased
through three published leads and does not reach one** (section 2.6), so it is
struck from the list rather than added to it: the file is superseded where a
published set overlaps it.

**What remains is ONE set of numbers that no paper carries, plus one
calculation that has not been done.**

| # | reference | bibcode / DOI | what must be checked in it | what is blocked without it |
|---|---|---|---|---|
| 1 | **The RATE TABLES of Bossion et al. (2018) AS EXTENDED FOR Nesterenok et al. (2019)**, i.e. the H2-H state-to-state (de-)excitation and three-body collisional dissociation rate coefficients over ALL bound states of H2 to the dissociation limit (~55100 K internal energy) and to 20000 K of gas temperature. Request from the four authors who hold it: A. V. Nesterenok (Ioffe), D. Bossion and Y. Scribano (Montpellier), F. Lique (Rennes/Le Havre). Neither Bossion et al. (2018) nor Nesterenok et al. (2019) distributes it (READ, both texts checked in full: no data-availability statement, supporting information, repository, fit formula or table of rates) | 2018MNRAS.480.3718B, 10.1093/mnras/sty2089; 2019MNRAS.489.4520N, 10.1093/mnras/stz2441 | the state-to-state `k(v,j -> v',j'; T)` table and, separately, the state-resolved dissociation rate `k_diss(v,j; T)`; and whether the quasi-bound contribution can be separated out of the dissociative channel, since both papers book it there and both state the assumption holds only for low- to moderate-density media | **the ONE remaining obstruction to a level-resolved cascade for M = H.** Its detailed-balance inverse is the nascent distribution as NUMBERS, and it is the only set that covers the levels the nascent molecule occupies. Everything else a cascade needs is now either in the tree or distributed (section 6.3) |
| 2 | **A modern H2-He set extended above v = 3.** This is a calculation that has not been done, not a paper to obtain. Jozwiak et al. (2024) is the current H2-He set and stops at v = 3 (1.8451 eV); the only set that reaches v = 14 is Cloudy's `coll_rates_He_ORNL.dat`, whose thermal table was never published and which is high by a median factor 3.7 to 4.2 against Jozwiak on the rovibrational channels (section 2.6). The natural request is to the Jozwiak or Lique group: is an extension above v = 3 planned or available? | Jozwiak et al.: 2024A&A...685A.113J | whether the modern surface has been or can be carried to the near-dissociation levels | the helium half of a level-resolved cascade. Note this is a WEAKER need than item 1: helium is 0.1 per cent of the collider sum at the base (MEASURED, section 2.2b) |
| *(closed 2026-09-16)* | the Lique (2015) supplementary file (now `references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat`, section 2.1); the Jozwiak et al. (2024) CDS tables (now `references/Jozwiak_2024J_A+A_685_A113/`, section 2.2b); and the hunt for a published source of the ORNL H2-He file, which is answered in the negative and needs no further chasing (section 2.6) | | | |

---

**Status 2026-09-17 (the user's decision).** The extended Bossion et al. (2018)
H2-H tables cannot be obtained now. The level-resolved cascade of item L7g
is therefore DEFERRED, not planned: the code keeps the bounded one-parameter
form `q(R15) f_quench` with the all-level critical density and the
collider-resolved rates of the model item, whose validity in this layer is
measured (`1 - f_quench` about 1e-6 at the base). The request stays on
record here; if the tables arrive, the cascade is the item to reopen, with
the radiative half already in the tree (Wolniewicz et al. 1998) and the
helium half still without a modern set above v = 3.

## 9. What this item changed

Nothing in `src/`. **Updated 2026-09-16**: the user supplied TWELVE published
papers, the whole of the original request list, filed in `references/` as
`Bossion_2018_MNRAS_480_3718.pdf`, `Lique_2015_MNRAS_453_810.pdf`,
`Balakrishnan_1999_ApJ_514_520.pdf`, `Orel_1987JCP_87_314.pdf`,
`Schwenke_1988JCP_89_2076.pdf`, `Esposito_2009_JPCA_113_15307.pdf` (renamed
from `Esposito_2009JCP_113_15307.pdf` after the PDF header confirmed
J. Phys. Chem. A), `Paolini_2011PRA_83_042713.pdf`,
`Kokoouline_2001Nature_412_891.pdf`, `Strasser_2001PhysRevLett_86_779.pdf`,
`Guberman_1994PRA_49_R4277.pdf`, `Giusti-Suzor_1983PRA_28_682.pdf` and
`Takagi_2002_PhST_96_52.pdf`. All are read in sections 1.6, 2 and 5 and
entered in the source table of section 7; sections 0, 1.3, 1.6, 2, 2.1, 2.2,
3, 5, 6.2, 6.3, 7 and 8 were revised on them. **The model recommendation of
section 6.1 is unchanged in form, and its f_trans parameter is now pinned by
data below 0.01 instead of bounded by transfer from grain formation.** An
attribution error was corrected: PRA 28, 682 (1983) is Giusti-Suzor, Bardsley
& Derkits, not Guberman.

**Third update 2026-09-16**: the user supplied the Lique (2015) supplementary
data file (`references/Lique_2015MNRAS_453_810_data/Rates_H_H2.dat`, which
closes request item 2), `references/Jozwiak_2024A&A_685_A113.pdf` and
`references/Nesterenok_2019MNRAS_489_4520.pdf`. Sections 0, 2, 2.1, 2.2a (new),
6.2, 6.3, 6.4, 7 and 8 were revised on them. The correction factor on the
code's `gamma_10_H` is now measured across 200 to 1600 K rather than at one
point (1.05 to 3.47); the helium collider has a current, distributed set
(Jozwiak et al. 2024, CDS J/A+A/685/A113); and **the request list is down to
one item**, the Bossion rates as extended to all bound states for Nesterenok
et al. (2019), which neither paper distributes. The user then also supplied the
Jozwiak CDS tables themselves (`references/Jozwiak_2024J_A+A_685_A113/`), read
and measured in the new section 2.2b: they cover 53 of the code's 302 levels to
1.8451 eV over 20-8000 K, **do not reach v = 11 to 14**, obey detailed balance
to four digits, revise Balakrishnan et al. (1999) downward by a factor 8.3 on
the (1,0) to (0,0) channel at 1000 K, and put helium at **0.106 per cent of the
collider sum the code forms at the base** (which is 99.7 per cent atomic
hydrogen). The single-channel estimate of section 2.2 was about an order of
magnitude low and is corrected there.

**Fourth and fifth updates 2026-09-16**: the Cloudy c25.00 H2 data directory
(`~/CLOUDY/c25.00/data/h2/`, read-only, unmodified) was read and measured, and
five further published papers were obtained or supplied for `references/`
(`LeBourlot_1999_MNRAS_305_802.pdf`, `Wan_2018_ApJ_862_132.pdf`,
`Lee_2008_ApJ_689_1105.pdf`, `Mack_2006PA_74_052718.pdf`,
`Ohlinger_2007PRA_76_042712.pdf`). New sections 2.3, 2.4, 2.5 and 2.6;
sections 0, 3, 6.3, 7 and 8 revised. Four things changed: the **radiative holes
are closed** (14 of 15, by Wolniewicz via `transprob_X.dat`, whose all-level
A_tot maximum reproduces the tree's to four digits); a **state-resolved H2-H2
set does exist** (Le Bourlot et al. 1999) and it shows the code's
`gamma_10_H2` to be **21 to 107 times too large**, settling the code's own open
12000-against-18100 question in favour of the Burton reading; the **H2-only
form of the f_quench bound is retired**, because with the measured H2-H2 rate a
pure-H2 collider would give f_quench = 0.958 at the thin end, so the claim now
rests on the atomic-hydrogen channel (99.5 per cent of the collider sum at the
base); and the **ORNL H2-He file, which alone reaches v = 14, is neither
published nor reliable** where a published set overlaps it. **The model
recommendation is unchanged throughout.** Six published PDFs were
added to `references/` by this item itself:
`Flower_2007_MNRAS_377_705.pdf`, `Glover_2008_MNRAS_388_1627.pdf`,
`Lepp_1983_ApJ_270_578.pdf`, `Palla_1983_ApJ_271_632.pdf`,
`Wolniewicz_1998_ApJS_115_293.pdf`, `Wrathmall_2007_MNRAS_382_133.pdf`.
Row L7g of `docs/PLAN_20260913_lhs_stationary.md` records the inventory.

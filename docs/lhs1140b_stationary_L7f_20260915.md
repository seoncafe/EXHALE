# LHS 1140 b, item L7f: which chemistry is right about the base layer

Item L7f of `docs/PLAN_20260913_lhs_stationary.md`. The item opened on a
finding of item L7e: that inside the base layer of
`molecular_scalar_gj1132_kzz1e9/HeH2.13` the lower-atmosphere handoff
(`q_H2_base` = 0.19011, i.e. x2 = 2 n(H2)/n_H = 0.99998) stands **7.24 to
21.00 times** the root of the wind network's own H2 carrier row, so that one
of the two chemistries had to be wrong about a layer both describe.

Nothing in `build/`, `EXHALE.x`, `backup/regression/` or any catalogue case
directory of `LHS1140b/models/` was written by this item: every binary is a
private build in `build_L7f/`, every run is under `LHS1140b/models/.L7f/`,
and the regression matrix was run on a copy of `backup/regression/` in the
scratch directory.

| build | what it is | md5 |
|---|---|---|
| `EXHALE_L7f_diag.x` | the tree as delivered plus the reaction-by-reaction record of the H2 row (section 2) and nothing else | `63266ad0eaf5d345c9d181dcc7bb3ff8` |
| `EXHALE_L7f.x`, the SEED build | the record, the solved chemical root and the seed's root-imbalance report; the control solve and the control regression of section 5 were run with it | `d153215676a0bb9697a2b5a50ad6dafd` |
| `EXHALE_L7f.x`, DELIVERED | the above plus the three rate repairs of section 3 (R20 retired, the Black channel, R17 rewritten) | `48db9e876031aa3cf39198ccdf933d24` |
| `EXHALE_L18_ctl.x` | the tree as delivered to this item, used only to attribute the `mol_metals` golden failure | `922dec0fa183...` (first 12, as the harness prints it) |
| `EXHALE.x` = `EXHALE_L7f.x`, the CAMPAIGN build | the delivered text plus item L21's base-boundary correction, which landed in the tree while this item ran; the tree binary and this item's are the same file | `db87b88d1ce53facf1d61084fa535ca5` |

The regression of section 5.3 was run on `48db9e876031`, i.e. BEFORE L21's
boundary correction, so its golden movement is this item's rate repairs alone
and is not contaminated by that change.

Every number below is MEASURED -- produced by a run of one of those builds,
or by arithmetic stated here on a file such a run wrote -- unless it is
marked READ (from a source file, a published paper or a document).

**SUPERSEDED IN ONE PART, 2026-09-15, and the pointer is here so the memo is
not read as current where it is not.** Every base-face and base-temperature
measurement below was made on the base boundary condition AS IT STOOD. Item
L21 then found that the boundary was choosing its branch on the SIGN OF CELL
1's velocity -- the odd-even artefact of the P44 note -- and was therefore
taking 97 percent of the face state from the interior isentrope instead of
from the reservoir. With that corrected, the same case on the same seed and
the same network gives a base ghost at **226.02 K** (the `base.inp` value, in
place of 873 K), a first cell at 808 K in place of 1023, a cell-1 mass row of
**4.92e-09** in place of 1.68e-07, and the state is **ACCEPTED at outer pass
11** where it was refused after forty. So section 6's proposal (a) is no
longer the open question it was: the composition step across the base face is
still there, and the row certifies anyway. The numbers are in
`docs/lhs1140b_stationary_L7g_20260915.md` section 6. **Sections 1 to 4 --
the root correction, the reaction-by-reaction record and the three rate
repairs -- are independent of the boundary and stand as written.**

## 1. Verdict

**Neither chemistry is wrong about the base layer. The 7 to 21 was an
arithmetic artefact of how the "root" was formed, and with the root actually
solved for, the wind network and the lower-atmosphere handoff agree at the
bottom of the column to 1.7 percent in x2.**

`carrier_h2_chemical_root` formed the root of the H2 carrier row as
production divided by loss RATE at the trial composition, on the stated
ground that "the loss of H2 is proportional to n(H2) and its production is
not". **That ground is false.** The dominant formation channel is the
three-body association R15, `k15 n(H I)^2`, and `carrier_source` closes
atomic hydrogen out of the element budget,

    n(H I) = n_H,available - 2 n(H2) - n(OH) - 2 n(H2O) - [n(H+)],

so the production FALLS as n(H2) rises and is exactly zero where every
hydrogen nucleus is already bound. The row is quadratic in its unknown, and
P/(L/n) at the trial is not its root. The seed evaluates it immediately
after applying the thermochemical fit q_H2(p, T), and that fit is SATURATED
at x2 = 1 over the whole base layer -- so the quotient was evaluated at
n(H I) ~ 0, where the production it divides has collapsed.

MEASURED, on the atomic state the case is seeded from
(`.L7d/src_atomic_HeH2.13`, base at 417.6 K and 1 microbar):

| | x2 at cell 1 | ratio handoff/root, cells 1-142 |
|---|---|---|
| P/(L/n), the quotient | 8.915e-02 | **7.24 to 21.00** |
| the row's actual root | **9.837e-01** | **1.02 to 5.87** |
| the handoff | 9.9998e-01 | - |

So at the base the network's own H2 balance sits 1.7 percent below the
composition the photochemical column hands over, and the layer they disagree
about is not the base at all but the top of the layer, where the
temperature has risen past 1200 K and the ratio reaches 5.9.

**What destroys H2 in that layer is ion chemistry, not heat and not
photons.** The row was written out reaction by reaction for the first time
(section 2). Thermal dissociation R12 stands **nine decades** below the
first cell's total H2 loss at 1023 K and **twenty-two decades** below it at
629 K, and nowhere in the layer does it reach 0.4 percent of that loss;
electron-impact dissociation R14 is further below still; Lyman-Werner
photodissociation, correctly attenuated
(tau from a star-ward column of 5.7e18 cm^-2, shielding factor 7.4e-4), is
2 to 4 percent; H2 photoionization, which survives to the base only as
X-rays, is 4 to 17 percent. The leading terms are the helium-ion channels
(R17/R20/R23, 20 to 39 percent) with R18 (HeH+ + H2) behind them, and the
proton channels R10/R13 (8 to 53 percent) -- and much of the latter recycles,
because the H3+ that R13 makes returns H2 through R6.

**The formation routes the item suspected are not needed.** The L7f row of
the plan named the H^- route (H + e -> H^- + hv, H^- + H -> H2 + e) and
grain-surface formation as absent from the network and asked whether their
absence accounts for a factor 7 to 21. It does not, because there is no
factor 7 to 21: the network's own three-body association already puts the
base at x2 = 0.984 against the handoff's 0.99998, and a further formation
channel would only push a root that is already at 98 percent of the element
ceiling closer to it. Neither route is adopted, and the reason is now a
measurement rather than a judgement.

**AND ONE ENTRY OF THE RATE TABLE IS WRONG, which is the item's second
finding.** R20, He+ + H2 -> HeH+ + H at 4.2e-13, carried 25 to 30 percent of
the H2 loss through the base of this helium-rich case. The paper it cites
was obtained on 2026-09-15 (`references/Schauer_1989JCP_91_4593.pdf`) and it
does not support that channel: it measures the radiative and the dissociative
charge transfer over 15 < T < 40 K, states that the HeH+ channel "apparently
does not become allowed until the collision energy approaches 9 eV", and
reports an H3+ signal that is the SUM of the radiative and the HeH+ channels
-- which bounds the HeH+ channel at <= 1.0e-14, **42 times below the value
Table 1 gives it**. R20 is retired, and the HeH+ source the network now
carries is the one Koskinen Table 1 omits, H2+ + He -> HeH+ + H at
3.0e-10 exp(-6717/T), read from the primary source (Black 1978). Section 3.2
has the argument, the measurement of what the repair does -- it rewires the
base layer's HeH+ and H+ budgets by decades and leaves the net H2 destruction
there within a percent, because the helium ions are sink-limited -- and the
one entry it leaves open, R17.

## 2. The H2 row, reaction by reaction

The code carried no record of the individual reactions of a carrier row: the
row dump item L7e added (`EXHALE_CARRIER_ROW_TERMS=1` ->
`output/carrier_row_terms.txt`) gives the chemical production and the
chemical loss as two numbers. Each is a sum of several reactions, and which
of them carries the row is the physics question this item asks, so the dump
now writes a second line under every H2 row:

```
  chan <cell> R15_3body R9_H2p_H R6_H3p_e R11_H3p_H P_H2_photo LW_photdis
              R10R13_Hp R12_thermal R14_edis R8_H2p_H2 R17R20R23_Hep
              R18_HeHp HeI23S_H2 oxy_prod oxy_loss
```

each a volumetric rate in cm^-3 s^-1. The channels are handed out by
`mol_heh_rows` at the one place the terms are written and are assembled from
the same named factors as the row itself, so the row is unchanged to the last
bit and the record cannot drift from it. Diagnostic only; nothing in the
solution reads it.

`docs/figures/lhs1140b_L7f_base_h2.pdf` carries both halves of the
measurement: on the left the base-layer x2 of the solved state against the
corrected seed, the handoff and the value the old quotient read at the first
cell, and on the right every H2 loss channel as a fraction of the loss
through the layer. The step at the base face and the eleven-decade gap
between the thermal dissociation and the helium-ion channels are what the
figure is for.

### 2.1 The solved state, `.L7e/fixed5/HeH2.13` as written

The state whose cell-1 mass row refuses the certification. Its base layer is
at 1023 K (cell 1) to 1710 K (cell 142) at 1 microbar, **not** at the 417.6 K
of the atomic state the case is seeded from: the same pressure level is 2.4
times hotter in the molecular solution than in the atomic one.

| | cell 1 | cell 25 | cell 100 | cell 142 |
|---|---|---|---|---|
| r [R_p] | 1.00019 | 1.00483 | 1.02496 | 1.05333 |
| T [K] | 1023.0 | 1544.5 | 1544.2 | 1710.3 |
| n(H2) [cm^-3] | 3.787e+11 | 1.976e+11 | 7.237e+10 | 1.707e+10 |
| **production** | **6.958e+04** | 8.734e+03 | 1.297e+03 | 9.224e+02 |
| R15 three-body | 99.43 % | 88.64 % | 33.04 % | 0.89 % |
| R9 H2+ + H | 0.14 % | 8.57 % | 53.42 % | 83.18 % |
| R6 H3+ + e | 0.43 % | 2.63 % | 13.41 % | 15.80 % |
| **loss** | **1.916e+03** | 2.693e+03 | 2.516e+03 | 2.467e+03 |
| loss rate [s^-1] | 5.06e-09 | 1.36e-08 | 3.48e-08 | 1.44e-07 |
| P_H2 photoionization | 6.74 % | 3.49 % | 4.64 % | 4.23 % |
| LW photodissociation | 2.66 % | 1.98 % | 1.77 % | 1.60 % |
| R10 + R13, H+ | 46.82 % | 52.56 % | 45.71 % | 46.59 % |
| **R12 thermal** | **1.4e-06 cm^-3 s^-1** | 0.37 % | 0.05 % | 0.09 % |
| R14 electron impact | 3.7e-12 cm^-3 s^-1 | 1.1e-04 | 9.5e-05 | 7.7e-04 |
| R8 H2+ + H2 | 4.07 % | 21.47 % | 20.32 % | 19.47 % |
| R17 + R20 + R23, He+ | 38.59 % | 20.04 % | 27.38 % | 27.93 % |
| R18 HeH+ + H2 | 1.11 % | 0.10 % | 0.13 % | 0.08 % |
| diffusive divergence | +5.696e+05 | +6.072e+03 | -1.132e+03 | -1.344e+03 |
| advective divergence | -5.020e+05 | -3.059e+01 | -8.728e+01 | -2.001e+02 |

Two things are visible at once.

**At cell 1 the chemistry is not what sets x2.** The transport terms stand an
order of magnitude above the net chemical source: the diffusive divergence
alone is 5.7e+05 against a production of 7.0e+04. The base face carries NO
diffusive flux by construction (the carrier operator's faces 0 and N are
zero-flux, a documented decision), so that 5.7e+05 is the 1-2 face: the layer
is eddy-mixed by `He_Kzz` = 1e9 with the dissociation region above it, and
x2 = 0.33 through cells 1 to 142 is a quenched, well-mixed value and not a
chemical root. The mixing time over the layer, (6.0e7 cm)^2/1e9 = 3.6e6 s,
against the H2 formation time n(H2)/P = 5.4e6 s at cell 1: the two compete,
which is why the quenched value sits between 0 and 1.

**Thermal dissociation is absent from the budget.** At 1023 K, R12 is
1.4e-06 cm^-3 s^-1 against a total loss of 1.9e+03 -- nine decades down --
and it does not reach a tenth of a percent anywhere in the layer. The H2 of
a helium-rich, 1 microbar, weakly ionized base is destroyed by He+, by HeH+,
by H+ and by X-rays, and the temperature enters only through the rate
coefficients of those channels.

### 2.2 The cold state, at the seed's own composition

The same record on the state the corrected seed writes, reloaded and measured
after one outer pass (`.L7f/dump_seed418`), whose base carries the seed's
x2 = 0.98 at 629 K:

| | cell 1 | cell 25 | cell 100 | cell 142 |
|---|---|---|---|---|
| T [K] | 628.9 | 739.9 | 1269.6 | 1828.7 |
| x(H2) = n(H2)/n_tot | 1.864e-01 | 1.714e-01 | 3.251e-02 | 4.179e-02 |
| **production** | 9.514e+02 | 1.347e+03 | 1.896e+02 | 7.262e+02 |
| R15 three-body | 30.75 % | 47.11 % | 34.86 % | 0.05 % |
| R6 H3+ + e | 68.43 % | 49.65 % | 9.38 % | 12.73 % |
| **loss** | 4.455e+03 | 4.549e+03 | 1.171e+03 | 1.994e+03 |
| P_H2 photoionization | 16.57 % | 16.01 % | 8.45 % | 4.14 % |
| LW photodissociation | 2.88 % | 2.05 % | 2.58 % | 4.46 % |
| R10 + R13, H+ | 8.02 % | 18.66 % | 3.06 % | 43.39 % |
| **R12 thermal** | **1.9e-19 cm^-3 s^-1** | 8.2e-15 | 1.2e-04 | 0.08 % |
| R8 H2+ + H2 | 16.56 % | 16.19 % | 3.41 % | 15.39 % |
| R17 + R20 + R23, He+ | **31.83 %** | **32.94 %** | **82.01 %** | 32.48 % |
| R18 HeH+ + H2 | **24.13 %** | 14.16 % | 0.49 % | 0.06 % |

At 629 K the thermal dissociation is **nineteen decades** below the leading
term. **The helium-ion channels and the HeH+ they make are 56 percent of the
H2 destruction at the first cell, and 82 percent at cell 100.** This is the
measurement that promotes the recorded He-dominated-limit caveat to a
load-bearing one.

## 3. Are the coefficients right at these temperatures?

Checked against their sources for the channels the record shows carrying the
row. The full transcription audit is
`docs/molecular_chemistry_audit_he_rich.md`; only what this item adds is here.

| channel | coefficient | validity at 400-1800 K | verdict |
|---|---|---|---|
| R15, H + H + M | 2.8e-31 T^-0.6 n_M, Cohen & Westberg (1983), recommended 50-5000 K (READ) | the layer is inside the range | **superseded 2026-09-17**: that is the M = H2 coefficient, and n_M was the TOTAL heavy-particle density; the rate is now the collider sum over H2, H and He with the same paper's per-collider coefficients (He given argon's as an upper bound), which lowers the rate 22 per cent at this base. Item L7g, `docs/lhs1140b_stationary_L7g_model_20260916.md` |
| R12, H2 + M | k(R15)/K_eq(T), detailed balance on the same channel | exact reverse of R15 by construction; irrelevant here at 1e-19 to 1e-06 cm^-3 s^-1 | sound and negligible |
| LW photodissociation | Draine & Bertoldi (1996) self-shielding on the star-ward H2 column | MEASURED at cell 1: N(H2) = 5.70e18 cm^-2, f_shield = 7.37e-04, k_LW = 1.35e-10 s^-1. The column is integrated from the star inward and is largest at the base (6.4e18 at the outer ghost falling to 3.9e17 at cell 142), which is the right sense; f_shield at that column agrees with the Draine & Bertoldi N^-0.75 branch | attenuation correct; 2 to 4 percent of the loss |
| P_H2, photoionization | the band integral of the run's own spectrum | 3.4e-13 s^-1 at cell 1 under an H I column of ~1e20 cm^-2, which is opaque above 15 eV; the surviving rate is the X-ray band, and 9.3 erg cm^-2 s^-1 at 1 keV times sigma(H2) ~ 2e-22 cm^2 gives ~1e-12 s^-1, the same order | physical, and it is X-rays |
| R13, H+ + H2 + M | 3.2e-29 n_M, Miller et al. (1968) | the loss it carries returns largely through R6, H3+ + e -> H2 + H, which is 9 to 68 percent of the production in the same cells; the pair is one cycle | recycling, not net destruction |
| R17, He+ + H2 -> H+ + H + He | 1e-9 exp(-5700/T), Moses & Bass (2000) | 1.2e-13 at 629 K, 3.8e-12 at 1023 K; the fit is an extrapolation below ~650 K and falls up to 2e4 under three direct measurements there | see 3.2, left as published and marked |
| **R20, He+ + H2 -> HeH+ + H** | **4.2e-13, Koskinen Table 1 citing Schauer et al. (1989)** | 77 percent of the He+ channel at 629 K and the largest single H2 sink in the layer | **RETIRED; the cited paper bounds it at <= 1.0e-14, see 3.2** |
| H2+ + He -> HeH+ + H (new) | 3.0e-10 exp(-6717/T), Black (1978) Eq. 12 and rate table, from the Chupka et al. (1969) cross sections | 4.4e-16 at 500 K, 1.0e-11 at 2000 K | adopted, default on; it is the network's only HeH+ source once R20 is gone |

### 3.1 The figure the LW and photoionization rows rest on

Both are read from the run's own output (`output/Lyman_Werner.txt` and the
row record), so they are measurements of this code and not of a source.

### 3.2 R20: the cited measurement does not support it

The paper was obtained (`references/Schauer_1989JCP_91_4593.pdf`, READ). It
is an ion-trap measurement of He+ + H2 at 15 to 40 K. Its Sec. I lists the
three energetically allowed channels -- (1) dissociative charge transfer,
He+ + H2 -> He + H + H+, dH = -6.51 eV; (2) radiative charge transfer,
-> He + H2+, dH = -9.16 eV; (3) -> H + HeH+, dH = -8.36 eV -- and then says
of the third:

> Reaction (3) apparently does not become allowed until the collision energy
> approaches 9 eV.

and of what the experiment can separate:

> In the experiments described here, we detect the formation of both H+ and
> H3+ . Both H2+ and HeH+ react rapidly with H2 to form H3+ , so observation
> of this product serves to monitor the sum of reactions (2) and (3).
> Because reaction (3) is apparently ruled out, as just discussed, we
> nominally attribute observation of H3+ to reaction (2), and refer to it as
> radiative charge transfer (RCT)

So the paper's H3+ channel is the SUM of the radiative charge transfer and
the HeH+ channel, and its Fig. 4 gives that sum as 7.5e-15 to 1.1e-14 over
15-40 K with no temperature dependence, plus an upper limit of 1.0e-14 over
70-320 K from Johnsen et al. (1980) drawn in the same figure (all READ from
the figure, whose axis is 1e-15 cm^3 s^-1). **The measurement therefore
bounds the HeH+ channel at <= 1.0e-14 and the network carried 4.2e-13, 42
times more.** The reason for the bound is structural rather than
statistical: the reactants and products of He+ + H2 -> HeH+ + H do not
correlate (their Fig. 1) and the channel opens near 9 eV of collision
energy, which is 1e5 K in relative kinetic energy -- decades above anything
this code integrates.

**What replaces it.** With R20 gone, Koskinen Table 1 has no HeH+ formation
channel at all, and a species with destruction and no formation is not a
model. The network now carries the route Table 1 omits and Garcia Munoz
(2025) has, traced to its own primary source: Black, J. H. 1978, ApJ 222,
125, "Molecules in planetary nebulae", Eq. (12) and the rate table,
H2+ + He -> HeH+ + H at 3.0e-10 exp(-6717/T) from the Chupka, Berkowitz &
Russell (1969) cross sections. Black states its validity himself -- "This
reaction is substantially endothermic for ground-state H2+; however, it is
rapid for H2+ with v > 3, states which may be well populated at T > 5000 K"
-- so the coefficient is a Maxwellian average over a THERMAL population of
H2+ vibrational levels, which is written at the code site together with what
that assumption is worth in a layer where the H2+ lifetime is shorter than
its vibrational relaxation. The reaction is endothermic by 0.805 eV on this
code's own formation-energy table, which is the published ground-state
endothermicity; its heat enters the one reaction-heat ledger by the same
table difference as every other channel, and `physics_probe` holds that
identity.

**What the repair does, measured** (first cell of
`molecular_scalar_gj1132_kzz1e9/HeH2.13`, one sweep on each side so the
ionization is re-solved, not frozen):

| | n(He+) | n(HeH+) | n(H2+) | n(H+) | total H2 loss [cm^-3 s^-1] |
|---|---|---|---|---|---|
| 1023 K, before | 4.617e+02 | 3.751e-02 | 1.030e-01 | 1.091e+07 | 1.916e+03 |
| 1023 K, after | 5.125e+02 | 7.704e-05 | 7.309e-02 | 1.129e+07 | 1.904e+03 |
| 629 K, before | 1.280e+03 | 3.512e-01 | 1.808e-01 | 5.000e+05 | 4.455e+03 |
| 629 K, after | 5.645e+03 | 4.364e-06 | 1.908e-01 | 1.935e+06 | 4.443e+03 |

**The share overstated the effect and the reason is physics.** R20 was 25 to
30 percent of the H2 loss at a frozen composition, but the helium ions of
that layer are SINK-limited: remove one He+ sink and n(He+) rises until the
others carry the same flux. So the net H2 destruction moves by less than a
percent, while n(HeH+) falls by 487 at 1023 K and by 8e4 at 629 K and n(H+)
rises by 3.9 at 629 K. The channel shares move accordingly -- at 629 K
R10+R13 goes 8.0 -> 31.1 percent and R18 24.1 -> 0.0 percent -- and the base
x2 root of the seed barely moves, 0.98372 -> 0.98323, so the item's verdict
in section 1 is unchanged by the repair.

### 3.3 R17: the Table-1 Arrhenius is one mechanism of two

Both drift-tube papers were supplied the same day and read
(`references/Johnsen_1980JCP_72_3085.pdf`,
`references/Bohringer_1986JCP_84_1459.pdf`).

| source | range | two-body | three-body |
|---|---|---|---|
| Boehringer & Arnold (1986), abstract and Results | 18-408 K | k2 = 1.1e-13 (300/T)^0.24(+/-0.04) | k3 = 1.6e-30 (100/T)^1.27(+/-0.4) |
| Johnsen, Chen & Biondi (1980), Table I | 78-330 K | (1.1+/-0.1)e-13 at 330 K, (1.5+/-0.15)e-13 at 78 K | (4.4+/-2.0)e-31 at 330 K, (1.8+/-0.4)e-30 at 78 K |
| Johnsen et al. Fig. 3, titled with this reaction (READ from the figure; above 330 K the abscissa is an EFFECTIVE temperature from their earlier elevated-ion-energy data) | to T_eff = 700 K | ~1.45e-13 at 100 K, minimum ~1.05e-13 at 300 K, ~1.35e-13 at 400 K, ~2.0e-13 at 500 K, ~2.5e-13 at 600 K, ~3.0e-13 at 700 K | |
| Schauer et al. (1989) | 15-40 K | 3.0e-14 to 4.9e-14 | |

**The coded Arrhenius alone was wrong at every temperature a molecular base
layer reaches**: 1e-9 exp(-5700/T) gives 5.6e-18 at 300 K against a measured
1.05e-13 and 6.5e-16 at 400 K against 1.35e-13, i.e. 2e4 and 200 low. The
measurements do not switch off because the reaction has TWO mechanisms and
the Arrhenius is only the second: the reactants and products do not
correlate adiabatically (Mahan 1971), so at thermal energy the reaction goes
by tunnelling out of a long-lived He+ -H2 complex (Preston et al. 1978;
Boehringer & Arnold 1986, Discussion), which is weakly and NEGATIVELY
temperature dependent, and only above ~400 K does an over-barrier branch take
over.

**The form now carried is the sum they are:**

    k(R17) = [ 1.1e-13 (300/T)^0.24 - k(R23) ] + 1e-9 exp(-5700/T)

The bracket is the measured two-body TOTAL of He+ + H2 less the radiative
branch this network carries separately as R23, so the coded channels sum to
the measured total rather than each being quoted on its own, and the
branching is stated once and cannot drift from R23. (Johnsen et al.
attribute at least 80 percent of the total to the dissociative channel and
Boehringer & Arnold decline to fix the branching more tightly, so
"total less the radiative branch" is the assignment the measurements
support.) Verification against the only measurement of the rising branch,
Johnsen et al.'s Fig. 3: the sum gives 1.03e-13 at 300 K (read 1.05e-13),
9.6e-14 at 400 K (1.35e-13), 1.01e-13 at 500 K (2.0e-13), 1.60e-13 at 600 K
(2.5e-13) and 3.7e-13 at 700 K (3.0e-13) -- **inside a factor 2 everywhere**,
against the 200 and 2e4 the Arrhenius alone was out by.

Validity, written at the code site: measured 18-408 K; above 408 K the
bracket is carried up with its own weak exponent, which is mild (0.82 from
300 K to 1800 K) and, above ~625 K where the Arrhenius branch passes it, no
longer controls the rate.

**What the Arrhenius branch is, now that its source has been read.** Moses &
Bass (2000) was supplied on 2026-09-15
(`references/Moses_2000JGR_105_7013.pdf`). Their Table A2 row R596 is

    (R596)  He+ + H2 -> H+ + H + He   1.0e-9 e^(-5700/T)
                                      Reference: "Estimate, see text"

so it is NOT a fit to any measurement of He+ + H2. The 5700 K is the
exponent their Sec. 3.2 discusses for a DIFFERENT reaction, their (7)
H+ + H2 -> H2+ + H, which is endothermic and goes only with vibrationally
excited H2: "Although reaction (7) is endothermic, McElroy (1973) first
pointed out that the reaction will be exothermic for vibrationally excited
H2 (for vibrational levels v = 4 or greater). The rate constant for the
reaction of H+ with H2 (v >= 4) has not been measured but is estimated to be
of the order of 1-2 x 10^-9 cm^-3 s^-1, near its maximum possible kinetic
rate", and their Figs. 7 and 8 compare 2e-9 exp(-4900/T), exp(-5700/T) and
exp(-6500/T) as "simplified expressions". So reading the source did NOT
upgrade the >700 K basis from a model uncertainty: it confirmed it, and the
code now says so in the source's own words. A rising term is still REQUIRED
by the data -- Johnsen et al. measure 1.05e-13 at 300 K rising to ~3.0e-13 at
700 K, and the two-body plateau alone gives 8.5e-14 there, a factor 3.5 low
-- so the estimate is carried because the measurements demand a rising term
and nothing measured exists to replace it.

**One thing the same table does settle**: its R597, He+ + H2 -> H2+ + He at
9.35e-15, is the radiative channel, and it is a third independent value for
R23 beside Barlow's 7.2e-15 and Schauer et al.'s measured 7.5e-15 to
1.1e-14. Three sources inside 30 percent.
The three-body channel both papers measure is left out, and since
2026-09-15 the omission is GUARDED CELL BY CELL rather than argued from a
representative density: `rk_R17_Hep_H2_diss` forms k3(T) n / k2 in the cell
it is called for whenever the caller supplies the density, counts the cells
above 0.1 -- ten percent, the size of the two-body measurement's own scatter
-- and warns once with the cell's own n and T if any is. MEASURED on the
certified LHS 1140 b state: no cell trips it, so the omission is a
measurement here and not an assertion.

**And it is the third confirmation that R20 had to go.** Boehringer &
Arnold, Discussion: "Inspection of the potential energy surface shows that it
is not likely that HeH+ is formed from reaction (1) at thermal energy. In ion
beam measurements by Jones et al. no indication of a HeH+ product was found
at a collision energy of 0.15 eV, but H2+ was detected. A very careful search
for the HeH+ product ion was performed by Schindler and a very low upper
limit could be derived for the cross section of reaction (1c),
sigma < 6e-21 cm^2, at energies between 0.13 and 7.5 eV."

### 3.4 Black 1978, checked against the published paper

The published paper was supplied (`references/Black_1978ApJ_222_125.pdf`) and
the coefficient, its derivation and its caveat are quoted verbatim at the
code site. His p. 126, on reaction (12):

> When these cross section data are integrated over a Maxwellian velocity
> distribution and averaged over a thermal distribution of vibrational
> populations of H2+, we find an approximate representation of the rate
> coefficient
>     k12 ~ 3 x 10^-10 exp(-6717/T) cm^3 s^-1.
> This result is in harmony with the rate measured by Neynaber and Magnuson
> (1973).

so the value read from the ADS text is the published one, and it carries an
explicitly THERMAL H2+ vibrational population ("In what follows, the
populations of vibrational states of H2+ are assumed to be thermalized at the
kinetic temperature"). **The direction of that assumption's error is Black's
own, and it is the opposite of what was guessed before the paper was read**:

> If, however, the rate of reaction (10) is as large as suggested recently by
> Bottcher (1976), then some excited vibrational states may have nonthermal
> populations, and the predicted abundances of H2+ and HeH+ may be
> overestimates.

i.e. where the dissociative recombination of H2+ depletes its excited
vibrational levels faster than collisions refill them -- the regime of this
layer -- the coefficient is an UPPER bound. The code comment now says so and
no longer reasons in both directions.

**The other two HeH+ sources of that paper.** Black's rate table lists three,
and the code carries two:

| Black | in the network? |
|---|---|
| (12) H2+ + He -> HeH+ + H, 3.0e-10 exp(-6717/T) | YES, adopted here |
| (13) He* + H2 -> HeH+ + H + e, 1e-9 | YES, as the associative branch (1 - f_penning) of the He(2^3S) + H2 ionization in row 7, with the Garcia Munoz (2025) coefficient and branching rather than Black's single 1e-9 |
| (14) H+ + He -> HeH+* -> HeH+ + hv, 1e-18 | NO |

**(14) is not adopted, and the reason is what it is rather than how small.**
Black introduces it as a proposal -- "Dabrowski and Herzberg (1977) have
proposed that there is a substantial probability of forming HeH+ by
vibrational inverse predissociation" -- argued from a level count, "The
number of quasi-bound levels of HeH+ and its large dipole moment suggest a
large rate coefficient". Adopting an unmeasured 1978 estimate into a channel
is exactly how R20 entered this network. MEASURED so that the omission is
bounded and not merely noted (first cell, 1023 K): 1e-18 n(H+) n(He) = 55
cm^-3 s^-1 against 0.15 for reaction (12), so if adopted it would be the
DOMINANT HeH+ source and would raise n(HeH+) from 7.7e-05 to about
3e-02 cm^-3. It would still change nothing observable -- HeH+ would remain
1e-13 of the helium, and the H2 loss it drives through R18 would be 0.8
percent of the total instead of 0.002 percent. Resolving it needs a computed
radiative-association rate for H+ + He, which exists in the modern literature
(the He+ + H channel of Courtney et al. 2021, ApJ 919, 70, is the nearest
thing in `references/`) and was not obtained here.

## 4. The changes

### 4.1 The seed's root is solved, not divided

**`carrier_h2_chemical_root` now solves the row instead of dividing it.**
g(n2) = production(n2) - loss(n2) is evaluated through `carrier_source`
itself, so the chemistry is still stated once. It is strictly decreasing --
every production term is non-increasing in n(H2) and every loss term is
proportional to it -- and it changes sign between n2 = 0 and the cell's
element ceiling (every hydrogen nucleus the carriers may hold, less the
nuclei the other carriers hold, halved), so a bisection on that bracket
converges to the one root. The bracket closes on a RELATIVE width, because a
root 1e-15 of the ceiling -- the far wind -- is not resolved at all by an
absolute fraction of it. Two cells are handled without a bisection: one
whose production is zero with every nucleus free (no root above zero) and one
whose ceiling still produces (the row wants every nucleus it can have).

The routine now also reports |P - L|/(P + L) at the density it returns, and
the seed's `EXHALE_CARRIER_DEBUG=1` table carries it as a fourth column. It
is the one number that tells a solved root from a divided one. MEASURED on
the `src/tests/molecular_seed` fixture: 6.8e-12 at worst over the twenty-one
radii the report prints, where the old quotient stood at 1.000 by
construction at every saturated cell. A new suite row holds it below 1e-8.

### 4.2 The He+ + H2 channels

Retired: Koskinen R20, He+ + H2 -> HeH+ + H at 4.2e-13 (section 3.2).
Added: H2+ + He -> HeH+ + H at 3.0e-10 exp(-6717/T), Black (1978), default on
because it is then the network's only HeH+ source (3.4); its reaction energy
enters the one reaction-heat ledger by the same formation-table difference as
every other channel, -0.805 eV, and `physics_probe` holds that identity.
Rewritten: R17, from the Arrhenius alone to the sum of the measured thermal
branch and the Arrhenius (3.3). `mk20` is gone as a name, so the compiler
found every site that used it; the coefficient's slot in the packed transfer
vector is reused by the new channel, which is why no index had to be
renumbered.

### 4.3 Files

Files: `src/modules/lower_atmosphere/diffusive_photochemistry.f90` (the root,
the channel record, the writer), `src/modules/nonlinear_system_solver/
System_HeH_mol.f90` (the channels handed out beside the row),
`src/modules/init/molecular_seed_from_atomic_state.f90` (the fourth column),
`src/modules/lower_atmosphere/mol_rates.f90` (the R20 note),
`src/tests/molecular_seed/run.sh` (the new row). Documentation at the
corrected definition: `docs/input_schema.md` appendix D,
`LHS1140b/models/run_case.sh`, `LHS1140b/models/README.md`,
`LHS1140b/MODELS.md` -- all four still said the base layer carried the
handoff, which the code has not done since item L7e withdrew that variant, so
that statement is corrected in the same edit.
`docs/EXHALE_physics_and_algorithms` does not describe the molecular seed, so
that part is unchanged; its molecular-chemistry section carries the two rate
repairs and two new bibliography entries.

Also: `src/modules/lower_atmosphere/mol_rates.f90` (R20 retired, the Black
channel, R17 rewritten, the R23 and Moses & Bass provenance),
`src/modules/nonlinear_system_solver/System_HeH_mol.f90` (rows 2, 4, 5, 7 and
the row scales), `src/modules/lower_atmosphere/molecular_reaction_heat.f90`
(the reaction table row and its rate), `src/modules/radiation/util_ion_eq.f90`
(the ground-singlet helium the new channel's heat needs),
`src/modules/radiation/ionization_equilibrium.f90` and
`src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90`
(the renamed coefficient), `src/tests/physics_probe/
species_formation_energy_table.f90` (the reaction table's test copy), and
`LHS1140b/models/run_case.sh` (a molecular case takes `SEED=`).

**One edit outside this item's file set**, reported rather than hidden:
`src/modules/files_IO/write_output.f90` used `rk_R20_Hep_H2_HeHp` in its
base-cell H2-loss budget and the build cannot link without it, so two lines
there were changed -- the `use` list and the `h2loss(7)` expression, which is
now R17 + R23. That file belongs to another item and the change is the
minimum the link needs.

The H2 row itself is untouched: the channels are assembled from the same
named factors in the same order, so `fvec(4)` is bitwise what it was.

## 5. Verification

### 5.1 Suites

On `EXHALE_L7f.x` / `build_L7f`, each in its own object or output directory
so the tree's `build/` is untouched:

| suite | rows | result |
|---|---|---|
| `carrier_retry` | 143 | pass |
| `carrier_reference_scales` | 14 | pass |
| `carrier_constraint_attribution` | 10 | pass |
| `carrier_returned_state_acceptance` | 36 | pass |
| `certification` | 84 | pass |
| `steady_species_rows` | 195 | pass |
| `molecular_seed` | **25** | pass (24 before, plus the new row) |
| `physics_probe` | 1491 | pass |

The new `molecular_seed` row is
`molecular_seed_local_root_balances_the_row`: the worst |P - L|/(P + L) over
the twenty-one radii the seed report prints, held below 1e-8. MEASURED on
the fixture: **6.811e-12**. It is the row that would have caught the
quotient, which reads a row imbalance of exactly 1 wherever the fit
saturates.

`physics_probe` was run although no rate coefficient changed, because
`mol_heh_rows` was touched; its 1491 rows include the H2 thermochemistry and
the heating-channel closures.

### 5.2 The case: the continuation recipe, and what it leaves

**Adopted 2026-09-15 by the user, and it is the campaign's own recipe:** where
a base layer's composition is set by eddy transport rather than by local
chemistry, a local-chemistry seed cannot start near the fixed point, and the
seed is then a CONTINUATION from the nearest solved state. `run_case.sh`'s
molecular branch now takes `SEED=` like any other case -- a molecular case
with a stated seed maps that state onto the grid instead of converting an
atomic one -- and the catalogue case
`molecular_scalar_gj1132_kzz1e9/HeH2.13` was run that way from
`.L7e/fixed5/HeH2.13/output` on the delivered build.

MEASURED, 40 outer passes:

| | pass 1 | pass 2 | pass 3 | pass 40 |
|---|---|---|---|---|
| carrier row | 1.43e-01 | 6.34e-03 | 1.74e-03 | **2.11e-08** |
| cell-1 mass row | 1.70e-07 | 1.68e-07 | 1.68e-07 | **1.642e-07** |
| momentum | | | | 3.03e-11 |
| energy | 7.46e-07 | | | 7.37e-07 |

**Every equation certifies except one.** The carrier balance, the elemental
transport, the momentum and the energy rows are all inside their tolerances;
the state is refused on the hydrodynamic mass row of cell 1 at 1.642e-07
against 3.7e-09. That is the base-face composition step of section 6, and
nothing in the chemistry repair or the seed touches it -- which is the
item's own conclusion arriving from the other side.

For the record, against the local-chemistry seed on the same recipe: that
one reached 1.77e-02 on the carrier row after forty passes and 1.01e-06 on
the mass row (section 5.2.1).

### 5.2.1 What the corrected seed costs, and why it is not an argument

MEASURED, `.L7f/HeH2.13`: the same recipe and the same 40 outer passes as
`.L7e/fixed5`, seeded with the corrected root (base x2 = 0.984 instead of the
quotient's 0.089). Against the like-for-like control, which is fixed5's own
FIRST solve (`run_dtau0_first.log`, same recipe, same binary but for the root
and the diagnostics, old seed):

| pass | 1 | 4 | 8 | 12 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|
| carrier row, control | 1.00 | 4.80e-01 | 1.17e-02 | 9.98e-04 (p10) | - | - | 5.21e-08 |
| carrier row, corrected seed | 1.00 | 4.09e-01 | 2.54e-01 | 2.26e-01 | 9.51e-02 | 1.89e-02 | **1.77e-02** |
| cell-1 mass row, control | 3.35e-07 | 4.77e-08 | 1.56e-07 | - | - | - | 1.55e-07 |
| cell-1 mass row, corrected seed | **5.20e-08** | **2.15e-08** | 1.01e-06 | 2.57e-06 | 2.54e-06 | 1.33e-06 | 1.01e-06 |

**The corrected seed starts better and ends worse, and the reason is the one
this memo measured.** Its first four passes hold the cell-1 mass row at 2 to
5e-08, six to fifteen times below anything the old seed reached -- because
the base layer is at the handoff's composition and the two sides of the base
face are then the same gas. But that is not where the solve is going: the
fixed point of the full system puts the layer at x2 = 0.33, held there by
eddy mixing with the dissociation region above (section 6.2), so the whole
layer has to travel from 0.98 to 0.33, the movement bound takes that at one
percent a pass, and the column is held while it does. At pass 40 the carrier
row stands at 1.77e-02 against the control's 5.21e-08 and the state is
refused on four entries (cell-1 mass 1.01e-06 of 4.8e-09, cell-2 energy
5.04e-06 of 1.0e-06, the carrier row, and the He/H elemental transport at
3.25e-04).

So: **the root correction is right and the seed built on it is further from
the answer than the wrong one was.** That is not an argument against the
correction -- P/(L/n) is not the row's root, and a seed that happened to sit
near the fixed point for the wrong reason is not a reason to keep an
arithmetic error. It is a statement about what `SEED_X2 = local` can be: a
local CHEMICAL statement is the right seed only where chemistry sets the
composition, and in an eddy-mixed base layer nothing local does.

### 5.3 Regression, and the golden movement

Run on a copy of `backup/regression/` in the scratch directory (the tree's own
matrix was not touched), `REGRESSION_EXE` naming the binary, five molecular
cases. THREE runs, so that what moved can be attributed:

| run | binary | what it carries |
|---|---|---|
| control | `EXHALE_L7f.x` `d153215676a0` | the seed root and the diagnostics only, no rate change |
| measured | `EXHALE_L7f.x` `48db9e876031` | the above plus the three rate repairs |
| attribution | `EXHALE_L18_ctl.x` `922dec0fa183` | the tree as delivered to this item, nothing of it |

**The control says the seed and the record move nothing.**
`mol_carrier`, `mol_base_handoff`, `mol_diffusion` and `mol_lyman_werner`
PASS within 1e-3, worst relative movement 2.5e-05 -- the tree's own distance
from its goldens, the same numbers item L7e reported. `mol_metals` FAILS, at
7.315e-01 on the `cool` column of `Hydro_ioniz.txt` row 32.

**And that failure is not this item's.** The attribution run, which carries
NONE of this item's changes, fails the same case in the same cell against the
same golden: 4.7520e-09 against 1.7700e-08, **7.315e-01**, where the measured
build gives 4.7239e-09, 7.331e-01. Same column, same row, same ratio to three
digits. `mol_metals` was already out of tolerance before this item began.

**The rate repairs move the goldens, and they were meant to.** Attributed by
differencing the measured run against the control run on the same case
(`mol_metals`, worst relative movement over all cells):

| column | control - golden (the tree's own) | measured - control (**this item**) |
|---|---|---|
| HeH+ | 1.27e-02 | **9.998e-01** |
| H3+ | 1.24e-02 | **9.871e-01** |
| H2+ | 4.11e-03 | **9.519e-01** |
| H2 | 1.53e-02 | **6.725e-01** |
| He+ | 1.36e-02 | **2.166e+00** |
| He++ | 7.58e-03 | **2.203e+00** |
| He 2^3S | 2.61e-02 | **2.078e+00** |
| cool | **7.315e-01** | 7.24e-02 |
| heat | 1.93e-02 | 3.00e-02 |
| rho | 1.60e-02 | 2.40e-02 |
| T | 1.19e-02 | 7.16e-03 |
| v | 9.56e-03 | 1.34e-02 |

The signature is exactly the physics: **HeH+ falls to nothing** (8.6e-06
against a golden 7.8e-02 in `mol_carrier`, a factor 9000) because the channel
that made it is gone; **He+ and He++ rise by a factor 2.2** because a helium
sink was removed; H2+ and H3+ follow through the chain; and the hydrodynamic
state moves 1 to 3 percent. The other four cases show the same pattern
(`mol_carrier` and `mol_lyman_werner` lead on HeH+ and H2+; the two handoff
cases lead on the helium ions).

**No golden is refreshed by this item**, and not because the movement is
small: `backup/regression/` is read-only to this item by its own brief, the
tree already fails `mol_metals` for another item's reasons, and item L21's
base-boundary correction landed in the same tree while this ran -- a refresh
now would bake three items into one snapshot. The movement is reported here
so that whoever refreshes knows what is in it.

## 5.4 The nine molecular cases: what certified and why the rest did not

The campaign ran on two machines, this host and `lart3`, over the same NFS
tree. Every case names its binary and its md5 in its own `REPRODUCE.md`.
**The final binary is `EXHALE.x` md5 `c2e9c9990b9f14f1be8cd77abca68945`**,
and every case that certified was re-solved on it from its own state so that
the record carries that md5 and not the one it first certified on.

| # | case | verdict | machine | passes | Mdot log10 [g/s] | what it stands on |
|---|---|---|---|---|---|---|
| 1 | `kzz1e9/HeH0.55` | **certified** | this host, final pass on lart3 | 4, then 3 | 7.64 | its own state |
| 2 | `kzz1e9/HeH2.13` | **certified** | this host, final pass on lart3 | 5, then 2 | 7.91 | its own state |
| 3 | `kzz1e9/HeH9.7` | **certified** | this host, final pass on lart3 | 13, then 2 | 7.95 | its own state |
| 4 | `kzz1e9/HeH0.083` | not solved (2.972e-03 at cell 500, tolerance 1.0e-05) | lart3 | 40, the one continuation allowed | -- | continuation from `wellmixed/HeH0.083` |
| 5 | `wellmixed/HeH0.083` | not solved | this host | 24, refused on the stagnation rule | -- | cold start |
| 6 | `wellmixed/HeH0.55` | not solved | this host | 40 (cap) | -- | cold start |
| 7 | `wellmixed/HeH2.13` | not started | -- | -- | -- | no certified molecular wellmixed state to seed it from |
| 8 | `photochem/HeH2.09` | waiting on its atomic rung | lart3 | -- | -- | the converged column and certified wind of `atomic_photochem_gj1132_kzzprofile/HeH2.09` |
| 9 | `photochem/HeH9` | waiting on its atomic rung | lart3 | -- | -- | the same, at `HeH9` |

**The four not solved, with the reason and nothing else.**

Cases 4, 5 and 6 are one finding, measured in item L7e section 23 and not
argued from: **one scalar movement bound over a column holding a slow H2
front and a far wind.** Every carrier relaxation of every pass of all three
ended on the composition movement bound and none on its own residual; the
bound was cut 5.0e-03 -> 2.5e-03 -> 1.25e-03 -> floor 1.0e-03 with the
movement of a pass falling in proportion; and the cell that ATTAINS the
bound is the H2 front (cell 231 at 1.2565 R_p, its own row measure 6.6e-03;
cell 225 at 1.2309 R_p, 5.1e-03) while the cell that REFUSES is four to five
times worse and decades away in abundance (cell 500 at 29.0 R_p, 2.46e-02;
cell 306 at 1.95 R_p with x2 = 2.7e-06, 2.48e-02). The trajectories say the
same thing from two sides: `wellmixed/HeH0.083` is STALLED (the refusing
cell reached the domain edge at pass 14 and the measure has been flat
between 2.41e-02 and 2.46e-02 since, the front static within one cell), and
`wellmixed/HeH0.55` has NOT ARRIVED (the refusing cell migrates 233 -> 306
at the pace of the x2 = 1e-02 radius, 1.4816 -> 1.5075 R_p over five
passes). `kzz1e9/HeH0.083`, continued from the first of them, is the same
shape one decade lower and it ran its forty passes out: the row at cell 500
fell to 2.07e-03 at pass 11, rose to 3.29e-03 by pass 24 and came back to
2.972e-03 at pass 40, wandering in a band and not descending, against a
tolerance of 1.0e-05. Every relaxation of it also ended on the movement
bound. What differs from the wellmixed pair is WHERE the bound is attained:
here it is the outer wind itself (cell 466 at 16.5 R_p), not the front, so
the same one bound is set by the outer wind in one case and by the front in
the other -- which is the statement that it is one number over a column that
holds two different things.

Case 7 has no admissible seed: the two molecular wellmixed states that would
have provided one are cases 5 and 6.

A regional or coupled relaxation is the obvious next thing and is NOT tried
here. It is opened as plan item L22, a proposal only.

**The two profile cases, and why their definition was changed** (2026-09-16).
They were defined on the 2026-08-30 stored Photochem columns. A molecular
case of that group is solved from the certified wind of its atomic pair, and
`load_IC` admits such a seed only when the two states agree on the reservoir
to a part in 1e6. Measured: the stored 2.09 column and the re-run 2.09
closure rung differ by 1.2e-04 in He/H, twenty times the threshold, and the
seed was refused outright. The group now stands on the converged column of
its OWN closure rung and carries the He/H that column holds at the matching
level (`make_models.py`, `profile_source`; `MODELS.md` section 3), so the
pair cannot drift apart again; the case names are the rungs' names,
`HeH2.09` and `HeH9`, and the former `HeH9.05` directory, which named the
stored column's value and held no result, was removed. The two cases are
armed behind their rungs, which are being re-run on the final binary.

Each of the four carries a `not_solved.md` in its own case directory with
the refusing cell, what its row is made of and what holds it. `MODELS.md`
sections 7 and 8 are written by `models/status.py` from the case directories
and carry the verdict alone, so the reason lives beside the case and in the
hand-written group table of `MODELS.md` section 3.

No case definition was loosened to make anything certify, no cell was
excluded from certification, and no local equilibrium was imposed anywhere.

## 6. What this leaves open, and the proposal

The item's own question is answered and the answer removes the reason the
base BC's range was to be reopened *for a chemistry disagreement*: there is
no such disagreement at the base. What remains at the base face is a
DISCRETIZATION question, and it is a real one. These are the measured facts
it has to be decided on; nothing here is implemented.

**1. The two sides of the base face are treated by two different rules, and
one of them has no mixing in it.** The ionization sweep imposes the handoff
partition on the GHOST rows alone (`ionization_equilibrium`, the
`base_h2_composition_imposed` block, `j <= 0`). The carrier operator then
takes the base face's ADVECTIVE flux with that ghost composition (the face
composition is reconstructed from the ghosts wherever the face mass flux
flows inward) and its DIFFUSIVE flux as exactly zero -- faces 0 and N carry
no diffusive flux, a documented decision whose stated reason is that "the
handoff supplies the composition of the gas that flows in, not a mixing rate
across a boundary whose gradient is set by the ghost spacing and by a K_zz
that describes an unresolved region".

**2. One cell above that face the same gas IS eddy-mixed, hard.** MEASURED at
cell 1 of `.L7e/fixed5/HeH2.13`: the diffusive divergence of the H2 row is
5.70e+05 cm^-3 s^-1 against a net chemical source of 6.77e+04, and with the
base face carrying none of it, all of that is the 1-2 face. The eddy time
across one cell, (2.18e5 cm)^2/1e9 = 48 s, is five decades below the H2
chemical time. So the model says "no mixing at 1 microbar" and "mixing 200 km
higher, on a 48 s time scale" at the same time, and a composition step across
the face is the only way both can hold.

The cell's whole H2 budget closes on those three terms and says where the
molecule comes from and where it goes (cm^-3 s^-1):

    advective supply from the handoff  5.020e+05
    net chemical production            0.677e+05
    ------------------------------------------- 
    diffusive export upward            5.696e+05

-- the boundary's advective inflow, carrying the handoff's x2 = 0.99998, is
already the dominant SOURCE of H2 in the first cell, and the eddy diffusion
into the dissociation region above is what removes it. The chemistry is a
tenth of the budget.

The layer is flat in composition, which is what a well-mixed layer looks
like: x2 = 0.3315, 0.3308, 0.3208, 0.2860 at cells 1, 25, 100 and 142. And
it is NOT where its own chemistry would put it: solving the quadratic
k15 (A - 2 n2)^2 + P_other = L_rate n2 on the channels the dump gives at
cell 1, with A = 2.2848e12 the hydrogen available and the ion background
frozen, gives x2 = 0.82 (ARITHMETIC on the dumped terms, not a run of
`carrier_h2_chemical_root`, which is called only by the seed). Eddy mixing
with the dissociation region above holds the layer a factor 2.5 below its own
chemical root, and it is the mixing and not the chemistry that decides how
far the base cell can sit from the boundary's composition.

**3. What the step costs.** MEASURED (item L7e section 21.2, re-read here):
ghost x2 = 0.99998 against cell 1 at 0.3315, so 2.63 against 2.964 particles
per hydrogen nucleus -- 12.8 percent in the mean molecular mass across one
face -- with rho stepping 26 percent and v from -0.827 to -0.137 cm/s. The
certified ATOMIC case of the same planet and He/H, which has no composition
step at all, steps 4.2 percent in rho and its worst mass row anywhere is
1.487e-09.

**The proposal, in the order the measurements rank it.**

*(a) Make the base face's two fluxes state the same thing.* Either the
handoff states the composition AT the boundary -- in which case the face
carries the eddy-diffusive flux the same way every other face does, the
composition is continuous across it, and the column relaxes to the handoff
over a resolved diffusive boundary layer -- or it states only the composition
of the gas that flows IN, in which case the ghost partition must not be
imposed either and the ghost takes the first cell's composition
(df/dr = 0), which is exactly what the operator already does when no handoff
states a partition. The present pairing, a Dirichlet ghost with a zero
diffusive flux, is the one combination in which the two sides of the face
are different gas by construction. The second option is the smaller change
and keeps the stated reason for the zero flux; its cost is that the
lower-atmosphere model's H2 at 1 microbar then informs only the element
reservoirs, the base density and the entropy, and not the partition.

The budget above argues for the FIRST. The base face already supplies the
first cell's H2, through the advective term and at the handoff composition,
and what carries it away is a diffusive flux the same face is not allowed to
have. Giving that face its eddy-diffusive flux adds a source on the side the
gas comes from, which is the side the physics puts it on: a K_zz of 1e9 does
not stop at 1 microbar, and `base.inp` states `Kzz_base` for that level
itself. The objection recorded against it -- that the gradient would be set
by the ghost spacing -- is a statement about accuracy, not about existence,
and it is answerable: the ghosts sit at known radii and the handoff states
the composition AT the base level, so the one-sided gradient is as well
defined as any other face's.

*(b) Before either, settle whether the two models describe the same gas at
the matching level at all -- and the answer is that they do not, for a reason
that is now traced.* MEASURED on the written state: the two base ghost rows
carry T = 879 and 873 K, cell 1 carries 1023 K, and `base.inp` states the
lower atmosphere's base at **226 K** and 1 microbar. The certified atomic
case of the same planet carries 419 and 414 K in the same rows.

**Why `T_base` does not set the run's base T** (traced in `base_boundary.f90`,
and the reason is the boundary condition and not a heating): `Base BC:
pressure` states a (p, s) reservoir, but the entropy reaches the gas only
through the ENTERING characteristic. `characteristic_base_face_state` forms
the face density as `rho_b = (1 - w_rev) rho_res + w_rev rho_rev` with
`w_rev = characteristic_branch_weight(M_i/base_face_mach_blend)`,
`base_face_mach_blend` = 1e-6 and `rho_rev` the INTERIOR isentrope continued
to the face; the ghosts are then the volume average of the hydrostatic
isentrope through that face state. The weight is the cubic smoothstep, 1 at
x <= -1 and 0.5 at x = 0, so ANY inflowing base face gives w_rev >= 0.5 and
|M_i| >= 1e-6 gives w_rev = 1 exactly. MEASURED: cell 1 carries
v = -0.137 cm/s at 1023 K (molecular) and -0.113 cm/s at 417.6 K (atomic),
i.e. |M_i| ~ 7 to 8e-7 and **w_rev ~ 0.95 to 0.98** -- the face state is 95
to 98 percent the interior's own isentrope and the 226 K reservoir is
blended almost entirely out. The module's design note says that window "must
be far BELOW the physical operating point" and quotes the hot-Uranus gate at
a face Mach number of 6e-6; LHS 1140 b's base velocity is ten times smaller,
so this planet sits INSIDE the window. Not determined here: whether the
negative base velocity is physical or the odd-even cell artefact recorded for
the hot Uranus, which would make the cell velocity a poor proxy for the face
one.

**And the heating that holds the molecular base at 1023 K is measured**
(plan item L7g, opened for it): at cell 1 the heating is 97.7 percent
`heat_mol_chem`, the molecular reaction heat, 5.205e-07 erg cm^-3 s^-1 of
5.327e-07, against a TOTAL cooling of 1.101e-08 -- a factor 48. Turning on
both infrared coolants (`Base IR field`, `Molecular IR bands`) raises the
cooling to 1.225e-08, +11 percent, and moves the base temperature by 0.2 K
and the base x2 by 0.002: **the molecular model is not incomplete in that
sense, and the IR channels are not what is missing.** The R15 three-body
association alone runs at 6.92e+04 cm^-3 s^-1 there and releases 4.48 eV per
H2 formed, 4.95e-07 erg cm^-3 s^-1, essentially the whole of that heating
term, while the H2 infrared bands radiate 1.17e-09 of it.

*(c) The helium-ion rates (section 3 and section 7) are the entry to resolve
before any molecular result of a helium-rich case is quoted.* They are 56
percent of the base H2 destruction, and R20, the largest of them, is not
verified against its source.

## 7. Publications not obtained

Five papers were asked for during this item and ALL FIVE were supplied the
same day; none of the He+ + H2 entries is blocked any longer. They are now
`references/Schauer_1989JCP_91_4593.pdf` (settled R20, section 3.2),
`references/Johnsen_1980JCP_72_3085.pdf` and
`references/Bohringer_1986JCP_84_1459.pdf` (settled R17, section 3.3) and
`references/Black_1978ApJ_222_125.pdf` (verified the replacement channel and
its caveat, section 3.4).

What remains unobtained is smaller and is named so it is not forgotten.

Moses & Bass (2000) was supplied the same day and is now
`references/Moses_2000JGR_105_7013.pdf`; what it says about R17's Arrhenius
branch is in section 3.3. One item remains.

**A computed radiative-association rate for H+ + He -> HeH+ + hv**, to
settle whether Black's reaction (14) belongs in the network (section 3.4).
Black's 1e-18 is a 1978 estimate from a level count, not a measurement; the
modern computed rates exist (Courtney et al. 2021, ApJ 919, 70, is in
`references/` but is the He+ + H channel, not this one). Nothing in this
item's results depends on it: adopted at 1e-18 it would move n(HeH+) from
7.7e-05 to 3e-02 cm^-3 and the H2 loss by 0.8 percent.

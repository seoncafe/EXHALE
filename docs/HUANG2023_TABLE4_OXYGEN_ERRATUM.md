# Huang et al. (2023) Table 4: the oxygen charge-exchange rows are exchanged

- Date: 2026-07-31
- Paper: Huang, C., et al. 2023, ApJ, 951, 123, "A Hydrodynamic Study of the
  Escape of Metal Species and Excited Hydrogen from the Atmosphere of the Hot
  Jupiter WASP-121b" (`references/Huang_2023_ApJ_951_123.pdf`)
- Table: Table 4, "Charge Exchange Rates", 65 reactions, column headed
  **Reactants**
- Affected rows: `O + H+` and `O+ + H`
- Found while auditing the charge-exchange data of MoCHII, which transcribes
  this table; recorded here because EXHALE draws on the same literature for
  metal charge exchange with hydrogen.

## Summary

The two oxygen rows carry the right rate coefficients with the reactant labels
**exchanged**. The Boltzmann factor `exp(-227/T)` is printed on the `O+ + H`
row; the physics puts it on `O + H+`. Read as printed, under the table's own
column header, the pair violates detailed balance by a factor of 1.47 at
8000 K. Exchanged, it satisfies detailed balance to 9% and matches Cloudy
c25.00's independent fits to 1-3% per direction; MOCASSIN confirms the
assignment (which direction is the barrierless, exothermic one) but not the rate
magnitude, since its oxygen recombination rate is a Kingdon & Ferland (1996)
fit, 0.54x Huang's at 8000 K.

The magnitude, 227 K, is the oxygen-hydrogen ionization-potential difference to
0.3%, so this is a placement error and not a numerical one.

## The rows as printed

```
O  + H+ :  2.08e-9 T4^0.405 + 1.11e-11 T4^-0.458
O+ + H  : (1.26e-9 T4^0.517 + 4.25e-10 T4^6.69e-3) exp(-227/T)
```

## Why this is the wrong way round

`IP(O I) = 13.6181 eV` against `IP(H I) = 13.5984 eV`. Oxygen holds its
electron more tightly, so

- `O0 + H+ -> O+ + H0` is **endothermic** by 0.0196 eV, i.e. `dE/k = 227.7 K`,
  and must carry `exp(-dE/kT)`;
- `O+ + H0 -> O0 + H+` is **exothermic** and must not.

## The paper's own helium row is the control

Three species in Table 4 exchange with hydrogen and hold their electron more
tightly than hydrogen does: helium, oxygen, and nitrogen. Nitrogen's
`N + H+ / N+ + H` pair (Lin et al. 2005) sits on the table's continuation page,
with `IP(N I) = 14.534 eV > IP(H I)`. Of the three, only helium carries an
explicit Boltzmann factor in its printed rate, so helium alone is the control
for the convention. Nitrogen's two directions are numerical `exp[polynomial in
lnT]` fits with no separable exponential factor, so no placement error is
possible there, and its rates are negligible in any case (`3.3e-19` and
`1.1e-14 cm^3 s^-1` at 8000 K). Helium is the same case as oxygen, and the table
gets it right:

```
He + H+ : 1.75e-11 (T/300)^-0.75 exp(-12.75/T4)
He+ + H : 1.25e-15 (T/300)^0...        [no exponential]
```

`IP(He I) - IP(H I) = 24.5874 - 13.5984 = 10.989 eV`, i.e. **12.752** in units
of `1e4 K` — against the printed 12.75. Exact, and on the `He + H+` row, which
is the endothermic reactant pair. Carbon is the opposite case,
`IP(C) = 11.26 < IP(H)`, and its exponential sits on `C+ + H`, again the
endothermic side.

So the table's convention is reactant labeling and it is applied correctly for
helium and for carbon. Oxygen is the one row inconsistent with it.

| row | IP(X) [eV] | endothermic reactant pair | exponential printed on | verdict |
|---|---|---|---|---|
| He | 24.5874 | `He + H+` | `He + H+`, 12.75 vs 12.752 | correct |
| N | 14.534 | `N + H+` | none (polynomial fit) | no explicit factor; negligible |
| **O** | **13.6181** | **`O + H+`** | **`O+ + H`, 227 K** | **exchanged** |
| C | 11.2603 | `C+ + H` | `C+ + H` | correct side |

## Two independent codes agree with the corrected assignment

**Cloudy c25.00**, `source/atmdat_char_tran.cpp:112-176`, carries TableCurve
fits for both directions and assigns them the other way from the printed table:

| T [K] | `k(O0+H+)` Cloudy | `k(O++H0)` Cloudy | ratio | detailed balance requires |
|---|---|---|---|---|
| 4000 | 1.112e-09 | 1.429e-09 | 0.778 | 0.839 |
| 8000 | 1.499e-09 | 1.874e-09 | 0.800 | 0.864 |
| 20000 | 2.230e-09 | 2.709e-09 | 0.823 | 0.879 |

Cloudy's pair satisfies detailed balance to 7% across the range, with the
ionizing direction the **smaller** of the two, as an endothermic channel must be.

**MOCASSIN 3.x**, `source/update_mod.f90:1671`, stores
`chex(8,1) = (1.04, 3.15e-2, -0.61, -9.73)` with **no Boltzmann factor**, and
its `chex` comments label the product ion, so that entry is `O+ + H0 -> O0 + H+`
— the exothermic direction, correctly barrierless. Its rate magnitude is a
Kingdon & Ferland (1996) fit, 0.54x Huang's value at 8000 K, so MOCASSIN
confirms only the assignment — which direction is barrierless — and not the
rate.

## The decisive numerical check

Huang's two coefficients match Cloudy's two coefficients when, and only when,
the labels are exchanged:

| T [K] | Huang `O + H+` | Cloudy `O+ + H0` | ratio | Huang `O+ + H` | Cloudy `O0 + H+` | ratio |
|---|---|---|---|---|---|---|
| 4000 | 1.452e-09 | 1.429e-09 | 1.016 | 1.140e-09 | 1.112e-09 | 1.026 |
| 8000 | 1.913e-09 | 1.874e-09 | 1.021 | 1.504e-09 | 1.499e-09 | 1.003 |
| 20000 | 2.762e-09 | 2.709e-09 | 1.019 | 2.205e-09 | 2.230e-09 | 0.989 |

Agreement is 1-3% under the exchange (the largest deviation is 2.6%, at
4000 K). The values are right; only the row labels are swapped.

The agreement is not two independent determinations converging. Huang's oxygen
rows carry footnote `g`, which is Stancil et al. (1999, A&AS 140, 225), and
Cloudy c25.00's two oxygen fits are TableCurve fits of that same Stancil
calculation (stated in the comments of `source/atmdat_char_tran.cpp`). So
Huang's two expressions and Cloudy's two fits descend from one calculation, and
the ratios above show that the two sources attach the labels the opposite way —
which is an independent confirmation of the exchange, not an independent
measurement of the rate.

## Effect of using the table as printed

Detailed balance requires
`k_i/k_r = [g(O+) g(H0)] / [g(O0) g(H+)] * exp(-dE/kT) = (4*2)/(9*1) * exp(-227.7/T)`,
which is at most `8/9 = 0.889` and cannot exceed it for an endothermic ionizing
direction.

| | measured `k_i/k_r` at 8000 K | required | error |
|---|---|---|---|
| as printed | 1.272 | 0.864 | **1.47x** |
| exchanged | 0.786 | 0.864 | 0.91 |

Using the rows as printed makes `O0 + H+` too fast by about 27% and
`O+ + H0` too slow by about the same, so it **overestimates the ionized fraction
of oxygen** wherever charge exchange controls it — which, because the reaction
is near-resonant, is wherever hydrogen is partially neutral. Oxygen's ionization
is locked to hydrogen's through this pair, so the error propagates into any
quantity that depends on `x(O+)/x(O0)`.

## A related caution: deriving one direction from the other

A code can avoid inheriting a table error like this by tabulating one direction
and **deriving** the other from detailed balance. MOCASSIN does exactly that
(`source/update_mod.f90:1741-1744`). This is worth stating precisely, because it
is easy to overstate: derivation does not make a code more or less error-prone.
It changes **which** errors are possible, and **whether** they can be found.

Derivation needs an energy defect, so a wrong one is possible — and since the
pair is then self-consistent by construction, no detailed-balance test can
reveal it. Independent fitting uses no energy defect, so that error cannot
occur; a wrong pair of rates can, and that is precisely what such a test finds.

Both arrangements failed on oxygen in practice. MOCASSIN tabulates two defects:
nitrogen is right (10863 K against `[IP(N)-IP(H)]/k = 10858.4`, 0.04%, which
also confirms the units) while **oxygen is 2205 K against a true 227.7, wrong by
a factor 9.7**, making its derived reverse rate 22% low at 8000 K and
undetectable by construction. The table documented above is the other failure,
1.47x, and it was found by the test that derivation would have disabled.
Neither method prevented an error; one made it findable.

The 2205 has no obvious provenance. The true 227.7 K coincides with the O I
`3P2-3P1` fine-structure interval at 227.4 K, but none of the natural oxygen
energies — `1D-3P` at 22826 K, O II `2D-4S` at 38573 K, the O/He defect at
127294 K — is 2205.

Two further limits on derivation. The Boltzmann factor carries only the
asymptotic energy difference, while real charge transfer can have a
curve-crossing barrier — Huang's own carbon row is 170000 K against an
ionization-potential difference of 27132 K, and deriving carbon's reverse from
its forward rate would have erased that. And ground-term statistical weights are
an approximation, since charge transfer often proceeds through particular
fine-structure or excited channels: adequate for checking at a tolerance,
weaker as a generator of rates.

The safer arrangement is the one that caught this: keep both directions from the
literature, and test them against each other.

## What to check in a code that uses this table

1. Whether the two oxygen rates are assigned as printed. If the code puts
   `exp(-227/T)` on the `O+ + H0` channel, it has inherited the error.
2. More generally, whether any pair of directions taken from a table satisfies
   `k_i/k_r = [g(X+) g(H0)]/[g(X0) g(H+)] exp(-dE/kT)`. This audit was done by
   hand and the error had stood undetected; it is cheap to make it automatic.
   MoCHII now runs it as a gate over all its pairs
   (`tests/charge_exchange/check_detailed_balance.f90`), which fails on the
   printed oxygen assignment at every temperature from 5000 to 20000 K.
3. Note that below about 5000 K the Kingdon & Ferland forms are usually
   evaluated with the temperature clamped into their validity range, which
   freezes one direction's Boltzmann factor while the other keeps running.
   Detailed balance is then broken by the clamp rather than by the data, so a
   check of this kind should not be applied below the clamp.

## Status in EXHALE

EXHALE had inherited the printed assignment verbatim: in
`src/modules/radiation/charge_exchange.f90`, case A13 (`O + H+`, the endothermic
ionizing channel) carried no barrier and case A14 (`O+ + H`) carried
`exp(-227/T)`. The two rate expressions have been swapped so the barrier sits on
the endothermic `O0 + H+` direction, with a comment at the site and the
correction propagated to `docs/charge_exchange_table4.md`. An audit of every
other `X0 <-> X+ + H` pair (Mg, Fe, Si, C, N, S, Na) and the `X + He` pairs
found the barrier on the correct side with `Ek` matching `|IP(X)-IP(H)|` to
within a percent (Mg/Fe/Si/S satisfy detailed balance to a fit-scatter factor of
~1.3; He/Si/C/O + He are correct). The N/Na/K rows are radiative
charge-transfer fits, not resonant two-body exchange, and are not subject to this
detailed-balance test; their rates are negligible. Oxygen was the only swap.

## Not verified here

Whether the error is in the published table or in its typesetting, and whether
the authors' own code used the correct assignment. Only the printed table was
read. Nothing here should be taken as a statement about the WASP-121b results
of that paper, which depend on how the rates were used rather than on how they
were tabulated.

Carbon deserves a separate look: its exponential is on the correct side, but the
printed argument is 17 in units of `1e4 K`, i.e. 170000 K, against an
ionization-potential difference of 27132 K. That may be a genuine
curve-crossing barrier rather than the asymptotic defect. The rate is
2.3e-23 cm^3 s^-1 at 8000 K either way, so nothing in a nebular or atmospheric
model turns on it.

# Huang et al. (2023) Table 4: Charge Exchange Rates (verified transcription)

This is the authoritative, image-verified transcription of Table 4 of Huang,
Koskinen, Lavvas & Fossati (2023, ApJ 951, 123), "Charge Exchange Rates," used
to implement charge exchange in `EXHALE` (Phase 1d). Every coefficient and
`exp()` argument below was read directly from the rendered PDF (pp. 25-26) at
high magnification, not from a lossy text extraction. Where the layout could be
read two ways, the resolved reading is noted.

Notation (exactly as in the paper):
- `T4 = T / 1e4 K`; `T` is in Kelvin; `lnT = ln(T[K])`.
- Bracketed factors use `exp(-a*T4)` (a *multiplicative* T4 term).
- Trailing Boltzmann factors use `exp(-a/T4)` (a *divided* T4 term, i.e. an
  activation barrier `a*1e4 K`). The two forms look identical in the typeset
  table; they are distinguished here by position and by physical role.
- Rates are in cm^3 s^-1.

A charge-exchange reaction transfers one electron. For a row written
`A + B -> ` the products are implied by single-electron transfer and overall
charge conservation, e.g. `X + H+ -> X+ + H` and `X+ + H -> X + H+`. The two
directions are listed as separate rows with independently fitted rates (the
reverse rates were obtained by Huang via microscopic balance, their Eq. 1-2,
where direct data were unavailable).

## Reference key (Table 4 footnotes)

| Tag | Source |
| --- | --- |
| GJ07 | Glover & Jappsen (2007) |
| KF96 | Kingdon & Ferland (1996) |
| RV72 | Rutherford & Vroom (1972) |
| ND87 | Neufeld & Dalgarno (1987) |
| Satta13 | Satta et al. (2013) |
| S99 | Stancil et al. (1999) |
| S98 | Stancil et al. (1998) |
| PH80 | Prasad & Huntress (1980) |
| Z05 | Zhao et al. (2005) |
| Chenel10 | Chenel et al. (2010) |
| Z04 | Zhao et al. (2004) |
| Kwolek19 | Kwolek et al. (2019) |
| L05 | Lin et al. (2005) |
| D01 | Dutta et al. (2001) |
| W02 | Watanabe et al. (2002) |
| Lavvas14 | Lavvas et al. (2014) |

## Group A: charge exchange with H / H+ (the Phase-1d default active set)

These couple the H ionization fraction to each metal's ionization balance and
are the only reactions enabled by default (`cx_full` off). Note **Ca has no
direct H charge exchange in Table 4**: it appears only in Group D.

> [2026-08-15: "the only reactions enabled by default" is no longer accurate.
> The He<->H pair of Group B (B1, B2) is also on by default, under its own key
> `He_H_charge_exchange` (default `.true.`) and independently of `cx_full`; it
> is applied by dedicated routines in every ionization system that contains He.
> `cx_full` still gates Groups C and D, and the *metal* part of Group A is
> still what it enables by default. See
> `src/modules/radiation/charge_exchange.f90`.]

| # | Reaction | Ref | Rate (cm^3 s^-1) |
| --- | --- | --- | --- |
| A1 | Mg + H+ | KF96 | `9.76e-12 * T4^3.14 * [1 + 55.54*exp(-1.12*T4)]` |
| A2 | Mg+ + H | KF96 | `2.95e-12 * T4^3.28 * [1 + 55.54*exp(-1.12*T4)] * exp(-6.91/T4)` |
| A3 | Mg+ + H+ | KF96 | `7.6e-14 * [1 - 1.97*exp(-4.32*T4)] * exp(-1.67/T4)` |
| A4 | Mg2+ + H | KF96 | `8.58e-14 * T4^(2.49e-3) * [1 + 0.0293*exp(-4.33*T4)]` |
| A5 | Fe + H+ | RV72 | `4.0e-9` |
| A6 | Fe+ + H | RV72 | `1.16e-9 * (T/300)^0.072 * exp(-6.61/T4)` |
| A7 | Fe+ + H+ | ND87 | `2.3e-9 * exp(-3.0/T4)` |
| A8 | Fe2+ + H | ND87 | `1.26e-9 * T4^0.0772 * [1 - 0.41*exp(-7.31*T4)]` |
| A9 | Si + H+ | GJ07 | `7.41e-11 * (T/300)^0.85` |
| A10 | Si+ + H | GJ07 | `4.71e-11 * (T/300)^0.95 * exp(-6.32/T4)` |
| A11 | Si+ + H+ | KF96 | `4.1e-10 * T4^0.24 * [1 + 3.17*exp(-4.18e-3*T4)] * exp(-3.18/T4)` |
| A12 | Si2+ + H | KF96 | `1.26e-9 * T4^0.24 * [1 + 3.17*exp(-4.18e-3*T4)]` |
| A13 | O + H+ | S99 | `(1.26e-9 * T4^0.517 + 4.25e-10 * T4^(6.69e-3)) * exp(-227/T)` |
| A14 | O+ + H | S99 | `2.08e-9 * T4^0.405 + 1.11e-11 * T4^(-0.458)` |
| A15 | C + H+ | S98 | `1.31e-15 * (T/300)^0.213` |
| A16 | C+ + H | S98 | `6.3e-17 * (T/300)^1.96 * exp(-17/T4)` |
| A17 | N + H+ | L05 | `exp[-35.4 + 1.94*lnT - 0.154*lnT^2 - 6.3e-3*lnT^3 - 1.16e-3*lnT^4]` |
| A18 | N+ + H | L05 | `exp[-40.1 + 6.4*lnT - 1.75*lnT^2 + 0.18*lnT^3 - 5.96e-3*lnT^4]` |
| A19 | S + H+ | Z05 | `exp[-50 + 13.3*lnT - 2.77*lnT^2 + 0.243*lnT^3 - 7.24e-3*lnT^4]` |
| A20 | S+ + H | Z05 | `exp[-50.14 + 13.3*lnT - 2.77*lnT^2 + 0.243*lnT^3 - 7.24e-3*lnT^4 - 3.76/T4]` |
| A21 | Na + H+ | D01,W02 | `exp[48.2 - 43.0*lnT + 8.74*lnT^2 - 0.77*lnT^3 - 0.0256*lnT^4]` |
| A22 | Na+ + H | D01 | `exp[46.6 - 41.6*lnT + 8.02*lnT^2 - 0.65*lnT^3 + 0.0197e-3*lnT^4 - 9.817/T4]` |
| A23 | K + H+ | W02 | `exp[-27.8 + 0.125*lnT + 0.0663*lnT^2 - 0.0237*lnT^3 - 1.36e-3*lnT^4]` |

Notes on Group A:
- **A13/A14 (O + H+ / O+ + H)**: the two rate coefficients are shown here with
  the reactant labels corrected relative to the printed Table 4, which exchanges
  them. `IP(O I) = 13.6181 eV > IP(H I) = 13.5984 eV`, so `O0 + H+ -> O+ + H0`
  is endothermic (`dE/k = 227.7 K`) and carries `exp(-227/T)`, while `O+ + H0`
  is exothermic and does not. As printed the pair violates detailed balance by
  1.47x at 8000 K; corrected it satisfies it to ~9% and matches Cloudy c25.00
  to 1-3% per direction. See `HUANG2023_TABLE4_OXYGEN_ERRATUM.md`.
- **A16 (C+ + H)** uses `exp(-17/T4)`, i.e. a 1.7e5 K activation barrier. This is
  the endothermic reverse of the fast `C + H+` charge transfer and is utterly
  negligible (~1e-24 at T4 = 1); ported verbatim per the paper.
- **A17-A23** are `exp[polynomial in lnT]` fits (Lin/Zhao/Dutta/Watanabe). The fits
  blow up only for T -> 1 K (lnT -> 0), which never occurs; the implementation
  caps the evaluated coefficient at a physical maximum to be safe.
  [2026-09-24: "tiny over the model range" was not a property of a rate but a
  sign that the rows are not rates: A17, A21, A22 and A23 fall with T faster
  than any Maxwellian average of a cross section can (the threshold law, second
  block below); A19/A20 are of ordinary size (1.8e-11 at 1e4 K).]
- **K has only the forward `K + H+`** in Table 4; no `K+ + H` reverse is tabulated.

> [2026-09-24: what the code now does with the Group A rows, each checked
> against the paper Huang et al. cite (and, where that paper is a fit, the
> calculation behind it). The Table 4 rows above stay as printed; this block
> records the departures. Full account:
> `md/charge_exchange_detailed_balance_20260924.md`.
> - A1 is Kingdon & Ferland (1996) Table 3 row Mg0, their direct fit to
>   Allan et al. (1988, MNRAS 235, 1245) Table 1 (READ: 1.867e-10 against the
>   tabulated 1.878e-10 cm^3 s^-1 at 1e4 K). Allan et al. form Mg+ in the
>   3p 2Po and 4s 2S states, so A2 is not the reverse of their channel from
>   the populated Mg+ 3s; **A2 is not carried** (Kingdon & Ferland give no
>   reverse, "negligible due to the energetics").
> - A3 is **computed from A4 by detailed balance** with internal partition
>   functions (it was Kingdon & Ferland's Table 3 row Mg+, their refit of the
>   same relation over 1e4 - 3e5 K). A4 is Kingdon & Ferland Table 1 row Mg+2
>   (Butler & Dalgarno 1980, ApJ 241, 838, Table 1A, READ).
> - A5/A6 are carried as printed: Rutherford & Vroom (1972) was not read.
>   [Superseded 2026-09-25: A5 from their Table I, A6 not carried.]
> - A7 is **computed from A8 by detailed balance**; A8 is Kingdon & Ferland
>   Table 1 row Fe+2 (Neufeld & Dalgarno 1987, not read). The printed A7,
>   `2.3e-9 * exp(-3.0/T4)` citing ND87, is not carried.
> - A9 is Glover & Jappsen (2007) R42 with both of its branches,
>   `5.88e-13 * T^0.848` (T <= 1e4 K) and `1.45e-13 * T` (T > 1e4 K); Table 4
>   carries the first alone. R42 cites Kingdon & Ferland (1996), whose tables
>   have no Si0 row; the fit reproduces Kimura et al. (1996, ApJ 473, 1114)
>   Table 1 (READ), whose channel ends in the metastable Si+ 3s3p2 4P term.
>   **A10 is not carried**: it is the reverse of that channel written for
>   ground-term Si+.
> - A12 is corrected to Kingdon & Ferland (1996) Table 1 row Si+2 and
>   Glover & Jappsen (2007) R46: `1.23e-9 * T4^0.24 * [1 + 3.17*exp(+4.18e-3*T4)]`
>   (Table 4 prints 1.26e-9 and the exponent with a minus sign). It is within
>   3% of Gargaud et al. (1982, A&A 106, 197) Table 3 (READ). A11 is
>   **computed from A12 by detailed balance**, as Gargaud et al. do.
> - A15 is the tabulated rate of Stancil et al. (1998, ApJ 502, 1006) Table 1,
>   reactions (2) + (3) (non-radiative plus RADIATIVE), interpolated in
>   log-log. Table 4's `1.31e-15 * (T/300)^0.213` is the first term of their
>   fit only; with the printed exp(-c2/T) the second term does not reproduce
>   their table above 1e4 K (3.2e-15 against 1.16e-14 at 2e4 K). A16 is their
>   reaction (1) fit, `6.08e-14 * T4^1.96 * exp(-1.7e5/T)` (Table 4's 6.3e-17
>   (T/300)^1.96 is its rounding), calculated for 3e4 - 1e7 K and
>   extrapolated below; the 17/T4 barrier is their fit parameter, not a
>   transcription error.
> - A19/A20 are carried as printed: Zhao et al. (2005, Phys. Rev. A 71,
>   062713) could not be obtained. [Superseded 2026-09-25: A19 from their
>   Table III, A20 by detailed balance.]]

> [2026-09-24, second round (`md/charge_exchange_detailed_balance_20260924.md`
> section 10):
> - **A18 is Kingdon & Ferland (1996) Table 1 row N+1** (READ from the ADS scan:
>   a = 1.01(-3), b = -0.29, c = -0.92, d = -8.38, 1e2 - 5e4 K), their fit to
>   Butler & Dalgarno (1979, ApJ 234, 765; READ: 1.0e-12 at 1e4 K, 1.2e-12 at
>   1e3 K; product N(4S), the only term the 0.94 eV reaches). **A17 is computed
>   from A18 by detailed balance**; Kingdon & Ferland's own Table 3 row N0 is
>   the same relation with ground-term weights (4.55(-3) = 4.5 x 1.01(-3)). The
>   printed rows (Lin et al. 2005, Phys. Rev. A 71, 062708, not obtained) are
>   not carried: the endothermic N + H+ fell from 4.1e-18 at 5e3 K to 8.1e-22
>   at 2e4 K, and at 1e4 K the two rows were 2e7 and 70 times below these.
> - [Superseded 2026-09-25, next block: A21 and A23 are carried from the
>   papers, A22 is not.] **A21, A22 and A23 are not carried.** As printed (ln T, T in K) A21 falls
>   from 1.3e-128 to 3.2e-224 over 5e3 - 2e4 K, A22 from 1.4e-64 to 4.8e-94,
>   A23 from 8.7e-16 at 1e3 K to 3.8e-22 at 1e4 K: faster than T^(-3/2), which
>   no Maxwellian average of a cross section can do. No other reading of the
>   variable (log10 T, ln(T/300), T in eV, T/1e4 K) gives a rate either. Dutta
>   et al. (2001, Phys. Rev. A 63, 022709) and Watanabe et al. (2002, Phys.
>   Rev. A 66, 044701) could not be obtained (APS paywall); until they are,
>   the code has no Na-H or K-H charge transfer (it had none in effect
>   before: A21 was 4.4e-171 at 1e4 K).]

> [2026-09-25, third round (`md/charge_exchange_detailed_balance_20260924.md`
> section 11), with the originals supplied:
> - **A5 is Rutherford & Vroom (1972) Table I** (READ: 73.5, 53.4, 39.7 x
>   1e-10 cm^3 s^-1 at 300, 600, 1200 K), log-log interpolated and held above
>   1200 K; Table 4's 4.0e-9 is the 1200 K entry used as a constant. The
>   product state "could not be determined"; their near-resonance argument
>   puts it in Fe+ 3d6(5D)4p z4F/z4D, which decay by allowed lines, so **A6 is
>   not carried** and each event heats the gas by 1.18 eV (not 5.70).
> - **A19 is Zhao et al. (2005) Table III** (READ; S(3P) and S(1D) totals,
>   20 K - 1e6 K) weighted by the populations within Q(S I); **A20 is
>   computed from A19 by detailed balance**.
> - **A21 is carried again**, as the sum of the non-radiative Na(3s) + H+ ->
>   Na+ + H(n=2) of Dutta et al. (2001) (the Maxwellian average of their Fig. 5
>   cross sections; their Fig. 7 rates break the threshold law below about
>   5000 K and are not used) and the radiative channel of Watanabe et al.
>   (2002) Fig. 5. **A22 is not carried** (the non-radiative product is
>   H(n=2)). **A23 is carried again**, Watanabe et al.'s radiative K + H+.
>   None of the three papers gives a fit; Table 4's exp[sum c_i (ln T)^i]
>   forms appear in neither.]

## Group B: charge exchange with He / He+ (only when `cx_full` is on)

> [2026-08-15: the heading is out of date, B1 and B2 are on by default under
> `He_H_charge_exchange`, not gated by `cx_full`.]

| # | Reaction | Ref | Rate (cm^3 s^-1) |
| --- | --- | --- | --- |
| B1 | He + H+ | GJ07 | `1.75e-11 * (T/300)^(-0.75) * exp(-12.75/T4)` |
| B2 | He+ + H | GJ07 | `1.25e-15 * (T/300)^0.25` |

> [2026-09-24: B1 is no longer the Glover & Jappsen (2007) R27 fit of Kimura
> et al. (1993). It is computed as the detailed-balance reverse of the
> non-radiative He+ + H channel of Kingdon & Ferland (1996) (from Zygelman
> et al. 1989), 4 exp(-10.989 eV/kT) times that rate; the two calculations
> disagree by 25x at 1e4 K. B2 is unchanged and is radiative.]
>
> [2026-09-26: the non-radiative He+ + H channel behind B1 is now the
> Maxwellian average of the H+ + He cross section of Loreau, Ryabchenko &
> Vaeck (2014) turned forward by reciprocity (4.4e-17 at 1e4 K, against
> 9.0e-15 from Kingdon & Ferland, which fits Zygelman et al.'s values without
> the 1/4 of the singlet approach), and B1 is its detailed-balance reverse.
> B2 is 1.20e-15 (T/300)^0.25, the value Stancil, Lepp & Dalgarno (1998)
> Table 1 row (19) print; Glover & Jappsen's 1.25e-15 was carried until then.
> See section 12.9 of the detailed-balance memo.]

## Group C: metal + He / He+ (only when `cx_full` is on)

Let `f_Si(T) = 3.32e-13*sqrt(T) + 1.2e-16*T + 4.2e-9/sqrt(T) - 7.9e-13`.

| # | Reaction | Ref | Rate (cm^3 s^-1) |
| --- | --- | --- | --- |
| C1 | Si + He+ | Satta13 | `f_Si(T)` |
| C2 | Si+ + He | Satta13 | `1.415 * f_Si(T) * T^0.103 * exp(-19.1/T4)` |
| C3 | C + He+ | GJ07 | `2.5e-15 * (T/300)^1.597` |
| C4 | C+ + He | GJ07 | `6.75e-15 * (T/300)^1.654 * exp(-15.5/T4)` |
| C5 | O + He+ | Z04 | `4.99e-15 * T4^0.379 + 2.78e-15 * T4^(-0.216) * exp(T4/81.97)` |
| C6 | O+ + He | Z04 | `3.2 * (5.0e-15*T4^0.38 + 2.78e-15*T4^(-0.22)*exp(T4/81.97)) * T^0.0377 * exp(-12.7/T4)` |

> [2026-09-24: C1/C2 are carried as printed (Satta et al. 2013 not read;
> superseded 2026-09-25, block below).
> C3 is the three-branch fit of Glover & Jappsen (2007) R41 to Kimura et al.
> (1993, ApJ 417, 812) Table 1 (READ); Table 4 carries only the T > 2000 K
> branch. **C4 is not carried**: Kimura et al.'s channel ends in C+(2D),
> which decays by an allowed line. C5 is Zhao et al. (2004, ApJ 615, 1063)
> eq. (19) with their Table 3 parameters (READ):
> `4.991e-15*T4^0.3794*exp(-T/1.121e6) + 2.780e-15*T4^(-0.2163)*exp(T/8.158e5)`;
> Table 4 drops the first exp and prints 81.97 for 81.58. The channel is
> RADIATIVE ("Total Rate Coefficients for Radiative Charge Transfer"), so
> **C6 is not carried**.
> Second round: C5 is evaluated with T held inside the calculated 10 - 1e6 K,
> since exp(+T/8.158e5) grows without bound above it.]

## Group D: metal + metal (only when `cx_full` is on)

Let `D_CSi(T) = 1.87e8 + 5.09e10 * T^(-0.527)`.

| # | Reaction | Ref | Rate (cm^3 s^-1) |
| --- | --- | --- | --- |
| D1 | C + Si+ | Satta13 | `(0.724 * T^0.0463 * exp(-3.61/T4)) / D_CSi(T)` |
| D2 | C+ + Si | Satta13 | `1 / D_CSi(T)` |
| D3 | Mg + Si+ | PH80 | `2.9e-9` |
| D4 | Mg+ + Si | PH80 | `9.8e-10 * (T/300)^(-0.0264) * exp(-0.59/T4)` |
| D5 | Fe + Si+ | PH80 | `1.9e-9` |
| D6 | Fe+ + Si | PH80 | `1.4e-9 * (T/300)^(-0.236) * exp(-0.299/T4)` |
| D7 | Na + Mg+ | PH80 | `1.0e-11` |
| D8 | Na+ + Mg | PH80 | `3.76e-11 * (T/300)^(-0.0302) * exp(-2.909/T4)` |
| D9 | Mg + C+ | PH80 | `1.1e-9` |
| D10 | Mg+ + C | PH80 | `3.01e-10 * (T/300)^(-0.0832) * exp(-4.19/T4)` |
| D11 | Mg + N+ | PH80 | `1.2e-9` |
| D12 | Mg+ + N | PH80 | `9.48e-10 * (T/300)^0.139 * exp(-7.99/T4)` |
| D13 | Na + Fe+ | PH80 | `1.0e-11` |
| D14 | Na+ + Fe | PH80 | `2.4e-11 * (T/300)^(-0.0987) * exp(-3.21/T4)` |
| D15 | Fe + C+ | PH80 | `2.6e-9` |
| D16 | Fe+ + C | PH80 | `1.12e-9 * (T/300)^0.0147 * exp(-3.90/T4)` |
| D17 | Fe + O+ | RV72 | `1.71e-9` |
| D18 | Fe+ + O | RV72 | `5.0e-10 * (T/300)^(-0.0337) * exp(-6.63/T4)` |
| D19 | Mg + S+ | PH80 | `2.8e-10` |
| D20 | Mg+ + S | PH80 | `5.35e-11 * (T/300)^0.121 * exp(-3.15/T4)` |
| D21 | Na + S+ | PH80 | `2.6e-10` |
| D22 | Na+ + S | PH80 | `1.86e-10 * (T/300)^0.151 * exp(-6.06/T4)` |
| D23 | Fe + S+ | PH80 | `1.8e-10` |
| D24 | Fe+ + S | PH80 | `5.38e-11 * (T/300)^0.052 * exp(-2.85/T4)` |
| D25 | S + C+ | Chenel10 | `3.26e-11 * (T/300)^0.289 * exp(-338/T)` |
| D26 | S+ + C | Chenel10 | `6.97e-11 * (T/300)^0.356 * exp(-1.06/T4)` |
| D27 | Fe + N+ | RV72 | `9.46e-10` |
| D28 | Fe+ + N | RV72 | `1.17e-9 * (T/300)^0.07 * exp(-7.696/T4)` |
| D29 | Na + Ca+ | Kwolek19 | `3.0e-9` |
| D30 | Na+ + Ca | Kwolek19 | `1.7e-8 * (T/300)^(-0.168) * exp(-1.13/T4)` |
| D31 | K + Na+ | Lavvas14 | `1.0e-11` |
| D32 | K+ + Na | Lavvas14 | `6.58e-12 * (T/300)^(-0.203) * exp(-0.927/T4)` |

> [2026-09-24: the Prasad & Huntress (1980) forwards D3, D5, D7, D9, D11, D13,
> D15, D19, D21 and D23 were checked against their Table 2 (READ; reactions
> 243, 241, 274, 138, 188, 273, 136, 254, 253, 252, all with alpha = beta =
> 0), and D31 against Lavvas et al. (2014) Table 2 R16 (READ, "Est"). Their
> reverses D4, D6, D8, D10, D12, D14, D16, D20, D22, D24 and D32 are now
> **computed by detailed balance** with
> internal partition functions; the printed power-law fits are not carried.
> D1/D2, D17/D18 and D27/D28 are carried as printed: Satta et al. (2013) and
> Rutherford & Vroom (1972) were not read (superseded 2026-09-25, block
> below). Second round: **D26 and D30 are
> computed by detailed balance** from D25 (Chenel et al. 2010 computed C+ + S)
> and D29 (Kwolek et al. 2019 measured Na + Ca+), both transcribed as Table 4
> prints them; their 0.90 and 0.97 eV defects reach no excited product term,
> so the reverse of the ground-term channel is the detailed-balance one.
> Table 4's own microscopic-balance fits departed from it by 0.25 - 0.52
> (D26) and 1.2 - 4.1 (D30) over 5e3 - 2e4 K. [D30: superseded 2026-09-25, not carried; Kwolek et al.'s
> channel ends in Ca 4s4p 3P.] Prasad & Huntress also list
> C+ + Na (137, 1.1e-9), Si+ + Na (242, 2.7e-9) and S+ + Si (251, 1.6e-9),
> which Table 4 does not carry.]

> [2026-09-25, third round: **C1 is Satta et al. (2013) Table 1 k1** (READ:
> 3.226 28(-13), 1.327 26(-16), 4.1533(-9) (printed "4.1.533(-9)"),
> -9.162 43(-13) in k = c0 T^0.5 + c1 T + c2 T^-0.5 + c3; Table 4's
> 3.32e-13, 1.2e-16, 4.2e-9, -7.9e-13 differ); the product is an excited
> Si+* that radiates, so **C2 is not carried**. Satta et al. contain no C +
> Si+ rate: **D1/D2 stay blocked**, their origin unknown. **D17 and D27 are
> Rutherford & Vroom Table I** (O+: 29.4, 22.2, 17.1; N+: 14.7, 11.4, 9.46 x
> 1e-10 at 300, 600, 1200 K); **D18 and D28 are not carried** (near-resonant
> excited Fe+ products). **D25 is Chenel et al. (2010) Table 2, wave-packet
> column** (1.3e-11 - 2.2e-10 over 500 - 1e5 K); D26 stays detailed balance.
> **D29 = 3.0e-9 is Kwolek et al.'s Langevin rate kL(g) = 2.995e-9**, which
> their trap measurement (relative collision energies about 0.015 - 0.09 eV)
> reaches; they infer the product Ca 4s4p 3P, so **D30 is not carried**.]

## Beyond Table 4: group E and He2+ + H (2026-09-26)

> [2026-09-26: two reactions Table 4 does not carry were added; the details
> and the checks are in `md/charge_exchange_detailed_balance_20260924.md`,
> section 12.]
>
> | # | Reaction | Source | Rate (cm^3 s^-1) | Product | Reverse |
> | --- | --- | --- | --- | --- | --- |
> | E1 | O2+ + H -> O+ + H+ | Barragan et al. (2006) Table 1 k1 | table, 1e2 - 1e5 K | O+ 2s2p4 4P (833 A decay) | not carried |
> | E2 | N2+ + H -> N+ + H+ | Barragan et al. (2006) Table 2 k4 | table, 1e2 - 1e5 K; 9.70e-10 at 1e4 K | N+ 2s2p3 3Do (1084 A decay; Barragan et al.'s text, Butler, Heil & Dalgarno 1980 Table 2) | not carried (the product decays by an allowed line; 16.0 eV endothermic) |
> | B3 | He2+ + H -> He+(1s) + H+ + photon | West, Lane & Cohen (1982), Maxwellian average of their radiative cross section (Figs. 4-6) | table, 1e2 - 1e5 K; 1.65e-13 at 1e4 K (Kingdon & Ferland's 1.0e-14 until 2026-09-26) | He+(1s) and a photon (radiative) | not carried |
>
> E1 and E2 are gated by `metals.inp: cx_O2p_H` and `cx_N2p_H` (scale, default
> 1, 0 leaves the row out); B3 is applied with the He <-> H pair under
> `He_H_charge_exchange`, in every system that carries He III.
>
> [2026-09-26: the detailed-balance reverses of the pairs whose sources name
> the product levels (A17 from A18, A20 from A19, D26 from D25) now depend on
> the electron density through the populations of N I, S II and C I; the
> other derived reverses do not (thermal branching). See section 12 of the
> detailed-balance memo.]

## Count

63 distinct rows transcribed above (23 in A, 2 in B, 6 in C, 32 in D). The paper
states "a total of 65 charge exchange reactions"; the two-row difference is in
the metal-metal/He set (Groups C-D) and is immaterial to the H ionization
balance and to WASP-121b (see below).

> [2026-09-26: the code carries 65 rows: the 63 above (ten of them not
> applied, sections 10-11 of the detailed-balance memo) and E1, E2 of group
> E; He2+ + H (B3) is applied with the He <-> H pair outside the row table.]

## Physical context for WASP-121b (Case A)

Huang Section 4.2 states explicitly that on WASP-121b "the charge exchange rate
between Fe I and H+ is lower than the H photoionization and recombination rates,"
and that the metal mechanisms (including Fe-H charge exchange) "do not play a
significant role on WASP-121b" because the atmosphere is highly ionized. Charge
exchange is expected to matter more for cooler, lower-XUV planets. The Phase-1d
validation therefore reports what the model shows (CX rates vs recombination,
Huang Fig. 11) rather than asserting CX dominance.

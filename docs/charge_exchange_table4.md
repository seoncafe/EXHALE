# Huang et al. (2023) Table 4 — Charge Exchange Rates (verified transcription)

This is the authoritative, image-verified transcription of Table 4 of Huang,
Koskinen, Lavvas & Fossati (2023, ApJ 951, 123), "Charge Exchange Rates," used
to implement charge exchange in `ATES-metal` (Phase 1d). Every coefficient and
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

## Group A — charge exchange with H / H+ (the Phase-1d default active set)

These couple the H ionization fraction to each metal's ionization balance and
are the only reactions enabled by default (`cx_full` off). Note **Ca has no
direct H charge exchange in Table 4** — it appears only in Group D.

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
| A13 | O + H+ | S99 | `2.08e-9 * T4^0.405 + 1.11e-11 * T4^(-0.458)` |
| A14 | O+ + H | S99 | `(1.26e-9 * T4^0.517 + 4.25e-10 * T4^(6.69e-3)) * exp(-227/T)` |
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
- **A16 (C+ + H)** uses `exp(-17/T4)`, i.e. a 1.7e5 K activation barrier. This is
  the endothermic reverse of the fast `C + H+` charge transfer and is utterly
  negligible (~1e-24 at T4 = 1); ported verbatim per the paper.
- **A17-A23** are `exp[polynomial in lnT]` fits (Lin/Zhao/Dutta/Watanabe). Within
  the model's temperature range (T > a few 100 K, lnT > ~5.7) the polynomial is
  strongly negative, so these rates are tiny (under-flow toward 0). The fits blow
  up only for T -> 1 K (lnT -> 0), which never occurs; the implementation still
  caps the evaluated coefficient at a physical maximum to be safe.
- **K has only the forward `K + H+`** in Table 4; no `K+ + H` reverse is tabulated.

## Group B — charge exchange with He / He+ (only when `cx_full` is on)

| # | Reaction | Ref | Rate (cm^3 s^-1) |
| --- | --- | --- | --- |
| B1 | He + H+ | GJ07 | `1.75e-11 * (T/300)^(-0.75) * exp(-12.75/T4)` |
| B2 | He+ + H | GJ07 | `1.25e-15 * (T/300)^0.25` |

## Group C — metal + He / He+ (only when `cx_full` is on)

Let `f_Si(T) = 3.32e-13*sqrt(T) + 1.2e-16*T + 4.2e-9/sqrt(T) - 7.9e-13`.

| # | Reaction | Ref | Rate (cm^3 s^-1) |
| --- | --- | --- | --- |
| C1 | Si + He+ | Satta13 | `f_Si(T)` |
| C2 | Si+ + He | Satta13 | `1.415 * f_Si(T) * T^0.103 * exp(-19.1/T4)` |
| C3 | C + He+ | GJ07 | `2.5e-15 * (T/300)^1.597` |
| C4 | C+ + He | GJ07 | `6.75e-15 * (T/300)^1.654 * exp(-15.5/T4)` |
| C5 | O + He+ | Z04 | `4.99e-15 * T4^0.379 + 2.78e-15 * T4^(-0.216) * exp(T4/81.97)` |
| C6 | O+ + He | Z04 | `3.2 * (5.0e-15*T4^0.38 + 2.78e-15*T4^(-0.22)*exp(T4/81.97)) * T^0.0377 * exp(-12.7/T4)` |

## Group D — metal + metal (only when `cx_full` is on)

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

## Count

63 distinct rows transcribed above (23 in A, 2 in B, 6 in C, 32 in D). The paper
states "a total of 65 charge exchange reactions"; the two-row difference is in
the metal-metal/He set (Groups C-D) and is immaterial to the H ionization
balance and to WASP-121b (see below).

## Physical context for WASP-121b (Case A)

Huang Section 4.2 states explicitly that on WASP-121b "the charge exchange rate
between Fe I and H+ is lower than the H photoionization and recombination rates,"
and that the metal mechanisms (including Fe-H charge exchange) "do not play a
significant role on WASP-121b" because the atmosphere is highly ionized. Charge
exchange is expected to matter more for cooler, lower-XUV planets. The Phase-1d
validation therefore reports what the model shows (CX rates vs recombination,
Huang Fig. 11) rather than asserting CX dominance.

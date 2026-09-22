# Line trapping in the thick metal resonance lines: measurement and verdict

2026-08-11. Earlier work introduced the escape probability `beta(tau)` for the two ground-term
fine-structure lines `[O I] 63um` and `[C II] 158um` and explicitly left the
thick resonance lines of a metal-rich wind untreated. This memo measures that
open item and closes it. (Those two lines were the whole trapped set at the
time; it is now the eight ground-term lines of C I, C II, N II and O I, which
does not affect the verdict below: that verdict is about permitted resonance
lines, a different class.)

**Verdict. The correction is a no-op and was not implemented.** The lines are
genuinely optically thick (Mg II h&k reaches `tau0 ~ 7.6e4` at the WASP-121 b
base and stays above `1e2` through the region where Mg II carries a third of
the local cooling), but trapping does not suppress their cooling, because a
trapped photon in a permitted resonance line still escapes long before a
collision can destroy the excitation. Applying the correct correction changes
the total radiative losses of WASP-121 b by **0.02%**, of HD 209458 b by
0.002% and of HD 189733 b by 0.0015%. That is a physical result, not a
tolerance argument: the suppression factor is `1` to within the measurement
because the relevant density threshold sits one to six decades above the
densities these winds reach. The finding is recorded as a validity note in
`Cool_coeff.f90` next to the `beta = 1` it justifies.

---

## 1. Why `tau0 > 1` is the wrong test

The intuition that a thick line must radiate less is imported from the
fine-structure case, where it is right for a reason that does not carry over.

Write the two-level statistical equilibrium of a line with an escape
probability `beta`, collisional excitation `C_lu = ne q_lu` and de-excitation
`C_ul = ne q_ul`:

```
n_l C_lu = n_u ( beta A_ul + C_ul )
Lambda   = h nu * n_u * beta A_ul
         = h nu * n_l C_lu * [ beta A_ul / ( beta A_ul + ne q_ul ) ]
```

The CHIANTI closed-form fits in `Cool_coeff.f90` (`cool_MgI_func`,
`cool_MgII_func`, `cool_CaII_func`, `cool_NaI_func`, and the resonance group of
the Fe II fit) are by construction the coronal limit `Lambda = h nu q_lu`, i.e.
every collisional excitation ends as an escaped photon. The correction to them
is therefore the bracket

```
S = beta A_ul / ( beta A_ul + ne q_ul ),    q_ul = 8.629e-6 Ups(T) / ( g_u sqrt(T) )
```

and **not** a factor `beta`. Multiplying by `beta` would be the answer for a
line whose photons are destroyed on absorption; a resonance line scatters, and
in the pure-scattering limit the energy is emitted anyway, just later and from
further out. `S` departs from 1 only when the ion is collisionally de-excited
faster than the trapped photon diffuses out, i.e. above the effective critical
density

```
n_crit,eff = beta A_ul / q_ul .
```

This is the same structure the code already uses for `[O I]` and `[C II]`,
where `beta` enters as `A_ul -> beta A_ul` *inside* `lambda_2level` rather than
multiplying its result. The present memo is that same statement carried to the
resonance lines; what differs is only the value of `A_ul`.

* forbidden fine structure: `A_ul ~ 1e-5 - 1e-3 s^-1`, `n_crit,eff ~ 1e0 - 1e5 cm^-3`,
  far below the base density, so `S` is small and `beta` matters;
* permitted resonance: `A_ul ~ 1e8 s^-1`, and even `beta ~ 1e-4` leaves
  `beta A_ul ~ 1e4 s^-1`, so `n_crit,eff ~ 1e10 - 1e15 cm^-3`.

The wind never gets near the second threshold. That, and not the size of
`tau0`, is what decides the question.

## 2. What was measured

Read-only from `WASP-121b/output/` (the Newton-converged production run,
`Ion_species.txt`, `Hydro_ioniz.txt`, `Cooling_breakdown.txt`), and the same
for `HD209458b/output/` and `HD189733b/output/`.

For each line, the line-center opacity uses the same expression as
`line_center_opacity_lte`, in the equivalent oscillator-strength form
`kappa_0 = (sqrt(pi) e^2/m_e c) f_lu n_l / Dnu_D` with a Doppler core
`v_th = sqrt(2kT/m)` and no turbulent broadening: an upper bound on
`kappa_0`, so a *lower* bound on `beta`. Two escape channels are combined
exactly as `lya_rt.f90` does, `beta = 1 - (1-beta_static)(1-beta_Sobolev)`:

* `beta_static` from the outward line-center column (cell center to the top of
  the domain), through `line_escape_probability_one_face`;
* `beta_Sobolev` from the radial velocity gradient,
  `tau_S = (pi e^2/m_e c) f_lu lambda n_l / |dv/dr|`.

Lower-level populations: the ground level of Mg I/Mg II/Ca II/Na I holds
essentially the whole ion (the upper terms lie at 2.1-4.3 eV). For Fe II the
`a6D 9/2` fraction is a Boltzmann sum over the 16 lowest NIST levels
(`a6D`, `a4F`, `a4D`, `a4P`, `b4P`), giving ~0.3 at `1e4 K`.

Collision strengths `Ups(T)` are the ones the cooling fits themselves carry,
split between the two doublet components in the ratio of their `f`-values. Fe II
has no single `Ups`; its UV resonance group is recovered from the fitted term
`6.536e-16 exp(-58939/T)/sqrt(T)` of `cool_FeII_func`, which gives an effective
`Ups ~ 93` for `g_l = 10`. That lumps hundreds of transitions onto one line and
so over-states `tau0` per line; it is a deliberately conservative proxy.

## 3. Optical depths reached

WASP-121 b, `static | Sobolev` line-center depth:

| `r/R_p` | `T` [K] | Mg II k 2796 | Ca II K 3934 | Na I D2 5890 | Fe II 2383 | Mg I 2853 |
|---|---|---|---|---|---|---|
| 1.000 | 2358 | 7.5e4 \| 2.7e9 | 7.5e3 \| 2.1e8 | 3.6e3 \| 1.3e8 | 1.7e4 \| 3.6e8 | 3.9e4 \| 1.2e9 |
| 1.050 | 3259 | 3.9e3 \| 2.8e4 | 3.0e2 \| 1.8e3 | 5.3e1 \| 8.5e2 | 4.8e2 \| 2.9e3 | 4.4e2 \| 4.2e3 |
| 1.100 | 6739 | 1.2e3 \| 3.6e3 | 8.2e1 \| 2.2e2 | 5.1e-1 \| 6.3e0 | 1.1e2 \| 2.8e2 | 1.8e2 \| 3.7e2 |
| 1.200 | 9759 | 3.1e2 \| 2.1e2 | 1.8e1 \| 1.1e1 | 2.4e-2 \| 1.9e-2 | 1.6e1 \| 1.3e1 | 5.5e1 \| 3.8e1 |
| 1.300 | 10283 | 1.4e2 \| 2.0e2 | 7.1e0 \| 8.9e0 | 1.1e-2 \| 1.6e-2 | 4.6e0 \| 6.2e0 | 2.6e1 \| 3.6e1 |
| 1.501 | 11111 | 2.6e1 \| 6.1e1 | 1.1e0 \| 2.1e0 | 1.9e-3 \| 4.7e-3 | 5.3e-1 \| 9.2e-1 | 4.5e0 \| 1.1e1 |

and the share each ion holds of the *local* total cooling:

| `r/R_p` | `T` [K] | Mg II | Ca II | Na I | Fe II | O I | C I | Ly-alpha |
|---|---|---|---|---|---|---|---|---|
| 1.000 | 2358 | 6.4e-5 | 6.4e-4 | 0.069 | 8.9e-4 | 0.265 | 0.619 | ~0 |
| 1.050 | 3259 | 6.1e-3 | 0.011 | 0.204 | 2.1e-3 | 0.246 | 0.259 | ~0 |
| 1.100 | 6739 | 0.292 | 0.065 | 9.8e-3 | 0.074 | 0.227 | 0.078 | 0.026 |
| 1.200 | 9759 | 0.277 | 0.032 | 1.4e-4 | 0.088 | 0.021 | 0.027 | 0.215 |
| 1.300 | 10283 | 0.343 | 0.030 | 1.4e-4 | 0.058 | 3.0e-3 | 0.012 | 0.086 |
| 1.501 | 11111 | 0.349 | 0.022 | 1.1e-4 | 0.032 | 6.2e-4 | 3.7e-3 | 0.067 |

So the literal form of the test (is the ion a significant coolant where the
line is thick?) is passed decisively. Over 99% of the integrated Mg II, Ca II,
Fe II and Mg I cooling, and 73-81% of the Na I cooling, comes from gas with
`tau0 > 1`; Mg II peaks at 35% of the local cooling at `tau0 ~ 1e2-1e3`. On
that criterion alone one would implement the correction.

## 4. The suppression factor, and why it is 1

`S = beta A / (beta A + ne q_ul)` on the same WASP-121 b profile:

| `r/R_p` | `ne` [cm^-3] | `beta`(Mg II k) | `beta A / ne q_ul` | `S`(Mg II k) | `S`(Ca II K) | `S`(Na I D2) | `S`(Fe II 2383) |
|---|---|---|---|---|---|---|---|
| 1.000 | 1.9e9 | 2.0e-6 | 0.87 | 0.462 | 0.852 | 0.698 | 0.534 |
| 1.020 | 1.3e9 | 8.4e-6 | 4.3 | 0.810 | 0.970 | 0.940 | - |
| 1.050 | 1.4e9 | 8.2e-5 | 56 | 0.982 | 0.998 | 0.998 | - |
| 1.100 | 2.0e9 | 4.5e-4 | 2.9e2 | 0.997 | 1.000 | 1.000 | - |
| 1.200 | 4.3e9 | 5.4e-3 | 1.9e3 | 0.999 | 1.000 | 1.000 | - |
| 1.300 | 3.0e9 | 6.6e-3 | 3.4e3 | 1.000 | 1.000 | 1.000 | - |
| 1.501 | 1.5e9 | 2.8e-2 | 2.9e4 | 1.000 | 1.000 | 1.000 | - |

`beta` falls to `2e-6` at the base (the trapping is extreme) and `S` still
only reaches 0.46 there, because `beta A_ul = 5e2 s^-1` remains comparable to
`ne q_ul`. Everywhere the lines actually cool, `beta A_ul` exceeds `ne q_ul` by
two to five decades.

Equivalently, in terms of the threshold density (WASP-121 b, `cm^-3`):

| `r/R_p` | `ne` | `n_crit,eff` Mg II k | Ca II K | Na I D2 | Mg I 2853 |
|---|---|---|---|---|---|
| 1.000 | 1.9e9 | 1.6e9 | 1.1e10 | 4.3e9 | 4.7e10 |
| 1.100 | 2.0e9 | 5.7e11 | 5.4e12 | 9.5e13 | 4.0e13 |
| 1.200 | 4.3e9 | 8.2e12 | 8.2e13 | 1.8e14 | 2.9e14 |
| 1.501 | 1.5e9 | 4.5e13 | 5.1e14 | 1.9e14 | 1.6e15 |

`ne` crosses `n_crit,eff` only in the first cell of the densest of the three
planets, and only for Mg II.

## 5. Integrated impact

`Delta Lambda = Lambda (1 - S)`, integrated as `4 pi r^2 dr`:

| line | WASP-121 b | HD 209458 b | HD 189733 b |
|---|---|---|---|
| Mg II k 2796 | 0.076% of Mg II | 0.001% | 0.007% |
| Mg II h 2804 | 0.019% | 0.000% | 0.002% |
| Ca II K 3934 | 0.012% of Ca II | 0.000% | 0.001% |
| Ca II H 3968 | 0.003% | 0.000% | 0.000% |
| Na I D2 5890 | 1.06% of Na I | 0.150% | 0.030% |
| Na I D1 5896 | 0.32% | 0.064% | 0.007% |
| Mg I 2853 | 0.002% of Mg I | 0.000% | 0.000% |
| Fe II 2383 | 0.022% of Fe II | 0.077% | 0.758% |
| **total cooling** | **0.0214%** | **0.0022%** | **0.0015%** |

The largest *local* effect is 1.65% of the cooling of the first cell of
WASP-121 b; the loss exceeds 0.1% in 36 of 504 cells, all below `r = 1.007
R_p`, and exceeds 0.01% nowhere above `r = 1.27 R_p`. The first two cells sit
at `T = T_eq` exactly, pinned by the lower boundary condition, so even that
1.65% has little room to feed back on the temperature structure.

Na I D is the one line where the effect is not utterly negligible in relative
terms (1.4% of the Na I cooling), because it has the lowest `A_ul` of the set
(`6.2e7 s^-1`) and its cooling is concentrated at `r = 1.04-1.14` where the gas
is coldest and most neutral. Na I is 0.3% of the WASP-121 b radiative losses,
so this does not survive into the total.

## 6. Why the numbers above are a conservative bound

Three choices in the measurement all push `S` down, i.e. exaggerate the effect
that was found to be absent:

1. **Doppler-only escape.** `line_escape_probability_one_face` is the
   plane-parallel Doppler-core form of de Jong, Boland & Dalgarno (1980)
   eq. (B-7). These resonance lines are in the damping wings at the depths
   of section 3, where the published static-slab solution is Harrington
   (1973) eq. (40), `beta = 1/(1 + 0.909 tau0)` (see section 7). That is
   8x the Doppler form at `tau0 = 1e2`, 10x at `1e3` and 13x at the
   `7.5e4` of the WASP-121 b base, hence `S` closer to 1. The Doppler form
   is the right one for the fine-structure lines it was written for
   (`a ~ 1e-8`), not for these.
2. **No turbulent or bulk broadening in `kappa_0`**, which maximizes `tau0`.
3. **Electron collisions only in `q_ul`.** Adding neutral-H de-excitation would
   raise `ne q_ul` at the base. This is the one choice that goes the other way,
   and it is the reason the base cells are quoted as an upper bound rather than
   a measurement: quantum de-excitation rates for Na I D by H (Barklem et al.)
   are of order `1e-10 cm^3 s^-1` at 2000-3000 K, against `n_HI ~ 3e12`, so
   `n_HI k_ul ~ 3e2 s^-1`, comparable to `beta A_ul` there. It would deepen
   the ~1% first-cell effect and would not touch the wind, where H is ionized.
4. **The Fe II proxy** puts the whole UV resonance group on a single line,
   over-stating `tau0` per line.

A fourth destruction channel exists in principle and was not measured:
photoionization of the *upper* level while the photon is trapped. For Mg II the
`3p` level ionizes at 10.75 eV (1153 A), outside the 13.6 eV-1.24 keV bands the
model carries, so the code cannot supply the rate. An order-of-magnitude
estimate with `sigma ~ 1e-18 cm^2` and the FUV photon flux at 0.025 AU gives
`~1e-2 s^-1` against `beta A_ul ~ 1e5 s^-1`; on that estimate it appears
negligible, but it is an estimate, not a measurement.

## 7. Two related lines checked in passing

* **Fe II ground-term fine structure** (`a6D`, 25.99 and 35.35 um, `A ~ 2e-3
  s^-1`) is the one place in the metal set where an `A_ul` small enough to trap
  meets a real population. It is optically thin: `tau0 <= 0.06` in all three
  runs, because Fe II peaks near `1e5-1e6 cm^-3` and the `lambda^3` gain does
  not make up for it. The `cool_FeII_ne` 2-D statistical-equilibrium table,
  which does carry these lines and their LTE saturation, therefore needs no
  escape probability.
* **H I Ly-alpha** is not in scope but is the extreme case of the same physics
  and was evaluated for context: `tau0 = 1.3e8` at the WASP-121 b base,
  `beta ~ 2e-6` on the one-flight wing form `pi^(-1/4) sqrt(a/tau)` that
  `lya_rt.f90` uses, and yet `beta A / ne q_ul = 1.1e2` at the base rising to
  `1.8e5` in the wind, giving a cooling-weighted `S = 0.9994`. Whatever
  suppresses Ly-alpha cooling in a planetary wind, it is not two-level
  trapping; the candidate channels are destruction of the trapped photon
  (photoionization of `H(n=2)`, collisional `2p -> 2s` transfer to two-photon
  decay), which is a different calculation and is not made here.
  `Cool_coeff.f90` continues to treat the Ly-alpha collisional-excitation
  cooling as thin.

  Two corrections to that paragraph, from items REF-NEUFELD and LYA-BETA
  (2026-09-07). The one-flight wing form was **not** the Neufeld (1990) or
  Harrington (1973) result and was not the escape probability of a Ly-alpha
  photon; `lya_rt.f90` no longer uses it. The published static-slab,
  damping-wing solution is Harrington (1973), MNRAS 162, 43, eq. (40),
  `<N> = (4 sqrt(6)/pi^2) u_2 B = 0.909316 B` scatterings for a mid-plane
  source in a slab of line-centre optical half-thickness `B`, generalized to a
  source plane anywhere between the two faces by Neufeld (1990), ApJ 350, 216,
  eq. (3.27) at zero continuum destruction. `<N>` counts absorptions, so the
  escape chance per emission is `beta = 1/(1 + <N>)`, which is what the code
  now carries: proportional to `1/tau` in the wings, independent of the Voigt
  parameter, and tending to 1 in the thin limit. On the stored `wasp_full`
  column that is `9.1e-9` at the base against the `2.1e-6` used here, a factor
  234. With it, `beta A / ne q_ul = 0.50` at the base and `S = 0.33` rather
  than `~1`, so the base cells **are** trapping-suppressed; the wind is not
  (`beta A / ne q_ul = 70` at `1.20 R_p`, `4.4e2` at `1.28 R_p`), and the
  cooling-weighted conclusion of this section, which is set by the wind,
  stands.

## 8. What was changed, and what was not

Changed: the `SCOPE` comment of the fine-structure trapping block in
`src/modules/radiation/Cool_coeff.f90` now states the validity condition
`ne << beta A_ul / q_ul` and its measured margin, in place of the sentence
recording the resonance lines as untreated.

**No code path was modified and no result changed.** `beta = 1` for the
resonance lines is kept, now as a justified effective treatment rather than an
acknowledged omission.

Not checked: any planet other than the three above; a base an order of
magnitude denser in `ne`, where the Mg II threshold would be crossed over more
than the first cell; the neutral-H de-excitation channel of §6.3 as a rate
rather than an estimate; the effect of the upper-level pile-up `n_u = n_l C_lu
/ (beta A + C_ul)` on the ionization balance and on the transmission spectrum,
which is a separate question from cooling and is untouched by this memo.

## 9. Reproduction

`docs/` carries no script for this; the measurement is a direct read of the
three `output/` directories with the formulas of §1-§2 and takes a few seconds.
The inputs are `Ion_species.txt` (columns 8-34, canonical metal-ion order),
`Hydro_ioniz.txt` (`r`, `v`, `T`) and `Cooling_breakdown.txt` (`ne`,
`cool_total`, and one column per metal ion from column 12 --- column 11 is
`H3p_IR` since 2026-08-13; read the `# col...` header rather than the fixed
offset). Atomic data:
NIST ASD for `lambda`, `f_lu`, `A_ul`, `g_u` and the Fe II level list; `Ups(T)`
from the fits in `Cool_coeff.f90` itself.

---

## 10. The escape probability of the fine-structure lines: source and argument

2026-09-07, item LYA-BETA-METALS. Sections 1-9 settle the *resonance* lines
(`beta = 1`, justified). This section settles the function that carries the
`beta` of the other class, the eight ground-term fine-structure lines of C I,
C II, N II and O I, which are the only lines in the code that get one.

### 10.1 Where it is used

`line_escape_probability_one_face` in `Cool_coeff.f90` is called from exactly
one place, `fine_structure_line_transfer`, which fills `beta_fs` and
`nbar_fs` for the eight lines `[C I] 609/370um`, `[C II] 158um`,
`[N II] 205/122um` and `[O I] 63/145/44um`. They enter the cooling as
`A_ul -> beta A_ul` inside the statistical equilibrium of `cool_CI_ne_func`,
`cool_CII_ne_func`, `cool_NII_ne_func` and `cool_OI_ne_func` (and the
two-level legacy branch of `util_ion_eq.f90`). With `Base IR field` on, the
inward face value also sets the incident field `nbar_fs`.

No metal resonance line uses it: `cool_MgI_func`, `cool_MgII_func`,
`cool_CaII_func`, `cool_NaI_func`, the Fe I table and the Fe II 2-D
statistical-equilibrium table are all evaluated at `beta = 1`, for the reason
of sections 1-5. The transmission post-processor (`EXHALE_transit.py`,
`exhale_transit_lib.py`) does not use an escape probability at all; it
integrates Voigt profiles directly.

### 10.2 The published source

The expression **is** published, and the attribution can be given to the
equation. It is **de Jong, Boland & Dalgarno (1980), A&A 91, 68, appendix B,
eq. (B-7)**:

```
beta(tau) = [1 - exp(-2.34 tau)] / (4.68 tau)          for tau < 7
          = 1 / (4 tau [ln(tau/sqrt(pi))]^(1/2))       for tau >= 7
```

quoted there as approximating `beta(tau) = (1/2) int dx phi(x) E_2[tau phi(x)]`
"when photons escape only through the nearest cloud boundary", and "accurate
to within 10% for small and intermediate `tau` and exact at very large
`tau`". The code reproduces both branches to round-off, takes `beta(0) = 1/2`
as (B-7) does, and switches branches at `tau_c = sqrt(pi) exp(2.34^2/4) =
6.9676`, the point where the two expressions meet (the paper rounds it to 7),
so the switch loses only the `exp(-2.34 tau_c) = 8e-8` term.

This is a Doppler-core, complete-redistribution, single-flight result, which
is the right class for these lines: their Voigt parameter is `a ~ 1e-8`, so
nothing escapes through a damping wing and the Ly-alpha correction of item
LYA-BETA (Harrington 1973 eq. 40, section 7 above) does not transfer. Nothing
in these winds comes near `a tau > 1e3` in a fine-structure line.

### 10.3 The defect: the argument is short by `sqrt(pi)`

`phi` in the integral (B-7) approximates is a **normalized** profile,
`int phi dx = 1`, so `phi(0) = 1/sqrt(pi)` and the argument `tau` of (B-7) is
the frequency-integrated depth. The line-centre depth is `tau/sqrt(pi)`:

```
tau_dJ = sqrt(pi) * tau_line-centre
```

The independent confirmation is **Hollenbach & McKee (1979), ApJS 41, 555,
eq. (5.10)**, the same physical quantity written in the *line-centre*
convention (they state it explicitly: "we shall quote our final results in
terms of the line-center optical depth `tau = tau'/pi^(1/2)`"):

```
eps(tau) = 1 / (1 + tau [2 pi ln(2.13 + tau^2)]^(1/2))
```

for both faces of the slab. MEASURED: `2 beta_dJ(sqrt(pi) tau)` and
`eps_HM(tau)` agree to 3.4% at worst over `tau = 3` to `1e6` and to five
digits above `tau = 1e3` (`1.07331e-4` against `1.07320e-4` at `tau = 1e3`).
Two independently published fits collapsing onto one function is what fixes
the convention.

`fine_structure_line_transfer` integrates `line_center_opacity_lte`, i.e. it
builds a **line-centre** column and passes it straight in. The `beta` the
code returns is therefore evaluated at an argument `sqrt(pi)` too small and
is too large by (MEASURED, on a `tau` grid):

| line-centre `tau0` | `beta` as called | `beta` at the published argument | ratio |
|---|---|---|---|
| 0.01 | 9.884e-1 | 9.795e-1 | 1.009 |
| 0.1  | 8.916e-1 | 8.186e-1 | 1.089 |
| 1    | 3.862e-1 | 2.373e-1 | 1.627 |
| 7    | 6.095e-2 | 2.889e-2 | 2.110 |
| 1e2  | 2.490e-3 | 1.315e-3 | 1.894 |
| 1e3  | 1.986e-4 | 1.073e-4 | 1.851 |
| 1e5  | 1.512e-6 | 8.314e-7 | 1.818 |
| 1e6  | 1.374e-7 | 7.589e-8 | 1.810 |

The ratio rises to 2.11 at the branch point and tends to `sqrt(pi) = 1.772`.

### 10.4 What it costs in the cases that ship

MEASURED from the stored `backup/regression/*/output/` profiles of
`wasp_full`, `mol_metals`, `mol_ir_bands` and `lower_profile`, with the
line-centre columns built from `Ion_species.txt` by the same expressions
`fine_structure_line_opacity` uses. **The eight lines never get thick**:

| case | thickest line | max `tau0` | at `r/R_p` | `beta` ratio |
|---|---|---|---|---|
| `wasp_full` | `[O I] 63um` | 0.039 | 1.0002 | 1.035 |
| `mol_metals` | `[O I] 63um` | 0.140 | 1.0002 | 1.124 |
| `mol_ir_bands` | `[O I] 63um` | 0.140 | 1.0002 | 1.124 |
| `lower_profile` | `[O I] 63um` | 0.044 | 1.0002 | 1.040 |

Every other line is thinner: `[C I] 370um` peaks at `1.1e-2`, `[C II] 158um`
at `6.0e-3`, `[N II] 205/122um` below `1e-5`, `[O I] 44um` below `4e-8`. The
depths fall by four decades by `r = 1.13 R_p`.

In the two-level form the cooling carries, `Lambda ~ beta A/(beta A + C(1+x))`
with the H-atom de-excitation rate of the `[O I] 63um` channel, the base cell
of the hot Uranus gate is fully trapping-saturated, so the local `[O I] 63um`
rate moves by the full 11%. Integrated over the domain against the `[O I]`
column of `Cooling_breakdown.txt`, the change in the **total** radiative
losses from that channel is

| case | `[O I]` share of the losses | change of the total |
|---|---|---|
| `wasp_full` | 6.2% | 0.003% |
| `mol_metals` | 1.9% | 0.096% |
| `mol_ir_bands` | 1.2% | 0.056% |
| `lower_profile` | 8.4% | 0.003% |

### 10.5 The fix (applied 2026-09-07)

`fine_structure_line_transfer` now forms each cell's depth as `sqrt(pi)`
times the line-centre depth, so every argument handed to
`line_escape_probability_one_face` is the frequency-integrated depth (B-7)
is written in; the function itself stays the literal (B-7). One factor in
one routine. Expected movement: below 0.1% of the radiative losses in every
case that ships, largest in `mol_metals`; the cases are not byte-identical
and the goldens are refreshed at the series gate.

The reason it was made despite the size is that the error is a factor of two, not 12%,
the moment a run has a base thick in `[O I] 63um` -- a denser or more
oxygen-rich lower atmosphere than the gates carry, or the same gate with the
`Base IR field` on, where the same `beta` also scales the incident field
`nbar_fs` and so the heating the layer takes from below.

# Destruction channels of trapped Ly-alpha: measurement and verdict

2026-08-12. Follow-up to `docs/resonance_line_trapping.md` §7, which found that
pure two-level trapping does not explain any suppression of Ly-alpha cooling
(`tau0 = 1.3e8` and `beta ~ 2e-6` at the WASP-121 b base, yet a
cooling-weighted `S = 0.9994`) and left the destruction channels of the trapped
photon as the open question. This memo measures those channels and judges
whether the current treatment is physically consistent. Sections 1-10 are the
measurement, written before any code was touched, and describe the tree as it
then stood; section 11 records the fix that followed on the same day. Every
`file:line` citation in sections 1-10 is a line number of the 2026-08-12 tree
and no longer resolves; the routine and variable names still do.

---

## Verdict

**Three findings, in order of physical seriousness.**

1. **`n2_populations` in `src/modules/radiation/excited_hydrogen.f90` is
   physically wrong.** [2026-08-15: was wrong; fixed 2026-08-12 -- see
   section 11. The rate coefficients now have a single definition in
   `radiation/hydrogen_n2_rates.f90` and the dummy arguments that shadowed the
   statistical weights are gone. The finding is kept in the present tense as
   written, as the record of the defect.] Its dummy arguments `G2s`/`G2p` shadow the module
   parameters `g2s`/`g2p` — Fortran is case-insensitive — so lines 249-251
   evaluate the statistical weights as the Balmer photoionization rate
   `gamma2_bal` instead of 2 and 6. All three collisional rate coefficients
   this memo is about are affected: the `2s -> 1s` and `2p -> 1s` de-excitation
   rates come out 5 to 306 times too small, and the `2p -> 2s` l-mixing rate 3
   times too large. The error is not a constant; it scales with the stellar
   Balmer continuum, so it differs from planet to planet. Measured on the four
   production runs, recomputed from their own dumped state: `n(2s)` is
   over-predicted by 1.5-2.2x, `n(2p)` is unaffected to 2e-5, and the derived
   Balmer proton source and photoelectric heating are over-predicted by 29-41%.
   A re-converged run confirms it (§10). This is wrong independently of the
   effect size and should be fixed.

2. **No destruction channel suppresses Ly-alpha cooling by anything like an
   order of magnitude.** With the corrected weights, a photon trapped at
   Ly-alpha escapes with probability 0.982 (WASP-121 b) to 0.9999
   (HD 189733 b), weighted by where the Ly-alpha cooling is emitted.
   Photoionization of `H(n=2)` is the largest destruction channel, then
   collisional de-excitation, then the `2p -> 2s -> 2gamma` route. The
   distinction that decides the question is that **only collisional
   de-excitation returns energy to the electron gas**; the other two remove the
   Ly-alpha *photon* while the 10.2 eV still leaves the gas, as ionization
   potential that the recombination cascade later radiates or as two-photon
   continuum. They suppress the Ly-alpha line, not the cooling rate. The
   collisionally driven Ly-alpha cooling is suppressed by at most 0.06%.

3. **A real term is missing, but it is a heating source rather than a
   suppression factor.** 58-86% of the `n=2` population is maintained by
   radiative pumping, not by collisional excitation. Collisional de-excitation
   of that pumped population returns energy to the electrons at a rate
   (`Hdx_arr`, gated by `incl_deexc_heat`, default off [2026-08-15: default
   `.true.` since 2026-08-12; section 11 reversed it, for the reason this
   paragraph goes on to give]) equal to 1.8-92% of the
   Ly-alpha cooling on three planets and 5x it on HD 209458 b. The comment in
   `parameters.f90` justifying the default — that it "overlaps the existing HI
   coex cooling" — appears to be mistaken: the Cen (1992) coefficient in
   `Cool_coeff.f90:843` is a one-way, Boltzmann-suppressed excitation rate with
   no density-dependent de-excitation term in it, i.e. the coronal limit in
   which every excitation escapes. Subtracting a de-excitation term is the
   correction to that limit, not a double count.

The historical expectation of a ~10x Ly-alpha cooling suppression is not
supported by any of these channels.

---

## 1. Where each piece lives in the code

| quantity | code location | status |
|---|---|---|
| Ly-alpha cooling `Lambda_coex` | `Cool_coeff.f90:843` `coex_rate_HI`; assembled at `util_ion_eq.f90:640`; dumped as `cool_chan(:,3)` at `:800` | coronal limit; no escape probability, no destruction correction |
| escape probability `beta` | `lya_rt.f90:103-113` (Neufeld/Harrington wing plus Sobolev, combined as `1-(1-b_st)(1-b_sob)`; 2026-08-15: with `Lya absorbing bottom: True` a third factor `(1-b_bot)` joins the product -- see `lya_rt.f90`) | used *only* to build `J_int`; never multiplies the cooling |
| `J_lya` closure | `lya_rt.f90:127-128`, `J_int = Jpref (g1s/g2p) P (1-beta)/(A beta n1s)` | implies `n2p = P/(A beta)`; no destruction in the denominator (§5) |
| (i) `H(n=2)` photoionization | rate `gamma2_bal` at `excited_hydrogen.f90:126`; enters `L2p`/`L2s` at `:258-259`; proton source `:191-195`; heating `:198` | present, and active whenever `use_excited_H` |
| (ii) `2p -> 2s` then two-photon decay | `C2s2p`/`C2p2s` at `:246,251`; `A_2s1s` in `L2s` at `:259` | present, but `C2p2s` is 3x too large (§2) |
| (iii) collisional de-excitation `2 -> 1` | populations: `C2s1s`/`C2p1s` at `:249-250`; heating: `Hdx_arr` at `:200-205` | population rates wrong by 1-3 decades (§2); heating term default off |

So all three channels are present in the `n=2` balance — the Phase-3a
photoionization channel in particular is in `L2p`, `L2s`, the proton source and
the photoelectric heating, and is not double counted. What is missing is not a
channel but the correct rate coefficients for two of them.

Two independent sets of de-excitation coefficients exist inside the same
module. The standalone functions `c2s1s_rate`/`c2p1s_rate` (`:342-357`) take
only `T` as a dummy argument, so their `g1s/g2s` and `g1s/g2p` are the correct
1 and 1/3, and they are what `Hdx_arr` uses. The copies inside
`n2_populations` are the shadowed ones. The heating term and the population
solve therefore disagree with each other about the same rate coefficient.

`n2_populations` has a single caller, so all three `jlya_mode` values behave
identically in this respect.

## 2. The `g2s`/`g2p` shadowing

`excited_hydrogen.f90:40` declares

```
real*8, parameter :: g1s = 2.0d0, g2s = 2.0d0, g2p = 6.0d0
```

and `:223` opens

```
subroutine n2_populations(T, n1s, ne_l, Jlya, G2s, G2p, n2s, n2p)
```

Fortran does not distinguish `G2s` from `g2s`, so inside the subroutine both
spellings refer to the dummy argument, which the caller fills with
`gamma2_bal` twice (`:178-179`). Lines 249-251 then evaluate

| line | intended factor | actually evaluated | error |
|---|---|---|---|
| `C2s1s = 1.21e-8 (1/t4)^0.455 * (g1s/g2s)` | 1 | `2/gamma2_bal` | too small by `gamma2_bal/2` |
| `C2p1s = 1.71e-8 (1/t4)^0.077 * (g1s/g2p)` | 1/3 | `2/gamma2_bal` | too small by `gamma2_bal/6` |
| `C2p2s = C2s2p * (g2s/g2p)` | 1/3 | 1 | 3x too large |

`B12_lya = (g2p/g1s)*B21_lya` at `:43` is at module scope and is unaffected.
`gamma2_bal` still enters `L2p` and `L2s` correctly, because there the dummy
argument is what is wanted.

**Verification.** An instrumented scratch build was made to print the solver's
own coefficients. `M12` and `M21` came out numerically equal
(`8677.9501244847615` both) where they should differ by exactly 3, and a direct
print gave `g2s = g2p = 612.900334...`, which is the WASP-121 b `gamma2_bal`.
Reconstructing the solve in Python with the shadowed factors reproduces the
`Excited_H.txt` populations of the production runs to 5e-6 - 1.5e-5; with the
intended factors it does not (40-60% apart).

Per planet, with `Gamma_2` from each run's `Excited_H.txt` header:

| planet | `Gamma_2` [s^-1] | `C2s1s` too small by | `C2p1s` too small by | `C2p2s` too large by |
|---|---|---|---|---|
| HD 189733 b | 10.7 | 5.3x | 1.8x | 3x |
| WASP-52 b | 21.0 | 10.5x | 3.5x | 3x |
| HD 209458 b | 68.7 | 34x | 11.4x | 3x |
| WASP-121 b | 612.9 | 306x | 102x | 3x |

Effect on the converged production outputs, recomputed cell by cell from the
state each run dumped, volume-integrated over `r >= 1`:

| planet | `n(2s)` corrected/current | `n(2p)` corrected/current | `n(2)` total, hence `S_proton` and `H_pe` |
|---|---|---|---|
| HD 189733 b | 0.481 | 1.0000 | 0.709 |
| HD 209458 b | 0.460 | 1.0000 | 0.713 |
| WASP-52 b | 0.664 | 1.0000 | 0.778 |
| WASP-121 b | 0.526 | 1.0000 | 0.727 |

`n(2p)` is insensitive because `L2p` is dominated by `A_2p1s = 6.3e8`, far
above any of the corrupted collisional terms. `n(2s)` is not, because `M21`
sets the `2p -> 2s` feed and l-mixing supplies 94-99% of `L2s`. The clearest
signature is at the WASP-121 b base, where the current output gives
`n2s/n2p = 1.007` against `0.363` with the corrected weights — the latter being
the statistical ratio `g2s/g2p = 1/3` that strong l-mixing should enforce.

The Python transmission tool carries the same algorithm
(`exhale_transit_lib.py:138-173`) with the same variable names, but Python is
case-sensitive, so its `g2s`/`g2p` are the statistical weights and its
populations are correct. Only the Fortran transliteration is affected.

## 3. Branching ratios of a trapped Ly-alpha photon

For an atom in `2p` the competing rates are

```
escape         beta_tot * A_2p1s
coll 2p->1s    ne * q_2p1s
photoion n=2   Gamma_2 (Balmer continuum)
2p->2s->2gam   ne * C2p2s * A_2s1s/(A_2s1s + ne q_2s1s + ne C2s2p + Gamma_2)
```

with `beta_tot` rebuilt exactly as `lya_rt.f90` does (wing escape from the
dumped `tau_Lya`, combined with the Sobolev channel from the velocity
gradient). Weighting by where the Ly-alpha cooling is emitted:

| planet | escape | coll `2p->1s` | photoion `n=2` | `2p->2s->2gamma` |
|---|---|---|---|---|
| HD 189733 b | 0.999854 | 3.0e-5 | 9.2e-5 | 2.4e-5 |
| HD 209458 b | 0.999483 | 4.4e-6 | 5.0e-4 | 1.8e-5 |
| WASP-52 b | 0.999372 | 5.5e-5 | 5.1e-4 | 6.5e-5 |
| WASP-121 b | 0.982140 | 5.7e-4 | 1.7e-2 | 7.5e-5 |

Destruction is confined to the dense base, where `beta` is smallest. At
`r = 1.000`:

| planet | `beta` | escape | photoion `n=2` | coll `2p->1s` |
|---|---|---|---|---|
| WASP-121 b | 2.0e-6 | 0.671 | 0.321 | 6.2e-3 |
| HD 209458 b | 5.4e-7 | 0.824 | 0.168 | 2.0e-3 |
| WASP-52 b | 2.0e-6 | 0.979 | 1.6e-2 | 2.0e-3 |
| HD 189733 b | 1.7e-6 | 0.982 | 9.8e-3 | 5.3e-3 |

Photoionization dominates the destruction everywhere, by one to two decades
over collisional de-excitation, and it tracks `Gamma_2`: WASP-121 b, orbiting
the hottest star of the four, is the only case where destruction reaches the
percent level in the integrated budget.

**Why this does not become cooling suppression.** Of the three channels, only
collisional de-excitation returns the 10.2 eV to the electron gas. A
photoionization from `n=2` converts the excitation energy into ionization
potential, which the recombination cascade later radiates; the two-photon route
radiates it as continuum. Both leave the electron gas exactly as an escaping
Ly-alpha photon would. The branching table above therefore bounds the
suppression of the Ly-alpha *line*; the suppression of the *cooling* is bounded
by the collisional column alone. Measured as `1 - f_coll`, the collisionally
driven Ly-alpha cooling retains

| planet | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| `1 - f_coll` | 0.999970 | 0.999996 | 0.999945 | 0.999430 |

i.e. a suppression of 0.06% at worst. That is consistent with, and about a
decade below, the `1 - S = 6e-4` two-level estimate of the earlier memo.

## 4. The energy ledger, and the term that is switched off

`coex_HI` is a one-way rate. The net loss of the electron gas through the
`n=2` manifold is

```
Lambda_net = E21 ne [ n1s (q_1s2s + q_1s2p) - n2s q_2s1s - n2p q_2p1s ] - n2 hpe2
```

The first bracket term is `coex_HI`, the second is `Hdx_arr`, the last is
`Hpe_arr`. The code has all three; only `coex_HI` and `Hpe_arr` are active by
default. Recomputed with corrected populations, integrated over `r >= 1`:

| planet | `coex_HI` / total cooling | `Hdx`/`coex` | `Hpe`/`coex` | `1 - (Hdx+Hpe)/coex` |
|---|---|---|---|---|
| HD 189733 b | 0.150 | 0.031 | 0.0023 | 0.966 |
| WASP-121 b | 0.150 | 0.018 | 0.022 | 0.960 |
| WASP-52 b | 9.8e-3 | 0.915 | 0.289 | -0.20 |
| HD 209458 b | 4.9e-4 | 5.16 | 12.1 | -16.2 |

For HD 189733 b and WASP-121 b, where Ly-alpha carries 15% of the radiative
losses, the whole correction amounts to 4%. For the other two, Ly-alpha cooling
is 0.05-1% of the losses and the pumped `n=2` manifold is a net heating term;
that is not a suppression of anything, and `Hpe`, which is already active, is
the larger half of it. The omitted `Hdx` alone is 0.25-0.90% of the total
radiative losses of each planet, the same order as the `Hpe` term the code
already carries.

The reason `Hdx/coex` is orders of magnitude above the `f_coll` of §3 is that
`Hdx` de-excites the whole `n=2` population, most of which the Ly-alpha field
pumped rather than a collision:

| planet | HD 189733 b | HD 209458 b | WASP-52 b | WASP-121 b |
|---|---|---|---|---|
| pumped fraction of `n(2)` | 0.82 | 0.86 | 0.58 | 0.77 |

So `Hdx` is best read as an independent heating channel — absorbed Ly-alpha
thermalized by a collision — rather than as a suppression factor on the
Ly-alpha cooling.

Two smaller consistency notes on this ledger:

* The two sides use different atomic data. `E21*(C1s2s + C1s2p)` from the
  Christie coefficients is 0.78-0.91 of the Cen (1992) `coex_HI` over
  5e3-5e4 K. The residual is plausibly the `n >= 3` excitation that the Cen fit
  lumps in, so the two look consistent at the 10-30% level, but the loss and the
  return term are not derived from one set of rates.
* Nothing in the current treatment is double counted. A collisional excitation
  costs the electrons 10.2 eV once (`coex_HI`); re-absorption of the trapped
  photon costs nothing; a collisional de-excitation returns it once (`Hdx`).
  The loop closes.

Note that the regression cases `wasp_full` and `wasp_he23off` do set
`Deexc heat: True`, so `Hdx` is exercised there even though the four production
planet runs leave it off.
[2026-08-15: `Deexc heat` is on by default since 2026-08-12, so the explicit key
in those cases now restates the default and the production runs no longer leave
the term off.]

## 5. The `J_int` closure omits the destruction channels

`lya_rt.f90:127` builds the internal field from `n2p = P/(A beta)`, with `P`
the recombination plus collisional Ly-alpha production. The trapped-photon
balance should carry the destruction rates in the denominator,
`n2p = P/(A beta + D)`. The resulting over-estimate of `J_int` is exactly
`1/f_escape` from §3:

| planet | Ly-alpha-weighted `f_esc` | `J_int` over-estimate, weighted | worst local |
|---|---|---|---|
| HD 189733 b | 0.9998 | 1.0002 | 1.018 |
| WASP-52 b | 0.9993 | 1.0007 | 1.021 |
| HD 209458 b | 0.9994 | 1.0006 | 1.214 |
| WASP-121 b | 0.9788 | 1.022 | 1.489 |

So the closure is self-consistent to better than 3% integrated, and up to 1.5x
optimistic in the first cells of WASP-121 b. That the closure is otherwise
consistent with what `n2_populations` then solves was checked directly: where
the internal source dominates, at the WASP-121 b base, the `P/(A beta)` value
and the solved `n2p` agree to 1%. Elsewhere the solved `n2p` is larger because
the stellar beam `J_star` dominates the pumping, which is expected and not an
inconsistency.

Separately, `J_star` carries the tuned parameter `lya_star_boost` (default 5,
described in `parameters.f90` as tuned to Huang et al. 2023 Fig. 11). It sets
`n2p` wherever the stellar beam penetrates, and hence the Balmer proton source,
the de-excitation heating and the H-alpha opacity. It is documented as a
tuning; it is worth noting that finding 1 changes the populations it was tuned
against.

## 6. A smaller point: the recombination source uses `ne^2`

`excited_hydrogen.f90:260-261` writes the cascade source as `a2p*ne_l**2` and
`a2s*ne_l**2`. The physical rate is `alpha_2 ne n_HII`. The two agree only in a
pure-hydrogen, fully ionized gas. At the base the electrons come largely from
helium and metals while hydrogen is still neutral, so `ne/n_HII` reaches 2.0
(WASP-121 b), 6.2e2 (WASP-52 b), 3.7e3 (HD 189733 b) and 5.1e6 (HD 209458 b).

Measured effect of using `ne n_HII` instead: `n(2p)` unchanged to four digits,
because the pump dominates, and `Hdx` lower by 3.4-5.4% integrated. Small, but
wrong as written.

## 7. LaRT import mode

`jlya_mode = 1` reads an external `J_lya(r)` and feeds it to the same
`n2_populations`, so it inherits finding 1 unchanged; the escape-probability
closure of §5 is bypassed, which is the point of the mode. Consuming an
external mean intensity through a statistical equilibrium that uses the full
`A_2p1s` plus explicit `B12 J` and `B21 J` terms is the right way to do it, so
that part is consistent. One cosmetic consequence: `taulya` is forced to zero
in this mode (`:148`), so the `tau_Lya` column of `Excited_H.txt` is blank for
LaRT runs. Nothing reads `taulya` outside `jlya_mode = 0` and the output file.

## 8. Suggested fixes, in order

All four were implemented on 2026-08-12; section 11 records what was changed and
what it measured.

1. Rename the `n2_populations` dummy arguments (`G2s`/`G2p` -> `gam_2s`/
   `gam_2p`) so the statistical weights resolve to the module parameters. This
   changes results: expect the Balmer proton source and photoelectric heating
   to drop by roughly 30% and `n(2s)` by a factor ~2. `wasp_full` and
   `wasp_he23off` goldens will need refreshing; `mol_base_handoff` sets no
   stellar `T_eff` and should be unaffected. The H-alpha and H-beta comparisons
   in the paper runs are downstream of `n(2s)`.
2. Decide `incl_deexc_heat` on physics rather than on the overlap argument. If
   it is turned on, correct the justification comment in `parameters.f90` at the
   same time.
3. Put the destruction rates into the `J_int` denominator of `lya_rt.f90:127`,
   `n2p = P/(A beta + ne q_2p1s + Gamma_2 + ne C2p2s P_2gamma)`. Worth about 2%
   integrated and 1.5x in the innermost WASP-121 b cells.
4. Replace `ne_l**2` by `ne_l*n_HII` in the cascade source, which means passing
   `n_HII` into `n2_populations`.

Expected magnitudes if all four are done: Ly-alpha cooling changed by a few
percent, not 10x; Balmer proton source and `n=2` photoelectric heating down
~30%; the H-alpha transit depth moved by an amount set by `n(2s)` and `n(2p)`
together, which is the quantity to check first because `n(2s)` moves by ~2x
while `n(2p)` barely moves. Mass-loss rates are expected to move little (§10).

## 9. Reproduction

The measurement was run from four throwaway scripts in the session scratch
directory (`.../scratchpad/lya_destroy/`): `destroy.py` (branching ratios and
the energy ledger, read-only from each planet's `output/`), `verify.py`
(transcription check and the `ne^2` test), `impact.py` (shadowed vs. corrected
populations), `cmp.py` (paired scratch runs). **None of them was kept**: a
session scratch directory does not survive the session, and they were never
copied into the repository, so the tables above cannot be regenerated by
re-running them. What they did is described in enough detail here to redo, and
every input they read is a committed `output/` directory. The shadowing itself
was confirmed by an instrumented build that prints `g2s`, `g2p`, `C2s2p`,
`C2p2s` and the two off-diagonal terms from inside `n2_populations`.

One trap worth recording: `rsync -a` preserves source mtimes, so restoring an
unpatched file over a patched one leaves it *older* than the object built from
the patch and `make` silently skips the rebuild. Two binaries built that way
came out byte-identical. Touch the file after restoring it.

## 10. Paired runs with the shadowing fixed

Built from the current tree in a scratch copy, unpatched against
`G2s`/`G2p` renamed, same input, same step budget (4000) with the Newton finish,
`OMP_NUM_THREADS=8`.

| planet | log10 Mdot current | log10 Mdot fixed | max local `dT` | `n(2s)` f/c | `n(2)` f/c | `coex_HI` f/c |
|---|---|---|---|---|---|---|
| WASP-121 b | 13.340 | 13.340 | 0.30% | 0.555 | 0.732 | 1.037 |
| WASP-52 b | 9.610 | 9.610 | 0.64% | 0.496 | 0.714 | 1.029 |
| HD 189733 b | 9.130 | 9.140 | 5.5% | 0.491 | 0.714 | 1.026 |
| HD 209458 b | 9.380 | 9.400 | 7.3% | 0.467 | 0.718 | 1.138 |

The `n=2` response is remarkably uniform and matches the post-hoc estimates of
§2 to within the feedback the post-hoc calculation does not carry: the total
`n=2` population falls to 0.714-0.732 on every planet and `n(2s)` to 0.47-0.56,
so the Balmer proton source and the photoelectric heating drop by a little
under 30% everywhere. The Ly-alpha cooling rises 2.6-13.8%.

The wind response is not uniform. WASP-121 b and WASP-52 b are unchanged in
`Mdot` to three decimals with temperatures moving less than 0.7%. HD 189733 b
and HD 209458 b both shift: `Mdot` up by 0.01 and 0.02 dex (2-5%) with local
temperature changes of 5.5% and 7.3%. A plausible reading is that the two that
move are the two where the `n=2` heating terms are largest relative to the
local energy balance — HD 209458 b in particular carries the largest
`Hpe`/`coex` ratio of the four (§4) — but this memo did not decompose the
energy balance of the re-converged runs to confirm that, so it is an
interpretation rather than a measurement.

The practical consequence is that the fix should be treated as capable of
moving the wind, not only the diagnostics. The paper runs will need
re-convergence, not merely a recomputed transit.

**Not checked:** the transmission spectra, which were not recomputed; which
term drives the temperature shift on HD 189733 b and HD 209458 b; the
neutral-hydrogen contribution to `q_2p1s` at the base; any planet outside the
four production runs.

## 11. Fix

Implemented 2026-08-12. All four items of section 8 were applied, plus the
`taulya` gap of section 7 and one structural change that removes the condition
that produced finding 1.

### What changed

**New module `src/modules/radiation/hydrogen_n2_rates.f90`.** The n=2 / Ly-alpha
atomic data and every collisional rate coefficient now have a single definition.
They previously existed in three copies -- written out inside `n2_populations`,
again as the standalone `c2s1s_rate`/`c2p1s_rate` a hundred lines below, and a
third time inside `lya_rt.f90`. Finding 1 is exactly what that duplication
costs: the two copies in the same file disagreed by up to two decades and
nothing in the code could notice. `n2_populations` and `jlya_escape_prob` now
call `c1s2s_rate`, `c1s2p_rate`, `c2s2p_rate`, `c2s1s_rate`, `c2p1s_rate`,
`c2p2s_rate`, `alpha_B_hydrogen`, `alpha_2s_hydrogen`, `alpha_2p_hydrogen`, so
the statistical weights are read at module scope and a dummy argument cannot
reach them.

1. **Shadowing.** The `n2_populations` dummy arguments are `gam_ion_2s` /
   `gam_ion_2p`, and the reverse-rate expressions they used to shadow are gone
   from the routine body (previous paragraph). A name-by-name audit of the two
   modules against the `global_parameters` scope found two more shadows of the
   same kind, both currently harmless and both renamed: `gamma_n2_balmer` and
   `heat_n2_balmer` declared `integer, parameter :: ng = 400` for their
   frequency-sample count, which shadows the global ghost-cell count `Ng`. They
   are now `n_nu`. After the renames the audit is clean for
   `excited_hydrogen.f90`, `lya_rt.f90` and the new module, checked against the
   union of `global_parameters`, the shared module and each file's own module
   scope.

   Outside these files the same audit flags `rr_badnell(T,A,B,T0,T1,Cc,T2)` in
   `Cool_coeff.f90`, whose `T0` shadows the global temperature normalization.
   That one is deliberate -- `T0` is the published Badnell fit parameter and the
   routine uses the dummy, which is what is wanted -- so it was left alone and
   is recorded here as a known trap rather than changed.

2. **Recombination source.** `n2_populations` takes `nHII_l` and the cascade
   source is `alpha_2l*ne*nHII`. `excited_H_update` passes `max(nhii(j),0)`.
   The same correction was made in the Python transmission tool
   (`exhale_transit_lib.py::n2_populations`, called from `EXHALE_transit.py`)
   and in `examples/tpm_halpha_lart2d.py`, which carries its own transcription:
   both had `ne**2` as well. Their signatures gained an `nHII` argument.

3. **`J_int` closure.** `hydrogen_n2_rates::n2p_destruction_rate(T, ne, Gamma_2s,
   Gamma_2p)` returns `ne q_2p1s + Gamma_2p + ne C_2p2s P_2gamma` with
   `P_2gamma = A_2s1s/(A_2s1s + ne q_2s1s + ne C_2s2p + Gamma_2s)`, and
   `lya_rt.f90` divides by `A_2p1s beta + D` instead of `A_2p1s beta`.

   `tau_Lya` is no longer forced to zero in the LaRT import mode. The top-down
   line-center optical depth was computed by two identical inline loops (one in
   `jlya_escape_prob`, one in the `jlya_mode = 0` branch); it is now the single
   routine `lya_rt::lya_line_center_optical_depth`, which all three modes call.
   The relocation is exact arithmetic -- same expression, same operands, same
   order -- so it changes no number in modes 0 and 2. Verified in mode 1 on a
   `wasp_full` variant fed a synthetic `jlya_rt.txt`: `tau_Lya` runs from
   `2.4e7` at the base to 0 at the outer boundary, where the column used to be
   blank.

4. **`Deexc heat` default.** Now on. The justification comment in
   `parameters.f90` ("overlaps the existing HI coex cooling") is replaced by the
   argument of section 4: the Cen (1992) coefficient is the one-way coronal
   excitation rate with no de-excitation term in it, so subtracting `Hdx` is the
   correction to that limit rather than a double count, and most of the n=2
   population it de-excites was pumped by Ly-alpha rather than by a collision.
   The key now also accepts `False`; it previously ignored it, since the
   defaults block reset the flag and only `True` was tested.

   Caveat worth keeping in view: with the term on, the energy budget depends
   more strongly on `lya_star_boost` (default 5, tuned to Huang et al. 2023
   Fig. 11), because `Hdx` scales with the pumped n=2 population. The larger
   `Hpe` term already carried that dependence, so this is a change of degree.

### Audit of the Ly-alpha cooling suppression factor (section 3's question, in code)

**The energy equation applies no escape-probability factor to the Ly-alpha
cooling, and never did.** `eval_cool` calls `coex_rate_HI`
(`Cool_coeff.f90`), the bare Cen (1992) coefficient
`7.5e-19/(1+sqrt(T/1e5)) exp(-118348/T)`, forms
`coex = coeff_coex_rate_HI*nhi + ...` and
`cool = ne*(brem + coex + reco + coio) + cool_M + cool_H3p`
(`util_ion_eq.f90`, `eval_cool`). Neither `beta` nor
`n2p` nor `Jlya` appears anywhere in that path; `beta_tot` exists only inside
`lya_rt.f90` and is read only at the two lines that build `J_int` and `J_star`.
The path does not branch on `jlya_mode`, so the LaRT import mode is identical.
[2026-08-15: still true of the cooling path. Inside `lya_rt.f90`, `beta_tot`
gains the downward-escape factor when `Lya absorbing bottom: True`, which
changes `J_int`/`J_star` but not this cooling term.]

That is the physically right form given the measurement of section 3: the
cooling is set by the excitation rate `P`, and the collisional destruction that
could suppress it is at most 0.06%. The historical expectation of a ~10x
Ly-alpha cooling suppression from radiative transfer never entered the code, and
should not: it would have been the error of multiplying a coronal excitation
rate by an escape probability. No change was made here.

One labeling nuance, unchanged: `cool_chan(:,3)` is called `coex_HI [Lya]`, but
the Cen fit lumps the `n >= 3` H I excitation in with Ly-alpha (section 4 puts
`E21 (C_1s2s + C_1s2p)` at 0.78-0.91 of it), while `Hdx` de-excites only n=2.
The two sides of the ledger therefore come from different atomic data at the
10-30% level.

### Measured effect

Regression cases, single-threaded, built cumulatively so each row isolates one
item. `mol_base_handoff` sets no stellar `T_eff`, so `use_excited_H` is false and
every item is inert there -- byte-identical at every stage, as predicted.
`wasp_full` and `wasp_he23off` pin `Deexc heat: True` in their `input.inp`, so
item 4 moves nothing in them.

| build | `wasp_full` log10 Mdot | `wasp_he23off` | `mol_base_handoff` |
|---|---|---|---|
| before | 13.23039 | 13.23106 | 10.58396 |
| + item 1 (shadowing) | 13.22963 | 13.23024 | byte-identical |
| + item 2 (`ne nHII`) | 13.22959 | 13.23018 | byte-identical |
| + item 3 (`J_int`) | 13.22966 | 13.23026 | byte-identical |
| + item 4 (default) | 13.22966 | 13.23026 | byte-identical |

Largest local temperature change introduced by each item on `wasp_full`: 0.57%
at `r = 1.009` (item 1), 0.081% at `r = 1.032` (item 2), 0.147% at `r = 1.042`
(item 3).

Volume-integrated `r >= 1` diagnostics of `wasp_full`, cumulative, relative to
the pre-fix run:

| build | `n(2s)` | `n(2p)` | `S_proton` = `H_pe` | `Hdx` | `J_int` | `coex_HI` |
|---|---|---|---|---|---|---|
| + item 1 | 0.505 | 0.933 | 0.687 | 0.596 | 0.988 | 1.030 |
| + item 2 | 0.474 | 0.930 | 0.668 | 0.575 | 0.987 | 1.032 |
| + item 3 | 0.460 | 0.872 | 0.635 | 0.558 | 0.981 | 1.033 |

Two production planets re-converged from their stored initial conditions in a
scratch copy (8 threads; the planet directories were not touched). `fix123` is
items 1-3 with `Deexc heat: False` appended, `fix` is the shipped default.
(**Dated 2026-08-13**: the planet folders themselves were subsequently
re-converged with all four items in, giving `log10 Mdot = 9.46` (HD 209458 b),
`9.14` (HD 189733 b), `11.70` (WASP-52 b) and `13.20` (WASP-121 b) in each
folder's `run_20260812_lyafix.log`. Those are the numbers `paper/ms.tex` and
`python/paper_data.py` carry; the scratch values below are the isolation A/B.)

| planet | log10 Mdot before | items 1-3 | + item 4 |
|---|---|---|---|
| WASP-121 b | 13.2075 | 13.2035 | 13.2030 |
| HD 189733 b | 9.1286 | 9.1372 | 9.1392 |
| HD 209458 b | 9.3790 | 9.4540 | 9.4581 |

HD 209458 b moves most, by 0.075 dex (19%) on items 1-3 alone, with the
temperature moving up to 19% at `r = 1.028`. That is the planet section 4
identifies as the extreme case: its Ly-alpha cooling is 5e-4 of the radiative
losses, so the `n=2` heating terms, not the cooling, dominate that channel, and
a 30-40% correction to them is felt directly. It is also the planet
`docs/lya_deexcitation_heating.md` recorded as going NaN at the base with
`Deexc heat: True`. That no longer reproduces: with the corrected populations it
re-converges from its stored initial condition in ~2000 steps and a cold start
survives a 20000-step cap, both without a NaN. Item 4 there moves the first cell
by 17% in temperature and the domain outside `r = 1.01` by up to 9%.

Temperature: items 1-3 move WASP-121 b by at most 3.0% (at `r = 1.075`) and
HD 189733 b by 10.4% at `r = 1.003`, 2.4% once `r > 1.01`. Item 4 moves
WASP-121 b by 0.89% (first cell) and 0.58% outside `r = 1.01`; on HD 189733 b it
raises the first cell from 551 to 745 K (35%) and moves the rest of the domain
by at most 0.89%. The de-excitation heating is largest exactly where the gas is
coldest and densest, so the base cell is where it shows.

Transit depths, `1 - min(transmission)`, before -> after all four items, from
`EXHALE_transit.py` run on the same scratch copies:

| planet | Ly-alpha | H-alpha | H-beta | He 10830 |
|---|---|---|---|---|
| WASP-121 b | 3.7531e-2 -> 3.7531e-2 | 3.3674e-2 -> 3.3543e-2 | 2.3504e-2 -> 2.3278e-2 | 3.7104e-2 -> 3.7083e-2 |
| HD 189733 b | 3.7667e-1 -> 3.7889e-1 | 5.9418e-3 -> 6.0221e-3 | 3.5779e-3 -> 3.6395e-3 | 2.8163e-2 -> 2.8599e-2 |

The Balmer lines move by 0.4-1.7% even though `n(2s)` halves, because their
opacity is carried by `n(2p)`, which barely moves -- the same insensitivity
section 2 predicted. The WASP-121 b Ly-alpha core is saturated out to the
geometric limit, so its depth is unchanged to every digit printed.

Runtime checks: a scratch build at `-O1 -fopenmp -g -fcheck=bounds,do,mem`
running the `wasp_full` configuration (excited H, `jlya_mode = 2`, metals) for
1200 steps is clean.

**Not checked:** planets other than the two above; `jlya_mode = 1` beyond the
`tau_Lya` column and a 60-step start; the regression goldens were left as they
were -- they were already stale against the pre-fix tree for the two `wasp`
cases, and `mol_base_handoff` still reproduces its golden bit for bit.

# Lyman-alpha de-excitation heating in EXHALE

This note documents the collisional de-excitation heating term controlled by the
`Deexc heat` key in `input.inp`: what it computes, the formula as coded, where it
lives, the sources it is built from, and its numerical behavior (it destabilizes
the low-irradiation base of HD 209458 b).

## Physical picture

Lyman-alpha photons trapped in the escaping atmosphere pump hydrogen from the
ground state into `n = 2` (mostly the 2p sublevel, with 2s fed by cascades and
mixing). Some of those excited atoms are collisionally de-excited back to the
ground state by electron impact before they radiate, converting the 10.2 eV
internal energy of the transition into thermal energy of the gas. This is the
gas-heating counterpart of Lyman-alpha "trapping": instead of the line photon
escaping and cooling the gas, the excitation is re-thermalized locally.

The `n = 2` populations (`n_2s`, `n_2p`) come from the two-level (2s/2p) steady
state of Christie, Arras & Li (2013), driven by the local Lyman-alpha mean
intensity `J_Lya`; the de-excitation heating is the super-elastic (detailed
balance) reverse of the 1s to 2s/2p electron-impact excitation used in that
population solve.

### Relation to Lyman-alpha cooling

This heating term is distinct from, and complementary to, the standard
Lyman-alpha *cooling* that EXHALE inherits from ATES (the `collisional
excitation` term of the cooling module, `Cool_coeff.f90`; the dominant coolant
of an atomic-hydrogen escaping wind). That cooling assumes the optically thin
limit: every collisional excitation of neutral H is followed by a radiative
decay whose Lyman-alpha photon escapes and removes the energy. The de-excitation
heating is the optically thick correction to that picture: where Lyman-alpha is
trapped, a fraction of the excited atoms are collisionally de-excited instead of
radiating, so their energy stays in the gas. The two act on the same `n = 2`
atoms in opposite directions — pure cooling in the thin limit, reduced net
cooling (or net heating) once trapping is included.

## Formula as coded

The heating rate per unit volume is

```
H_dx = n_e * dE_21 * ( k_2s1s(T) * n_2s + k_2p1s(T) * n_2p )      [erg cm^-3 s^-1]
```

with `dE_21 = 1.634e-11 erg` (the 10.2 eV 1s-2s/2p energy gap) and the
electron-collision de-excitation rate coefficients (t4 = T / 10^4 K)

```
k_2s1s(T) = 1.21e-8 * (1/t4)^0.455 * (g_1s / g_2s)      [cm^3 s^-1]
k_2p1s(T) = 1.71e-8 * (1/t4)^0.077 * (g_1s / g_2p)      [cm^3 s^-1]
```

The rate coefficients are the Christie et al. (2013) Table 2 forward
excitation rates (originally from Janev et al. 2003) put in detailed-balance
de-excitation form; the statistical weights `g` supply the reverse-rate ratio.
The heating expression itself is the standard collisional-thermalization
assembly of these rates and the transition energy — it is constructed in the
code rather than transcribed from a single numbered equation in any one paper.

## Implementation

- Heating term: [excited_hydrogen.f90:199-206](EXHALE/src/modules/radiation/excited_hydrogen.f90#L199-L206) (`excited_H_update`, fills `Hdx_arr`).
- Rate coefficients: [excited_hydrogen.f90:342-357](EXHALE/src/modules/radiation/excited_hydrogen.f90#L342-L357) (`c2s1s_rate`, `c2p1s_rate`); the energy gap `E21_erg` at line 50.
- `J_Lya` that populates `n = 2`: `J_Lya ~ 0.1 F_LyC / dnu_D` (Huang et al. 2017, Eq. 6), [excited_hydrogen.f90:157](EXHALE/src/modules/radiation/excited_hydrogen.f90#L157).
- Runtime toggle: `Deexc heat: True` sets `incl_deexc_heat` in [input_read.f90:344-346](EXHALE/src/modules/files_IO/input_read.f90#L344-L346); the flag is declared in `parameters.f90`.
- Injection into the energy equation: [ionization_equilibrium.f90:220-224](EXHALE/src/modules/radiation/ionization_equilibrium.f90#L220-L224) (`heat = heat + heat_balmer`, with `heat_balmer = Hpe_arr + Hdx_arr`).

## Source references

- **Christie, Arras & Li (2013), ApJ 772, 144** (`2013ApJ...772..144C`) — the
  `n = 2` (2s/2p) population model and the collisional de-excitation rate
  coefficients used here (Table 2, detailed-balance form; forward rates from
  Janev et al. 2003). Cited in the module header of `excited_hydrogen.f90`.
- **Huang, Arras, Christie & Li (2017), ApJ 851, 150** (`2017ApJ...851..150H`)
  — the Lyman-alpha mean intensity `J_Lya` that sets the `n = 2` population.
- **Huang, Koskinen, Lavvas & Fossati (2023), ApJ 951, 123**
  (`2023ApJ...951..123H`) — the framing that this collisional de-excitation
  channel realizes the re-thermalization of trapped Lyman-alpha, rather than
  applying an escape-probability factor to Lyman-alpha cooling.

The atomic data (rates and the 10.2 eV gap) are attributed explicitly in code
comments to Christie et al. (2013); the physical justification for treating the
channel as local heating follows Huang et al. (2023). Further documentation:
`docs/Update_EXHALE.tex` (physics/atomic-data section and the Huang 2023
trapping rationale) and `docs/transmission_spectrum.tex` (the n=2 hydrogen
population section and reference list).

## Domain of validity and numerical behavior

The term is physically appropriate where Lyman-alpha pumping populates `n = 2`
and the electron density is high enough for collisional de-excitation to compete
with radiative decay — i.e. the denser, partially ionized layers of a strongly
irradiated escaping atmosphere.

**Numerical caveat (HD 209458 b).** With `Deexc heat: True`, HD 209458 b (the
lowest-irradiation case of the paper set, `log L_X = 27.2`, `log L_EUV = 27.9`)
develops a base instability: the run reaches NaN at the innermost cell
(`r ~ 1.0002 R_p`) after a few thousand steps. Isolation runs show the
instability is specific to this heating term — with the in-situ Lyman-alpha
field on but `Deexc heat` off, the same setup is stable, and the more strongly
irradiated planets tolerate the term. It is not yet established whether this is a
physical heating runaway (Lyman-alpha heating outrunning the local cooling at the
cool HD 209458 b base) or a numerical stiffness of the term as currently coupled
to the energy update; it interacts with the separately documented base-breathing
behavior of HD 209458 b / HD 189733 b. Until this is resolved, the term is kept
off for the low-irradiation cases, with the in-situ Lyman-alpha field itself left
on, and the choice recorded per run.

# Transmission spectra in `TPM.py`: H-alpha/H-beta and the metal resonance doublets

Quick reference for the H-alpha (6562.8 Å, n=2→n=3) and H-beta
transmission spectra (non-LTE n=2 population, Christie+2013 Ly-alpha
pumping) and the metal resonance doublets (Mg II h&k, Ca II H&K, Na I D)
in `TPM.py`. Full physics and equations are in
`transmission_spectrum.tex` / `.pdf`.

## What it does

`TPM.py` already builds transit transmission spectra (impact-parameter
LOS Voigt integration + disk average + instrument/rotation convolution)
for **He I 10830 Å** and **Ly-alpha 1215.67 Å**. The new branch adds
**H-alpha**, which absorbs out of the `n=2` hydrogen level. The `n=2`
population is computed from the non-LTE balance of
**Christie, Arras & Li (2013, ApJ 772, 144)** — collisions,
recombination cascade, two-photon decay, 2s↔2p ℓ-mixing, and
**Ly-alpha radiative pumping** (1s↔2p) — solved as a 2×2 system for the
2s and 2p densities per radial cell.

## How to use

H-alpha is **always computed**; only the source of the Ly-alpha mean
intensity `J_lya(r)` (which ATES does **not** produce) depends on whether
you supply a file:

**Option 1 — supply `J_lya(r)` as a file.** Set `Jlya_file` in `TPM.py`:

```python
Jlya_file = path + '/Jlya.txt'
```

with a two-column text file:

| col 1 | col 2 |
|-------|-------|
| `r / R_p` (same radial coordinate as the ATES output) | `J_lya` = Ly-alpha mean intensity `J_nu` at line center, **cgs** `erg s^-1 cm^-2 Hz^-1 sr^-1` |

Template: `inputdata/Jlya.txt.example`.

**Option 2 — default simple estimate (no file).** With `Jlya_file = ''`
(the default), `J_lya` is estimated following **Huang et al. (2017,
ApJ 851, 150), Eq.(6) and related text**:

```
J_lya(r) ~ 0.1 * F_LyC / Dnu_D(r)
Dnu_D    = nu_Lya * sqrt(2 kB T / m_p) / c     # local Ly-alpha Doppler width
```

`F_LyC` is the **deposited (absorbed) Lyman-continuum flux** (each LyC
ionization balanced by a recombination -> one Ly-alpha photon). It is
**computed from the actual stellar input**, not assumed:

```
F_LyC = xi * [ 10^LEUV / (4*pi*a^2) ] * ( 1 - exp(-sigma_LyC * N_HI) )
        ^xi   ^ incident stellar LyC     ^ fraction absorbed
```

where `LEUV` (EUV-band luminosity, E>13.6 eV) and orbital distance `a`
are read from `input.inp`, `N_HI` is the vertical neutral-H column from
the ATES profile, and `sigma_LyC ~ 6.3e-18 cm^2`. `xi` is the day-night /
2D flux dilution matching ATES's "2D approximate method" (`Rate/2`->0.5,
`Rate/4`->0.25, else 1; read from `input.inp`; cf. the xi factor of
Christie+2013 / Huang+2017). For the optically-thick atomic layer the
absorbed fraction -> 1, so `F_LyC` -> `xi` x incident stellar LyC (e.g.
xi=0.5 -> ~5.2e2 erg cm^-2 s^-1 for the HD209458b input; Huang+2017 found
2.6e4 absorbed for HD189733b). Force values with `F_LyC_override > 0`
and/or `xi_override > 0`.

This is a uniform-illumination approximation — it omits the decline of
`J_lya` deep in the atmosphere (Huang+2017 Eq.7) and the direct stellar
Ly-alpha component; for a self-consistent profile, use Option 1.

Then run `python3 TPM.py` as usual: **H-alpha and H-beta** are added as
extra figures and printouts; He/Ly-alpha behavior is unchanged.

**H-beta (n=2 -> 4, 4861.4 A)** shares the same n=2 population and is
computed automatically alongside H-alpha (knobs `Instr_res_Hb`,
`lmin_Hb/lmax_Hb/number_lambda_Hb`, `fig_name_hb`). Both lines are
optically thick, so their depth ratio is far smaller than the
optical-depth ratio (~7.3).

**n=2 photoionization** `Gamma_2s/Gamma_2p` is estimated from a diluted
stellar-blackbody Balmer continuum when `T_star > 0` (stellar effective
temperature [K]; uses `R_star`, `a`). This reproduces Huang+2017's
~20-26 s^-1 for a K dwarf and is larger for a hot G dwarf (~140 s^-1 for
HD209458), suppressing n2. Set `T_star <= 0` to use the manual
`Gamma_2s/Gamma_2p` (default 0).

Other H-alpha knobs: `Instr_res_Ha`, `lmin_Ha/lmax_Ha/number_lambda_Ha`,
`fig_name_ha`.

A worked **HD209458b comparison with Jensen et al. (2012)** is in the
notebook `HD209458b/Halpha_compare_HD209458b.ipynb`.

## Ly-alpha pumping convention

The pump rate is `B_{1s->2p} * J_lya` with the Einstein coefficient in
the mean-intensity convention,
`B_{1s->2p} = (g_2p/g_1s)(c^2/2 h nu^3) A_{2p->1s} ≈ 8.55e9` (cgs). So
`J_lya` must be the **angle-averaged mean intensity `J_nu`** at line
center in cgs. A different convention requires rescaling the constant.

## Caveats

- **n=2 photoionization** `Gamma_2s/Gamma_2p` is now estimated from a
  single diluted stellar blackbody (Balmer continuum) via `T_star` — an
  approximation to the true stellar near-UV spectrum. `T_star <= 0`
  reverts to the manual constants (default 0).
- `n_e = n_HII + n_HeII + 2 n_HeIII` from ATES (Christie assume
  `n_e = n_p`; difference is small).
- H-alpha uses the air wavelength 6562.8 Å and `f_23 = 0.64`; the line is
  Doppler-dominated (natural width negligible).

## Implementation notes

- `n2_populations(T, n1s, ne, Jlya, G2s, G2p)` returns `(n2s, n2p, n2)`.
- **Compatibility fix:** the `Ion_species.txt` reader now uses
  `usecols=range(7)` (EXHALE's file has 16 columns including the
  trace metals); this is also correct for the old 13- and 7-column files.

## References

- Christie, D., Arras, P., & Li, Z.-Y. 2013, ApJ, 772, 144 —
  *Hα Absorption in Transiting Exoplanet Atmospheres*
  (`references/Christie_2013ApJ_772_144.pdf`). The n=2 level-population
  method.
- Huang, C., Arras, P., Christie, D., & Li, Z.-Y. 2017, ApJ, 851, 150 —
  *A Model of the Hα and Na Transmission Spectrum of HD 189733b*
  (`references/Huang_2017_ApJ_851_150.pdf`). The default `J_lya`
  estimate (Eq. 6 and related text) and the F_LyC / Gamma_2 prescriptions.
- Jensen, A. G., et al. 2012, ApJ, 751, 86 — transit Hα spectroscopy of
  HD 209458b / HD 189733b (comparison target; see
  `HD209458b/Halpha_compare_HD209458b.ipynb`).


## Metal resonance doublets (Mg II h&k, Ca II H&K, Na I D)

The metals absorb out of the ion ground state (excited fine-structure
levels are Boltzmann-negligible at ~1e4 K), so `n_lower = n_ion` read
directly from the metal block of `Ion_species_adv.txt`; a metals-off run
is skipped automatically. Both doublet components (NIST f/A values:
Mg II 2796.352/2803.531, Ca II 3933.663/3968.469, Na I 5889.951/5895.924)
are summed in one wavelength window and pushed through the same pipeline
as He/Ly-alpha: spherical-chord Voigt LOS integration, disk average,
instrument convolution, and planet-rotation convolution with the
depth-dependent effective radius. Figures (PNG + vector PDF) are saved
via `fig_name_mgii/caii/nai`.

The validated single-component depth table (line-center and 4 Å-band %)
is still printed for the Huang+2023 comparisons. Two conventions matter
there: Huang's table is the effective transit radius R_eff/R_star =
sqrt((Rp/R*)^2 + h) (not a percent), and Mg II is quoted in a 4 Å NUV
bin while the optical lines are line-center.

For Roche-lobe runs, `geometry = 'triaxial'` (with `roche_recon.py`)
replaces the spherical chords by the 3-D equipotential reconstruction
with wind + tidally-locked-rotation LOS velocity.

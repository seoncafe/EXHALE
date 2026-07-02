# Version comparison: EXHALE v1.0 vs current (He/H + metal diffusion)

**Purpose.** Quantify what the new diffusive separation changes in EXHALE's predictions, by
running the **same** planet on the two code versions, for **two contrasting planets**:

- **v1.0** — `ATES/EXHALE_v1.0/`: clean GitHub `main` + the temperature-dependent Penning
  rate. **No** diffusive separation (He/H and metal/H frozen at the input values).
- **current** — `ATES/EXHALE/`: v1.0 physics **plus** He/H and per-element metal diffusion
  (`He_diffusion` + `He_metal_diffusion`, ambipolar settling on).

Both versions share the temperature-dependent Penning rate, so the **only** difference is the
diffusive separation. Setup: He 2³S on, trace metals from `metals.inp` (HD 209458b: C/N/O;
WASP-121b: C/N/O/Mg/Ca/Na/Fe), identical `input.inp` per planet, 15000 relaxation steps.
Run dirs / figures / regenerate script: `docs/version_compare/`.

---

## Headline numbers

### HD 209458b (hot Jupiter, log₁₀ Ṁ ≈ 9.3 — gentler escape)

| quantity | v1.0 (no diff) | current (diff) | change |
|---|---|---|---|
| He/H at 3 R_p | 0.083 (flat) | **0.014** | −83% |
| peak n(He 2³S) [cm⁻³] | 177 | 102 | −42% |
| **He I 10830 line-center absorption** | **70.5%** | **27.1%** | **−2.6×** |
| log₁₀ Ṁ [g s⁻¹] | 9.31 | 9.04 | −0.27 dex |

### WASP-121b (ultra-hot Jupiter, log₁₀ Ṁ ≈ 13.3 — furious escape)

| quantity | v1.0 (no diff) | current (diff) | change |
|---|---|---|---|
| He/H (min, 1.05–1.5 R_p) | 0.083 | **0.053** | −36% (weaker) |
| peak n(He 2³S) [cm⁻³] | 3517 | 2383 | −32% |
| He I 10830 line-center absorption | 53.7% | 53.6% | ~unchanged |
| metal/H (element/base) | 1.0 | tracks He (~0.6–1.0; no extra fractionation) | none |
| log₁₀ Ṁ [g s⁻¹] | 13.34 | 13.32 | unchanged |

**Two regimes.** On the gentler HD 209458b the escape is slow enough that *helium itself*
diffusively separates, cutting the He I 10830 line by ~2.6×, and the heavier C/N/O deplete
more than He aloft (mass ordering He > C > N > O at 3 R_p). On the furious WASP-121b
(~4 dex stronger wind) advection dominates settling for **every** species — helium *and* the
metals are dragged out essentially unfractionated (10830 unchanged; metal/H tracks He/H).
This is the expected (settling speed)/(wind speed) scaling of Koskinen et al. (2013) / Xing
et al. (2023): separation is suppressed by fast escape and — where it operates at all — is
stronger for heavier species.

*(Correction note: an earlier version of this comparison showed WASP-121b metals collapsing
to zero by ~1.25 R_p, a "metal homopause." That was an artifact of two bugs in the metal
rescale — a `rX ≤ 1` clamp that made depletion a one-way ratchet, and a skip of exhausted
cells that made holes permanent — which locked in transient early-relaxation settling before
the wind developed. With the fix, metals recover with the wind and track He, and the v1-vs-v2
temperature structures agree; the earlier "~2500 K hotter" claim is also retracted.)*

---

## HD 209458b figures

![HD 209458b profiles](version_compare/fig_hd209_profiles.pdf)

**Figure 1.** HD 209458b, v1.0 (blue dashed) vs current (red). He/H falls to ~0.15–0.3× aloft
with diffusion; the metastable He 2³S reservoir contracts, weakening the line. Bottom-right:
C/N/O diffuse independently and, being heavier than He, deplete more aloft (mass ordering
He > C > N > O above ~2.5 R_p).

![HD 209458b He I 10830](version_compare/fig_hd209_He10830.pdf)

**Figure 2.** He I 10830 transit line. Diffusive separation lowers the line-center absorption
from 70.5% to 27.1% and narrows the profile.

## WASP-121b figures

![WASP-121b profiles](version_compare/fig_wasp_profiles.pdf)

**Figure 3.** WASP-121b (compact domain, escape radius 1.5 R_p). He/H separates only modestly
(~0.6–0.75× in the inner thermosphere, returning to ~1 above 1.4 R_p) because the
~4-dex-stronger wind drags helium out. The metals (O, Mg, Fe) overlap He almost exactly —
even Fe (mass 56) is advection-dominated here, so no additional fractionation develops. The
v1/v2 temperature and velocity structures are nearly identical.

![WASP-121b He I 10830](version_compare/fig_wasp_He10830.pdf)

**Figure 4.** He I 10830 line: the core is essentially unchanged (53.7%→53.6%) since helium is
not strongly separated in the line-forming region.

---

## Interpretation

- **He/H.** v1.0 carries the input He/H at every radius; the current version transports the
  He element with advection + molecular-diffusion settling. On HD 209458b settling competes
  with the (slower) wind, so He/H declines with altitude and the 10830 line weakens. On
  WASP-121b the wind is ~4 dex stronger, so helium is dragged out nearly unfractionated.
- **He 2³S / 10830.** The 2³S metastable tracks the He abundance; where He is depleted the
  reservoir shrinks and the line weakens (HD 209458b), where He stays mixed the line is
  preserved (WASP-121b).
- **Metals.** Each trace metal diffuses independently against H with its own mass, binary
  diffusion coefficient, and ambipolar-corrected settling. Where diffusion competes with the
  wind (HD 209458b aloft) heavier elements deplete faster (He > C > N > O); where advection
  dominates (WASP-121b) even Fe stays locked to He. Relevant for metal-line transmission
  diagnostics (Mg II / Ca II / Na I in `TPM.py`).

## Caveats

- **Absolute values are uncalibrated.** The 27–70% He 10830 depths far exceed observations
  (~1–2%): heating efficiency, XUV level, and He 2³S microphysics were not tuned to data. The
  **relative** v1.0-vs-current change (and the HD 209458b-vs-WASP-121b contrast) is the robust
  result, not the absolute depth.
- **Scheme / convergence sensitivity.** Both runs are step-capped (15000 steps, not
  Newton-finished), so Ṁ is indicative. The exact He-separation magnitude is somewhat
  sensitive to the settling discretization (a Peclet-hybrid central/upwind scheme is used —
  central where well-resolved, upwind for the stiff heavy-metal settling); the qualitative
  trends are robust.
- **Phase-1/2 numerics.** The He/H cap and the not-fully-developed wind leave mild
  non-monotonic structure; see `design_hehe_diffusion.md` §7d–§7e. Diffusion flags are
  **default OFF**, so standard runs are unaffected.

---

*Regenerate:* re-run the cases in `docs/version_compare/{v1_nodiff,v2_diff,wasp_v1_nodiff,
wasp_v2_diff}/`, run `TPM.py` (TPM_PATH / TPM_SAVE_PREFIX) in each — for the WASP-121b runs
add `TPM_HE_LMIN=10827.5 TPM_HE_LMAX=10832.5 TPM_HE_N=335` (the broad line overflows the
default 10828.2–10831.2 Å window) — then `python3 docs/version_compare/plot_compare.py`.

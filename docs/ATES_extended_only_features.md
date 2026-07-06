# Features in ATES_extended but NOT in EXHALE

This lists everything present in `ATES/ATES_extended/` that has no
counterpart in `ATES/EXHALE/` (the latter being this session's
fully-coupled metal-cooling work). It is the "what would be lost if we
discarded ATES_extended" inventory.

See also `ATES_versions_diff.md` for the full side-by-side comparison.

---

## 1. Opacity model dispatcher (entire Phase 1 — no equivalent at all)

> **Update (2026):** No longer extended-only. The dispatcher
> (`opacity_models.f90`, A/C/P/T, `.atesopa` tables, `opacity.inp`) has
> been ported to EXHALE, and the Robinson & Catling
> pressure-broadening — a dormant hook in ATES_extended — is now
> **actually applied per-cell** (weights the opacity column density). See
> `Update_EXHALE_early_phase`, Part II. The description below is the original
> ATES_extended implementation.

ATES_extended adds a complete pluggable opacity layer:

* **`src/modules/radiation/opacity_models.f90`** (~232 LOC) — dispatcher
  `photoion_sigma` selecting among four models:
  * `'A'` analytic (current ATES default)
  * `'C'` constant (threshold cross section x user factor)
  * `'P'` physical (constant + Robinson & Catling pressure broadening)
  * `'T'` tabulated (per-species `.atesopa` file)
  Returns in ATES's internal `1e-18 cm^2` unit so downstream code is
  unchanged. Includes the `opacity_pT_factor(p)` per-cell multiplier and
  `opa_table` derived type with `load_opacity_tables` / `free_opacity_tables`.
* **`src/modules/files_IO/opacity_input_read.f90`** (~116 LOC) — parser
  for the new `opacity.inp` key=value file.
* **`.atesopa` table format** — two-column (E[eV], sigma[1e-18 cm^2])
  per-species tabulated cross-section file, with analytic fallback when a
  path is empty.

## 2. Nitrogen as a coolant/ion species

> **Update (2026):** No longer extended-only. EXHALE now also
> carries nitrogen (NI/NII/NIII) inside its MINPACK 9-equation system,
> with Verner+1996 cross sections, Badnell RR+DR recombination, and
> Voronov collisional ionization; like ATES_extended, N contributes no
> line cooling. See `Update_EXHALE_early_phase`, Part II. The description below is
> the original ATES_extended implementation.

ATES_extended carries:

* **NI, NII** cross sections (Verner+1996) in `cross_sec_metals.f90`.
* **NI, NII recombination** fit blocks (`rec_NI`, `rec_NII`) in
  `metals_solve.f90`.
* NI/NII are recognized ions in `metals.inp` (solar N/H = 6.76e-5 in the
  example file).
* Note: their *cooling* coefficient `lambda_X` returns 0 for N (no AIOLOS
  reference); the memo flags Hollenbach & McKee 1989 fits as the upgrade
  path. So N is fully wired for ionization/opacity but not yet for cooling.

## 3. Dielectronic recombination (DR)

> **Update (2026):** No longer extended-only. EXHALE now uses the
> **Badnell 2006 RR + adf48 DR** total-recombination fits for all metal
> stages (the same `rec_fit` table, ported from ATES_extended), replacing
> the Aldrovandi & Pequignot 1973 power-law. See `Update_EXHALE_early_phase`,
> Part II. The description below is the original ATES_extended implementation.

ATES_extended adds a second recombination channel beyond radiative RR:

* **Shull & Van Steenberg (1982) Burgess-form DR**:
  `alpha_DR(T) = A_DR * T^(-3/2) * exp(-T0/T) * [1 + B_DR*exp(-T1/T)]`,
  summed with RR as `alpha_rec = alpha_RR + alpha_DR`.
* Implemented as `alpha_dr_metal()` with 4-parameter `rec_fit` blocks for
  all six ions (CI, CII, NI, NII, OI, OII).

## 4. Runtime-configurable abundances (no recompile)

> **Update (2026):** No longer extended-only. EXHALE now also
> reads a runtime **`metals.inp`** (`metals_input_read.f90`): `CI/NI/OI`
> lines set `X_C/X_N/X_O` with no recompile; absent file leaves metals
> off. See `Update_EXHALE_early_phase`, Part II.

* ATES_extended reads abundances at runtime from an optional
  **`metals.inp`** file (`<ION> <abundance>` per line), parsed by
  `read_metals_input` in `metals.f90`. Absent file => `n_metals = 0` =>
  all metal code paths no-op. Supports up to `max_metals = 8` ions via a
  `metal_ion` registry derived type.

## 5. Separate metal solver architecture (coronal balance)

The whole module set below has no analog in EXHALE, which instead
folded metals into the existing MINPACK system (`System_HeHCO.f90`):

* **`src/modules/radiation/metals.f90`** (~162 LOC) — abundance registry,
  per-ion cross-section grids on `e_v(Nl)`, `photoion_rate_metal`.
* **`src/modules/radiation/metals_solve.f90`** (~209 LOC) — `coronal_ratio`
  helper solving `n_X^(k+1)/n_X^k = Gamma_k / (alpha_rec * n_e)` per ion
  pair independently, plus the RR/DR coefficient functions.
* **`src/modules/radiation/metals_drive.f90`** (~92 LOC) — per-cell driver
  `solve_metals_post` invoked from post-processing.
* **`src/modules/radiation/metals_cool.f90`** (~102 LOC) — `lambda_X`
  cooling fits and `eval_metal_cooling`; couples to `T_equation` via the
  `paramsT(13)` channel.

## 6. Dedicated metal output file

* **`src/modules/files_IO/write_metals_output.f90`** (~62 LOC) writes a
  separate **`Metals_ioniz_adv.txt`**.
* EXHALE instead appends six columns to the existing
  `Ion_species.txt` (no new file, no separate writer).

## 7. Example/documentation input files

In `ATES_extended/inputdata/`, none of which exist in EXHALE:

* `metals.inp.example` — annotated abundance template.
* `opacity.inp.example` — annotated opacity-model key documentation.
* `HI_sample.atesopa` — 13-row sample tabulated cross-section table.
* `README.opacity` — `.atesopa` / opacity format documentation.

## 8. Companion design memo

* **`ATES_extended/docs/aiolos_port_memo.pdf`** (and `.tex`) — Phase 1/2/3
  design and status document. EXHALE's analog is
  `EXHALE/docs/Update_EXHALE_early_phase` (Part II), but the two cover
  different designs; the ATES_extended memo additionally discusses the
  opacity dispatcher and the explicit Phase 3 (multi-fluid) justification,
  neither of which appears in the EXHALE work.

---

## Summary table

| Capability | In ATES_extended | In EXHALE |
|---|---|---|
| Opacity model dispatcher (A/C/P/T) | Yes | **Yes** (ported) |
| `.atesopa` tabulated cross sections | Yes | **Yes** |
| Pressure-broadening (Robinson & Catling) | Yes (Phase 1; per-cell **deferred**) | **Yes — applied per-cell** (completed) |
| Nitrogen (NI/NII/NIII) ionization + opacity | Yes | **Yes** (in the MINPACK 9-eq system) |
| Dielectronic recombination | Yes (Badnell, despite stale SVS-1982 docstring) | **Yes** (Badnell 2006 RR + adf48 DR) |
| Runtime abundance file (`metals.inp`) | Yes | **Yes** (CI/NI/OI → X_C/X_N/X_O) |
| `opacity.inp` runtime config | Yes | **Yes** |
| Separate `Metals_ioniz_adv.txt` output | Yes | No (appends to `Ion_species.txt`) |
| Coronal-balance metal solver | Yes | No (uses coupled MINPACK system) |
| **Charge transfer with H** (Kingdon & Ferland 1996) | **No** | **Yes** (couples metal & H ionization) |
| **H-alpha transmission** (`EXHALE_transit.py`, Christie+2013 Lyα pumping) | **No** (TPM has He 10830 + Lyα only) | **Yes** (n=2 population + J_lya input; see `transmission_spectrum.{tex,md}`) |

## Caveat — these are *additive* features, not strict supersets

ATES_extended is NOT a superset of EXHALE. The fully-coupled
physics in EXHALE (metals feeding back on electron density,
a beta line-escape probability formula carried in-source — though
currently overridden by `beta_esc = 1`, i.e. the same 100%-escape
assumption ATES_extended uses — CIII/OIII third stages, metals active
during time evolution rather than only post-processing) is absent from
ATES_extended. The two trees are complementary; see
`ATES_versions_diff.md` Section 8.

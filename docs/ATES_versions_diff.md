# ATES_extended vs. ATES-metal: Metal Cooling Implementation Comparison

Two independent attempts at adding AIOLOS-style metal cooling to ATES coexist
in this workspace. They differ in design philosophy, coupling strength,
data sources, and footprint. This document compares them side by side.

* `ATES/ATES_extended/` — earlier work (last edit 2026-04-29), accompanied
  by `docs/aiolos_port_memo.pdf` describing Phase 1 (opacity dispatcher)
  and Phase 2 (trace metals).
* `ATES/ATES-metal/` — recent work (this session, 2026-05-28),
  documented in `docs/Update_ATES_early_phase` (Part II).

---

## 1. Coupling Strategy

| Aspect | ATES_extended | ATES-metal |
|---|---|---|
| Call site | Inside `post_process_adv` only (`solve_metals_post`) | Inside `ionization_equilibrium` every step (`ion_system_HeHCO`) |
| Equilibrium method | **Coronal equilibrium** — each metal solved independently as `Gamma * n_X = alpha * n_e * n_(X+1)` via `coronal_ratio()` | **Full coupled 9-equation MINPACK system** — HII, HeII, HeIII, CII, CIII, OII, OIII, **NII, NIII** solved simultaneously |
| Metal feedback on electron density | None (n_e from H/He only; metals assumed trace) | Included: `n_e = n_HII + n_HeII + 2*n_HeIII + n_CII + 2*n_CIII + n_OII + 2*n_OIII` |
| `f_sp` array dimension | Unchanged (metals carried as separate state) | Extended 6 → 15 (H/He 1-6, C 7-9, O 10-12, **N 13-15**) |

---

## 2. Data Sources

| Aspect | ATES_extended | ATES-metal |
|---|---|---|
| Cross sections | Verner+1996 for CI, CII, NI, NII, OI, OII | Verner+1996 for CI, CII, **NI, NII**, OI, OII |
| Recombination | Badnell 2006 RR **+ Badnell adf48 DR** (total fits; the module docstring's mention of AP1973/SVS1982 is stale) | **Badnell 2006 RR + adf48 DR** (total fits, `rec_fit` table in `Cool_coeff.f90`); modernized from the original Aldrovandi & Pequignot 1973 power-law (RR only) |
| Cooling rates | **Black 1981 fits** (`lambda_X`) + two-level fine-structure ([OI]63, [CII]158) | **Black/AIOLOS fits** (identical formulas) + optional two-level fine-structure behind `use_2lev_cool` (default off) |
| Escape probability beta | **None** (optically thin assumed) | AIOLOS's beta formula retained but commented out (tagged `!To Be Checked/AIOLOS tuning?`); production code now sets `beta_esc = 1` (100% escape). AIOLOS multiplies the escape probability by a `1e8` boost, driving beta→~0 (lines fully trapped, cooling off); we instead assume 100% escape, matching ATES_extended |
| Elements covered | C, N, O | **C, N, O** |
| Ionization stages | CI/CII/CIII, NI/NII/NIII, OI/OII/OIII (coronal) | CI/CII/CIII, **NI/NII/NIII**, OI/OII/OIII (MINPACK 9-eq) |

---

## 3. Input / Output

| Aspect | ATES_extended | ATES-metal |
|---|---|---|
| Input files | New `metals.inp` (ion abundances), new `opacity.inp` (opacity-model selection) | **Same runtime `metals.inp` (CI/NI/OI → X_C/X_N/X_O) and `opacity.inp`** now supported (no recompile); `input_read.f90` defaults to metals-off |
| Output files | New `Metals_ioniz_adv.txt` (separate file) | Six metal columns appended to existing `Ion_species.txt` |

---

## 4. Auxiliary Features

| Aspect | ATES_extended | ATES-metal |
|---|---|---|
| Opacity model dispatcher (Phase 1) | **Yes** — A/C/P/T models in `opacity_models.f90`, `.atesopa` table format, pressure-broadening **hook (deferred)** | **Yes** — same dispatcher ported; pressure broadening **applied per-cell** (completed) |
| Charge transfer with H (Kingdon & Ferland 1996) | No | **Yes** — O/N/C ↔ H in `System_HeHCO.f90`, couples metal & H ionization |
| Module organization | Highly modular (8 new files plus 5 modifications) | Less modular (1 new file plus 14 modifications) |
| Regression test specification | Documented in memo: with no `metals.inp` present, results must be bit-identical to upstream | Only a smoke test (executable initializes correctly) |

---

## 5. Code Footprint

| Metric | ATES_extended | ATES-metal |
|---|---|---|
| New code | ~1,200 LOC across 12 new files | ~250 LOC in 1 new file (`System_HeHCO.f90`) |
| Modified existing code | ~150 LOC across 5 files | ~450 LOC across 14 files |
| Total | ~1,350 LOC | ~700 LOC |

---

## 6. Design Trade-offs

**ATES_extended** (conservative, modular):

* + Leaves the existing ATES core (H/He solver, `f_sp` array) untouched —
  low regression risk.
* + The opacity dispatcher (Phase 1) is independently useful even without
  the metals work.
* + Black 1981 cooling fits have a clear published reference.
* - Metals cannot feed back on the electron density (coronal approximation
  limitation).
* - Metals participate only in post-processing, not during time evolution.

**ATES-metal** (aggressive, fully coupled):

* + Metals contribute self-consistently to equilibrium, energy balance,
  and opacity — faithfully reproduces AIOLOS behavior.
* + The beta escape probability formula for line trapping is available
  in-source (currently overridden by `beta_esc = 1`, the 100%-escape /
  optically-thin assumption; AIOLOS's `1e8` boost is commented out).
* - Extending `f_sp` from 6 to 15 columns touches many files; higher
  regression risk.
* - AIOLOS analytic cooling fits have unclear provenance (see
  `ATES-metal/docs/Update_ATES_early_phase`, Part II).

---

## 7. Outstanding Work

### From the ATES_extended memo, "Recommended next steps":

Code work still pending in ATES_extended (testing and validation items
omitted):

1. **Time-loop integration** — currently the metals solver runs only in
   `post_process_adv.f90`; should be called from
   `ionization_equilibrium.f90` every step. (Most impactful unfinished item.)
2. **CIII / NIII / OIII** — extend the Verner table; generalize
   `coronal_ratio` to a three-element system.
3. **Recombination upgrade** — replace A&P 1973 RR parameter blocks with
   Verner & Ferland 1996 4-parameter values (~5x improvement for CI, OI);
   replace SVS 1982 DR with the Badnell 2003+ DR series. Local edits to
   parameter blocks; all callers remain unchanged.
4. **Cooling fidelity** — add the n_HI collider channel and the
   critical-density saturation factor to `lambda_X`; add NI/NII fits.
5. **Phase 1d** — pressure-broadening multiplier per cell. The current
   Phase-1 implementation only replaces the prebaked global vector;
   per-cell application is an explicit deferred item.
6. **Python plot helper** — extend `ATES_plots.py` to overlay
   `Metals_ioniz_adv.txt` on the same coordinates as
   `Hydro_ioniz_adv.txt`.

### From ATES-metal:

See `ATES/ATES-metal/docs/Update_ATES_early_phase`, Part II
("Caveats and Verification Items") for the corresponding list. The most
important verification items are:

1. Verner cross-section parameters cross-checked against Table 1 of
   Verner+1996.
2. Voronov ionization parameters cross-checked against Table 1 of
   Voronov 1997.
3. AIOLOS cooling fit provenance (likely SB23 paper or Schulik 2022
   thesis).
4. `post_process_adv.f90` currently passes zero metal density vectors —
   self-consistent post-processing with metals is deferred.
5. HeITR and metals are presently mutually exclusive in the dispatch
   logic.

---

## 8. Conclusion

The two implementations are **complementary rather than mutually
exclusive**.

* **For physical accuracy**, ATES-metal's full coupling (metals in
  the MINPACK system, electron-density feedback) is closer to AIOLOS's
  actual behavior. Note that both trees now apply the metal-line cooling
  with **100% escape** (`beta = 1`): ATES-metal overrides AIOLOS's
  `1e8`-boosted beta, and ATES_extended never carried one. With AIOLOS's
  `1e8` factor the metal cooling would instead be almost fully
  suppressed.
* **For maintainability and verifiability**, ATES_extended's modular
  layout (with regression testing built in) and clear data citations
  (Black 1981, Verner+1996) are stronger.

A reasonable next step is to merge the strongest elements of each into a
single tree: keep ATES_extended's opacity dispatcher (Phase 1) and the
Black 1981 cooling fits (which include nitrogen and have a clear
reference), then implement next-step #7 ("time-loop integration") in the
ATES-metal style (full coupling with the MINPACK 9-equation system).
This requires a decision on which branch to designate as the active
trunk.

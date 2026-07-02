# EXHALE Update Log

A running log of incremental updates to `EXHALE` as the code is extended
toward reproducing Huang et al. (2023) WASP-121b (see `Huang_update_plan.md`).
Each update is a self-contained section: what changed, why, which files were
touched, and the quantitative check that was run. New updates are appended as
new top-level sections. Full typeset version with equations and figures:
`Update_EXHALE.{tex,pdf}`.

---

## 1. Runtime domain-extent option: Roche vs Spherical

*Added 2026-06-04.*

A runtime option that selects how far the radial domain extends and which
gravitational potential is used, chosen from `input.inp` without recompiling.
It closes the Phase-0 geometry gap toward Huang et al. (2023) WASP-121b: the
paper's spherical Case A extends to ~19 R_p, whereas the standard ATES domain
truncates at the Roche (Hill/L1) radius (~2 R_p for WASP-121b).

### The two modes

Units are code-normalized: radius `r` in planetary radii `R0 = R_p`;
`b0 = G M_p mu / (k_B T0 R0)` is the surface Jeans parameter;
`M_r = M_star / M_p`; `atilde = a_orb / R0`.

- **Roche** (default): full Roche potential
  `phi(r) = -b0/r - b0*M_r/(atilde - r) - b0*(1+M_r)/(2*atilde^3) * (atilde*M_r/(1+M_r) - r)^2`
  (planetary gravity + stellar tidal + centrifugal), domain truncated at the
  Hill radius `r_max = (3*M_r)^(-1/3) * atilde`. Identical to previous ATES
  behavior.
- **Spherical**: pure planetary potential `phi(r) = -b0/r` (tidal and
  centrifugal terms dropped), domain extended to a user-set outer radius.
  Matches Huang Case A geometry.

**Why spherical drops the tidal terms.** Beyond L1 the full Roche potential
turns over (effective gravity reverses toward the star), so a steady outward
wind is undefined there and the sonic-point machinery breaks down. Dropping the
tidal/centrifugal terms restores a monotonic potential, which is what allows the
domain to extend arbitrarily far.

### Input interface

Append to `input.inp` after the "Force start:" line (omit for default Roche):

```
Domain mode: Spherical
Outer radius [R_p]: 19.0
```

- "Domain mode:" accepts `Spherical`; any other value or a missing line means
  `Roche`.
- "Outer radius [R_p]:" is read only in Spherical mode (R_p units, since the
  grid is normalized to R0 = R_p); ignored in Roche mode.

The trailing lines are scanned by keyword (`index`), not by fixed position, so
blank lines and ordering do not matter and older `input.inp` files (lacking both
lines) keep the default Roche behavior. Selecting Spherical without a valid
outer radius (> 1) stops the run with an explicit error. The active mode is
echoed in `ATES.out`.

### Touched files

- `src/modules/init/parameters.f90` — globals `spherical_domain` (logical) and
  `r_out_user` (outer radius in R_p).
- `src/modules/functions/grav_field.f90` — `phi` and `Dphi` branch on the mode.
- `src/modules/files_IO/input_read.f90` — keyword scan of trailing lines and the
  `r_max` branch with validity guard.
- `src/modules/files_IO/write_setup_report.f90` — prints the active mode.

### Validation on WASP-121b

Two runs differing only in this option (otherwise identical: power-law SED,
He/H = 0.0851, C/N/O metals on, "Rate/2 + Mdot/2"). Both reached steady state
via the proper criterion (`du < du_th`), not a stall or false convergence.

| mode | r_max [R_p] | log10 Mdot | sonic [R_p] | r(50% HII) | T_peak [K] |
|---|---:|---:|---:|---:|---:|
| Roche (default) | 2.01 | 12.68 | 1.825 | 1.178 | 9662 |
| Spherical | 19.00 | 11.81 | 5.270 | 1.152 | 10752 |
| Huang Case A (ref.) | ~19 | 12.57 | — | — | ~12000 |

Mdot in log10 g/s. Roche took ~5,392 iterations, Spherical ~154,183.

Observations (for evaluation, not pass/fail):

- Roche reproduces the prior Phase-0 Case A exactly (12.68) — the optional block
  does not perturb default behavior.
- Spherical T_peak (10,752 K) is higher than Roche (9,662 K) and closer to
  Huang's ~12,000 K — the extended domain develops the profile past the
  2-R_p L1 clip.
- Spherical Mdot (11.81) is below both Roche (12.68) and Huang (12.57),
  consistent with loss of tidal enhancement, but ~0.76 dex under Huang.
- The sonic point moves out from 1.83 to 5.27 R_p.

Caveats: two confounders remain before the spherical Mdot/T can be compared to
Huang cleanly — (a) only C/N/O metals here, not Huang's full
Mg/Fe/Si/Ca/Na/K/S set (Phase 1); and (b) the same N=500 cells now span
1 -> 19 R_p, so inner resolution is much coarser than the Roche run.

Scratch runs: `WASP-121b/domain_test/{roche,spherical}/`; comparison notebook
`WASP-121b/domain_test/domain_compare.ipynb`.

---

## 2. Ionization network: add magnesium (Mg I / II / III)

*Added 2026-06-05.*

The first trace metal beyond C/N/O, following the Huang et al. (2023) WASP-121b
species set (Mg, Fe, Si, O, C, N, S, Ca, Na, K). Per the Phase-1 plan, elements
are added one at a time (order Mg -> Fe -> Si -> Ca -> Na -> K -> S); this update
covers Mg only. Mg is solved to second ionization (Mg I / Mg II / Mg III) inside
the existing coupled MINPACK ionization system, raising the species count to 18
and the ratio-unknown count to 11.

### Atomic data and rates (sources, not fabricated)

- **Photoionization cross sections:** Verner et al. (1996) outer-shell fits.
  Mg I threshold 7.646 eV, Mg II threshold 15.035 eV (`cross_sec.f90`:
  `sigma_MgI`, `sigma_MgII`).
- **Recombination:** Badnell RR+DR (the same source already used for C/N/O),
  via `lookup_rec` cases `MgI`/`MgII` and `rec_MgII` (-> Mg I) / `rec_MgIII`
  (-> Mg II) in `Cool_coeff.f90`. **This is a deliberate deviation** from the
  validation-rate prescription in `Huang_update_plan.md` §1.3 (Shull & van
  Steenberg 1982 for Mg). It corresponds to the plan's Phase-6 *production*
  rates, applied early to match the existing C/N/O treatment. The Mg ion/neutral
  balance below should therefore be read as the Badnell result, not a
  Huang-faithful one.
- **Collisional ionization:** Voronov (1997) fits, `ion_coeff_MgI` (dE = 7.646,
  A = 6.21e-7, P = 0, X = 0.592, K = 0.39) and `ion_coeff_MgII` (dE = 15.035,
  A = 1.92e-8, P = 0, X = 0.0027, K = 0.85).
- **Line cooling:** **not yet ported** -- `cool_MgI = cool_MgII = 0`. No
  published Mg two-level/CHIANTI cooling fit is implemented (that is Phase 2).
  Magnesium therefore acts only as an electron donor and a photoelectric-heating
  channel here, not as a coolant.
- **Charge exchange Mg + H+:** omitted in this first cut.

### Grid extension below 13.6 eV

Because Mg I ionizes at 7.646 eV, turning Mg on triggers the existing
low-IP-metal branch (`thereis_lowIP_metal`), which extends the energy grid and
the stellar SED down to 7.646 eV (`set_energy_vectors.f90`). This is required to
photoionize Mg I at all, but it has a side effect described under Caveats.

### Touched files

Ionization/radiation plumbing only -- no hydro, mass, or `calc_ne`/`calc_rho`
changes (metal electrons enter solely inside the MINPACK system):

- `init/parameters.f90` -- `n_species = 18`; `e_th_MgI`/`e_th_MgII`; `X_Mg`;
  `s_mgi`/`s_mgii` cross-section arrays; `thereis_lowIP_metal`.
- `functions/cross_sec.f90` -- `sigma_MgI`, `sigma_MgII`.
- `radiation/Cool_coeff.f90` -- Badnell fits, `lookup_rec`, `rec_MgII`/
  `rec_MgIII`, Voronov `ion_coeff_MgI`/`MgII`, `cool_MgI`/`MgII` (= 0, flagged).
- `radiation/util_ion_eq.f90` -- `PH_heat_HHe` and `eval_cool` carry Mg.
- `radiation/ionization_equilibrium.f90` -- Mg densities, params slots 49-55,
  MINPACK unknowns 10/11 (Mg II/Mg III ratios), solution unpack.
- `init/set_energy_vectors.f90` -- 7.646 eV sub-grid; allocate/fill `s_mgi`/
  `s_mgii`.
- `files_IO/input_read.f90` -- `X_Mg` default; `N_eq = 11`; metals/low-IP flags.
- `files_IO/metals_input_read.f90` -- `MgI` abundance key.
- `files_IO/write_output.f90` -- Mg I/II/III as `Ion_species.txt` columns 17-19.
- `init/set_IC.f90`, `files_IO/load_IC.f90`, `init/init.f90`, `EXHALE_main.f90`,
  `time_step/energy_semi_implicit.f90`, `functions/utilities.f90`,
  `post_process/post_process_adv.f90` -- thread Mg through arrays, I/O, and
  post-processing. (`post_process_adv.f90` also got a latent-bug fix: its
  advection systems fill only the H/He equations, so they now use a local
  `Neq_adv` instead of the global `N_eq = 11`, which would feed `hybrd1`
  uninitialized metal entries.)

### Validation on WASP-121b

Four runs, identical `input.inp` (power-law SED, He/H = 0.0851, Roche domain),
differing only in which metals are present. All reached steady state via
`du < du_th` (no stall). Mg/H is conserved at the input X_Mg = 3.98e-5 across the
domain.

| run | metals | grid floor | log10 Mdot | T_peak [K] | r(T_peak) [R_p] |
|---|---|---:|---:|---:|---:|
| no metals | H/He only | 13.6 eV | 12.92 | 13722 | 1.703 |
| Mg only | Mg | 7.646 eV | 12.92 | 13682 | 1.697 |
| C/N/O | C, N, O | 13.6 eV | 12.68 | 9662 | 1.199 |
| C/N/O + Mg | C, N, O, Mg | 7.646 eV | 12.76 | 9828 | 1.224 |
| Huang Case A (ref.) | full set | — | 12.57 | ~12000 | — |

Mdot in log10 g/s.

Observations (for evaluation, not pass/fail):

- **Clean Mg isolation (Mg only vs no metals):** dlog10(Mdot) = +0.000 dex,
  dT_peak = -40 K. Adding magnesium alone leaves the wind essentially unchanged,
  as expected for a 4e-5 trace species with no line cooling: it neither heats nor
  cools appreciably, and the soft-UV band it adds is optically thin in the
  absence of other absorbers.
- **Mg ionization structure:** Mg is predominantly Mg II (Mg+) through the lower
  thermosphere and predominantly Mg III (Mg2+) in the upper thermosphere (Mg III
  fraction ~0.47 at 1.4 R_p rising to ~0.81 at 1.9 R_p); neutral Mg I stays
  subdominant because of its low 7.646 eV threshold. This is the same qualitative
  ordering as Huang+2023 Fig. 12.
- **MINPACK convergence:** the 11-unknown system converges without NaNs across
  the full temperature range in all runs.
- Mg-only Mdot (12.92) is ~2.2x the Huang Case A spherical value (12.57), in the
  same ballpark as the prior C/N/O Roche run; the metal set and geometry
  confounds from Phase 0/the gate still apply.

### Caveats and gaps

- **C I grid-extension confound (C/N/O + Mg vs C/N/O).** The Mdot/T offset
  between these two runs (+0.080 dex, +165 K) is *not* a magnesium effect.
  Extending the grid to 7.646 eV exposes **C I** (threshold 11.26 eV) to the new
  11.26-13.6 eV photons that the 13.6 eV-floor C/N/O run never sees, adding C I
  photoionization heating. This is why the clean Mg measurement uses Mg-only vs
  no-metals, where C is absent.
- **No Mg line cooling** (`cool_MgI = cool_MgII = 0`). Huang's headline Mg result
  -- Mg II h&k emission dominating the cooling at 1.15-1.4 R_p -- cannot appear
  here by construction. That is Phase 2.
- **No Mg + H+ charge exchange** yet.
- **Recombination is Badnell, not the §1.3 validation rate** (see above).

### Status and relation to the Huang plan

This completes the **Mg** step of Phase 1 (network expansion). The ionization
structure satisfies the qualitative Fig. 12 ordering and the "Mdot within a
factor ~2 of Case A" part of the gate; the thermal-balance part of Phase 1 is
deferred to Phase 2 (Mg cooling). Next element in the planned order: **Fe**.

Scratch runs: `mg_validation/{nometals,mg_only,cno,cno_mg}/`; comparison notebook
`mg_validation/Phase1_Mg_validation.ipynb`.

## 3. Phase 1a: interface / data-structure refactor (no new physics)

*Added 2026-06-05.*

Each metal added so far (C/N/O, then Mg in §2) widened the named-argument lists
of `PH_heat_HHe`, `eval_cool`, `calc_column_dens_metals`, `write_output`, the
MINPACK system `ion_system_HeH_metals`, and the driver `ionization_equilibrium.f90`
by several per-ion arguments, and threaded the per-cell MINPACK coefficients
through a hand-indexed flat `params(60)`. This step converts that plumbing to a
**species-metadata table + array-indexed (2D) data flow** so that adding an
element becomes "add rows to a table" instead of "thread N new arguments through
M routines." **No physics changes** -- the network is still C/N/O + Mg, and the
result is gated bit-for-bit against §2 (below).

### What changed (Steps A-F)

- **A -- Metadata module** (`init/species_table.f90`, new). One canonical table
  holds per-ion metadata (length `n_mion = 12`): f_sp column (`mion_fsp`), parent
  element (`mion_elem`), stage / charge^2 (`mion_stage`, `mion_z2`),
  photo-ionizable flag and photo-table column (`mion_isphot`, `mion_iphot`),
  threshold (`mion_ethr`), cooling flag (`mion_iscool`); and per-element metadata
  (length `n_melem = 4`): nuclear charge (`melem_Z`), neutral-ion index
  (`melem_i0`), top stage (`melem_top`). The canonical ion order
  (CI..MgIII = `f_sp` columns 7-18) is defined here once. Element indices
  `iel_C/O/N/Mg`.
- **B -- `write_output` 2D.** Metal ion densities pass as `nm(:, 1:n_mion)` and
  are written with an implied-do over the table (`(nm(j,i)*n0, i = 1,n_mion)`),
  not per-ion named writes.
- **C -- `PH_heat_HHe` + column density 2D.** Optical depth, photoheating, and
  the photoionization rates loop over photo-ionizable ions through
  `sigma_tab(:,mion_iphot(i))` and `Nm_col`, returning a 2D rate array
  `P_m(:, 1:n_mion)` (inert top stages stay zero) instead of per-ion `P_CI..P_MgII`
  arguments.
- **D -- `eval_cool` metal channels 2D.** Recombination and collisional-ionization
  coefficients return as 2D arrays (`rec_m`, `aion_m`); bremsstrahlung sums
  `mion_z2(i) * GF_elem(:,mion_elem(i)) * nm(:,i)` over ions. No per-ion
  coefficient arguments.
- **E -- Generalized MINPACK system.** `ion_system_HeH_metals` keeps hybrd1's fixed
  `(N_eq, x, fvec, iflag, params)` signature but (i) receives the per-cell metal
  coefficients through a **module-level block** in `System_HeH_metals`
  (`set_metal_coeffs` stores `met_ntot/g0/g1/b0/b1/a1/a2`), set by the driver before
  each `hybrd1` call -- safe because the ionization-equilibrium cell loop is
  serial -- and (ii) assembles `fvec` by **looping over elements** (two balance
  equations each, force-zeroing absent elements). `params` now carries only H/He
  (1-11), HeITR (12-18), and the C/N/O charge-transfer rates (40-48); the old
  per-element slots 49-55 are gone. The driver builds metal densities, initial
  guesses, and the solution unpack by looping over `n_melem`/`n_mion`. The system
  size is now the formula `N_eq = 3 + 2*n_melem` (= 11 for the current four
  elements), set in `input_read.f90`.
- **F -- Generalized grid floor.** `set_energy_vectors.f90` no longer hardwires
  the sub-13.6 eV floor to Mg's 7.646 eV; it computes
  `e_sub_low = min` over active elements of the neutral threshold
  `mion_ethr(melem_i0(e))` that lies below `e_th_HI`. The trigger flag
  `thereis_lowIP_metal` is set generically in `input_read.f90` (any active element
  with a sub-13.6 eV neutral threshold), driven by a new canonical-order abundance
  array `melem_ab` (declared in `parameters.f90`, populated in `input_read.f90`).
  For the current set only C I (11.26 eV) and Mg I (7.646 eV) lie below 13.6 eV
  and the minimum is Mg I, so the floor is unchanged at 7.646 eV when Mg is
  present. (Consequence: a future C/N/O-without-Mg run would now extend the grid
  to 11.26 eV, which it did not before -- and once Phase 1b adds Na/K/Ca, whose
  neutral thresholds 5.14/4.34/6.11 eV sit below Mg's, the floor drops further.)

### Gate: bit-for-bit reproduction of the C/N/O + Mg run

With C/N/O + Mg active and the same `input.inp`/`metals.inp` as §2, the
refactored code reproduces the `mg_validation/cno_mg` reference to round-off:
log10 Mdot = 12.76; `Hydro_ioniz.txt`, `Hydro_ioniz_adv.txt`, and
`Ion_species.txt` are byte-identical to the frozen reference.
`Ion_species_adv.txt` differs at exactly one number -- line 1 (innermost ghost
cell), column 7 (the HeITR fraction): the reference holds 4.0048e-5 left in
uninitialized stack memory while the deterministic build writes 0.0; every other
column of that line and every other line match. This is a pre-existing reference
artifact (a HeITR ghost-cell slot never used in the H/He + metals run), not an
effect of the refactor.

### What is deliberately *not* generalized yet

- **Charge exchange with H** is still hardcoded for C/N/O only (the `fvec(1)`
  C,N,O sum and the `cxlo`/`cxup` mapping via `iel_C/O/N` in `System_HeH_metals`).
  Generalizing it is Phase 1d.
- **Per-ion atomic data and per-element abundance remain manual entries** -- the
  `sigma_tab(:,1..8)` assignments and the `sigma_XX` cross-section functions in
  `set_energy_vectors.f90`, the abundance scalars `X_C/X_O/X_N/X_Mg`, the
  `metals.inp` keys, and the `melem_ab` population block. These are the
  atomic-data / abundance rows expected to grow per element in Phase 1b/1c; they
  are data, not interface plumbing.

### Touched files

- **New:** `init/species_table.f90` (Step A).
- `init/parameters.f90` -- `melem_ab(:)` abundance array; generalized
  `thereis_lowIP_metal` comment.
- `files_IO/input_read.f90` -- `N_eq = 3 + 2*n_melem`; `melem_ab` population;
  generalized low-IP-metal trigger.
- `files_IO/write_output.f90` -- 2D metal-column writer (B).
- `radiation/util_ion_eq.f90` -- `PH_heat_HHe` + column density 2D (C);
  `eval_cool` metal channels 2D (D).
- `nonlinear_system_solver/System_HeH_metals.f90` -- module-level per-cell metal block
  + element-loop residual (E).
- `radiation/ionization_equilibrium.f90` -- metadata-driven density build, params
  packing (H/He + charge transfer only), initial guesses, solution unpack (E).
- `init/set_energy_vectors.f90` -- min-over-active-metals grid floor (F).

### Status and relation to the Huang plan

Phase 1a is complete: the per-element plumbing of §2 is now table-driven. Adding
the Phase-1b metals should require only metadata rows plus atomic/abundance data,
with **no changes to `System_HeH_metals.f90` or `ionization_equilibrium.f90`**. Next:
**Phase 1b** -- Si, Ca, Na, K, S added together as a uniform-template batch.

Gate scratch run: `/tmp/ates_gate/cno_mg/` vs frozen reference
`mg_validation/cno_mg/output/`.

## 4. Phase 1b: correct C/N/O cross sections and add Si, Ca, Na, K, S

*Added 2026-06-05.*

Two changes are bundled in this update, both enabled by the Phase 1a (§3)
table-driven plumbing:

1. **Corrected the C/N/O photoionization cross sections** to the full
   Verner, Ferland, Korista & Yakovlev (1996, hereafter VFKY96) Eq. (1) form.
   The previous code used the simpler Verner & Yakovlev (1995) inner-shell form
   (no `y_0`/`y_1` offset/asymptote parameters) for C/N/O; that form is exact
   only when `y_0 = y_1 = 0`, which is not the case for the C/N/O outer shells.
   This is the dominant *physics* change here.
2. **Added the five uniform-template metals Si, Ca, Na, K, S** in one batch
   (Phase 1b proper). With Mg (§2) this brings the in-code element set to
   C, N, O, Mg, Si, Ca, Na, K, S -- the full Huang+2023 §2.4 set **except Fe**
   (Phase 1c). Species count `n_species = 30`; metal-ion stages `n_mion = 24`;
   elements `n_melem = 9`.

The metal MINPACK solver was also **renamed** (see below).

### Variable-stage solver (2- and 3-stage elements in one packing)

Huang+2023 carries Si, Ca (like Mg) to second ionization but Na, K, S only to
first. The solver keeps a uniform **two-unknowns-per-element** packing
(`N_eq = 3 + 2*n_melem = 21` with all nine elements in the table) and a
per-element `melem_top` flag:

- `melem_top = 2` (C, O, N, Mg, Si, Ca): solve both X0<->X+ and X+<->X++.
- `melem_top = 1` (Na, K, S): solve only X0<->X+ and **pin the unused upper
  unknown to zero** in the residual (`fvec(ix+1) = x(ix+1)`), which keeps the
  Jacobian non-singular. Every `ix+1` array access is guarded by `melem_top` so
  a two-stage element never reads or writes its neighbor's column.

Elements with zero abundance are still force-zeroed, so the same 21-unknown
system solves C-only, C/N/O+Mg, or all nine without renumbering.

### Atomic data and rates (sources, not fabricated)

- **Photoionization cross sections** (`functions/cross_sec.f90`):
  - Full VFKY96 Eq. (1) form (`sigma_VFKY96`) for C I/II, N I/II, O I/II,
    **Si I/II, Ca I/II, S I** -- Verner+1996 Table 1 `[E_0,sigma_0,y_a,P,y_w,y_0,y_1]`.
  - Simpler Verner & Yakovlev 1995 form (`sigma_Verner96`, `y_0=y_1=0`) for
    Mg I/II, **Na I**, and **K I** (their published fits have `y_0=y_1=0`; the
    two forms then agree exactly).
  - Neutral thresholds [eV]: Si I 8.152, Ca I 6.113, Na I 5.139, K I 4.341,
    S I 10.36.
- **Collisional ionization:** Voronov (1997) fits for all new ions.
- **Recombination:** Badnell RR+DR via the swappable dispatcher, as for C/N/O
  and Mg. **This is the same deliberate deviation** from the `Huang_update_plan.md`
  §1.3 validation rates (Shull & van Steenberg 1982 for Mg/Si/S/Ca; Verner &
  Ferland 1996 for Na/K) already accepted for Mg -- the Phase-6 production rate
  applied early. Read the ion/neutral balances below as the Badnell result.
- **Line cooling:** **not ported** for any of the five (`mion_iscool = .false.`),
  flagged To Be Checked. Like Mg, they act only as electron donors and
  photoelectric-heating channels, not coolants (Phase 2).
- **Charge exchange with H:** **not added** for the new elements; still C/N/O
  only (Phase 1d).

### Rename: `System_HeHCO` -> `System_HeH_metals`

The metal MINPACK module/routine/setter were named after the AIOLOS heritage
(AIOLOS tracked only C and O as metals -> "HCO" = H, C, O). The solver now
carries H + He + nine metals, so the C/O-specific name was misleading. Renamed
to parallel the existing `System_HeH` / `System_HeH_TR` family:

| old | new |
|---|---|
| `module System_HeHCO` / `System_HeHCO.f90` | `module System_HeH_metals` / `System_HeH_metals.f90` |
| `ion_system_HeHCO` (residual) | `ion_system_HeH_metals` |
| `set_hco_metals` (per-cell setter) | `set_metal_coeffs` |
| `hco_ntot/g0/g1/b0/b1/a1/a2/top` | `met_ntot/g0/g1/b0/b1/a1/a2/top` |

Pure rename, no behavior change. (`§3` above has been updated to the new names.)

### Touched files

Data/metadata rows + the variable-stage guard; no hydro or `calc_ne`/`calc_rho`
changes (metal electrons still enter only inside the MINPACK system):

- `init/species_table.f90` -- `n_mion=24`, `n_melem=9`, `n_mphot=15`; rows for
  Si/Ca/Na/K/S; `melem_top = [2,2,2,2,2,2,1,1,1]`.
- `functions/cross_sec.f90` -- C/N/O switched to `sigma_VFKY96`; new
  `sigma_SiI/SiII/CaI/CaII/NaI/KI/SI/SII`.
- `radiation/Cool_coeff.f90` -- Badnell `rec_fit`/`lookup_rec` + Voronov for the
  new ions.
- `nonlinear_system_solver/System_HeH_metals.f90` -- `met_top` array; variable-
  stage residual (pin unused upper unknown); rename.
- `radiation/ionization_equilibrium.f90` -- `melem_top` guards on density build,
  initial guess, and solution unpack; `set_metal_coeffs` call; rename.
- `init/parameters.f90`, `files_IO/input_read.f90`,
  `files_IO/metals_input_read.f90` -- `n_species=30`; abundances `X_Si/Ca/Na/K/S`;
  `metals.inp` keys; `melem_ab` rows.
- `init/set_IC.f90`, `files_IO/load_IC.f90` -- neutral-metal IC for the new
  `f_sp` columns 19-30.
- `files_IO/write_output.f90` -- already 2D (§3); now writes 24 metal columns.
- `Makefile` -- `System_HeHCO.f90` -> `System_HeH_metals.f90` in `SRC`.

### Validation on WASP-121b

Three converged runs, identical planet / SED / grid (power-law SED, He/H=0.0851,
Roche domain, `Rate/2 + Mdot/2`); they differ only in code version and metal set.

| key | code | metals | C/N/O sigma | grid floor [eV] | log10 Mdot | T_peak [K] | r(T_peak) [R_p] | steps |
|---|---|---|---|---:|---:|---:|---:|---:|
| `ref` | old (pre-1b) | C,N,O,Mg | Verner96 (VY95) | 7.646 (Mg I) | 12.76 | 9828 | 1.224 | 5003 |
| `A`   | new (1b)     | C,N,O,Mg | VFKY96          | 7.646 (Mg I) | 12.85 | 12197 | 1.899 | 15382 |
| `B`   | new (1b)     | + Si,Ca,Na,K,S | VFKY96     | 4.341 (K I)  | 12.85 | 12135 | 1.892 | 15505 |

Mdot in log10 g/s. Huang+2023 Case A (spherical) reference: log10 Mdot = 12.57.
`ref` and `A` share the 7.646 eV grid floor (Mg sets the minimum in the C/N/O+Mg
set), so `ref` vs `A` isolates the C/N/O cross-section change with no grid-floor
confound. `A` vs `B` changes both the element set and the grid floor (K I pulls
the floor to 4.341 eV).

Observations (for evaluation, not pass/fail):

- **C/N/O cross-section correction (`ref` -> `A`):** dlog10(Mdot) = **+0.090 dex**;
  T_peak rises from 9828 K to **12197 K** and the peak moves outward from 1.22 to
  ~1.90 R_p. The corrected VFKY96 C/N/O opacity deposits more photoheating, and
  the hotter peak now sits near Huang's quoted ~12000 K. This is a genuine model
  change, not numerical.
- **Adding Si/Ca/Na/K/S (`A` -> `B`):** dlog10(Mdot) = **+0.000 dex**, dT_peak
  = -62 K, despite the grid floor dropping to 4.341 eV. At solar abundance and
  with no line cooling these five are radiatively near-passive electron donors.
- **Magnesium moves indirectly.** Mg's own cross sections are byte-identical
  between `ref` and `A`, yet its ionization fractions shift by up to ~0.1: an
  *indirect* response to the hotter, more extended profile (changed n_e, optical
  depth, photoionizing flux), not a change in Mg atomic data. Mg is coupled to
  the global state, so it is **not** an isolated control here.
- **Numerics.** Zero NaN/Inf in both new runs; the 21-unknown system converges
  across the full T range; per-element X/H is conserved to ~1e-13 for **all nine
  elements including the two-stage Na/K/S** -- direct evidence the variable-stage
  packing does not leak between element columns.
- **New-metal ionization structure (run B).** Si and Ca are mostly **doubly**
  ionized in the upper thermosphere; Na, K, S are mostly **singly** ionized
  (no second stage carried) -- the same qualitative ordering as Huang+2023 Fig. 12.

### Caveats and gaps

- **No metal line cooling** for Si/Ca/Na/K/S (`mion_iscool=.false.`) -- Phase 2.
  The thermal structure will not match Huang until Fe II / Mg II / Ca II cooling
  is in.
- **No charge exchange with H** for the new elements (or Mg) -- Phase 1d.
- **Recombination is Badnell, not the §1.3 validation rate** (deliberate, as Mg).
- **Mdot is ~1.9x Huang Case A** (12.85 vs 12.57); the metal-set (Fe still
  missing) and geometry confounds from Phase 0 still apply.
- **Not done here:** Huang Fig. 12 is not digitized/overlaid (only qualitative
  ordering checked); the atomic-data values themselves are not re-derived beyond
  conservation and smoothness.

### Regression reference for Phase 1c

The accepted run `A` output is frozen at
`phase1b_validation/.gate_ref_cno_mg_vfky96/output/` as the post-Phase-1b
C/N/O+Mg baseline. When Phase 1c adds Fe, running with Fe **off** must reproduce
this baseline. Note that adding Fe bumps `n_species`/`n_mion` and therefore the
`Ion_species.txt` column count, so the Phase-1c check is a **round-off match on
the shared physical columns**, not a byte-for-byte file compare (the old
`mg_validation/cno_mg` byte-gate applied only to the physics-neutral Phase 1a
refactor, where `n_species` did not change).

### Status and relation to the Huang plan

This completes **Phase 1b**: the ionization network now carries C/N/O/Mg/Si/Ca/Na/K/S,
the C/N/O cross sections are corrected to VFKY96, and the variable-stage solver
handles 1st- and 2nd-ionization elements uniformly. The Phase-1b gate's
"converges without NaNs, conserves X/H, Fig. 12 ionization ordering, Mdot within
~2x Case A" items are addressed by the numbers above; the thermal-balance part
remains deferred to Phase 2 (metal line cooling). Next: **Phase 1c** -- iron
(Fe I/II/III), which needs non-uniform cross-section (Zatsarinny+2019) and
recombination data.

Scratch runs: `phase1b_validation/{.gate_ref_cno_mg (old), cno_mg_new (A), all_metals (B)}/`;
comparison notebook `phase1b_validation/Phase1b_validation.ipynb`.

## 5. Phase 1c: add iron (Fe I / II / III)

*Added 2026-06-05.*

Iron is added as the third three-stage trace metal, completing the Huang+2023
§2.4 element set (C, N, O, Mg, Si, Ca, Na, K, S, **Fe**). No solver, cooling-loop,
heating, output, or post-processing code changed: the variable-stage MINPACK
packing and every metal loop are generic over `n_melem`, so adding iron is purely
a matter of new metadata rows and four new rate functions. Counts move to
`n_species = 33`, `n_mion = 27`, `n_melem = 10`, `n_mphot = 17`; the MINPACK
system grows to `N_eq = 3 + 2*n_melem = 23`. Iron carries `melem_top = 2` (three
stages, like C/O/N/Mg/Si/Ca).

### Atomic data and rates (sources, not fabricated)

Iron is the one element whose data are **not** the uniform-template choice. Each
quantity below was checked against its primary source before encoding; in
particular the Voronov fit was cross-checked because an earlier abandoned attempt
had its `P`/`X`/`K` columns transposed.

- **Photoionization cross sections** (`functions/cross_sec.f90`): single-shell
  VFKY96 form (`sigma_VFKY96`) for Fe I and Fe II, parameters from the
  Verner+1996 `phfit2.f` `PH2` block. Neutral thresholds [eV]: Fe I **7.902**
  (NIST ionization potential), Fe II **16.199**.
  - **To Be Checked / known limitation:** Huang+2023 use the multi-shell
    Zatsarinny+2019 Fe I cross section. The single-shell VFKY96 fit here
    under-represents inner-shell absorption at high XUV energy, so the iron
    photoionization rate is a lower bound at the hardest photons. Flagged inline
    in `cross_sec.f90`. A multi-shell upgrade is left for later refinement.
- **Collisional ionization:** Voronov (1997) `cfit.dat` fits.
  Fe I (`dE=7.9, A=2.52e-7, P=0, X=0.701, K=0.25`),
  Fe II (`dE=16.2, A=2.21e-8, P=1, X=0.033, K=0.45`).
- **Recombination:** the Huang+2023 Eqs. (5)-(6) total fits (a dielectronic
  Arrhenius term `A T^{-1.5} e^{-T0/T}(1 + B e^{-T1/T})` plus a radiative
  power-law `C (T/10^4)^{-eta}`) used **directly**, not the Badnell dispatcher
  used for the lighter metals. This is a deliberate choice: the Huang iron rates
  are the validation target for this phase, so the model is run on exactly those
  rates. The Badnell-swap for iron is deferred to Phase 6.
  - Fe II $+\,e\rightarrow$ Fe I: `A=2.833e-8, T0=5.731e4, B=1.383e4, T1=120.4;  C=1.248e-12, eta=0.485`.
  - Fe III $+\,e\rightarrow$ Fe II: `A=1.094e-5, T0=1.490e4, B=36.74, T1=1.153e5;  C=1.728e-12, eta=0.618`.
- **Line cooling:** **not** ported (`mion_iscool = .false.` for all Fe stages) --
  Phase 2. This matters more for iron than for the lighter metals: Fe II is the
  **dominant low-altitude line coolant** in Huang+2023, so a metals-on thermal
  profile is **expected to run hotter than Huang** until Phase 2.
- **Charge exchange with H:** **not** added; still C/N/O only (Phase 1d). Iron's
  Fe$^+$+H$^+\!\leftrightarrow$Fe$^{2+}$+H reactions (Huang Table 4) come in Phase 1d.

### Touched files

Metadata rows + four rate functions + abundance/IC wiring; no solver, cooling,
heating, output, or post-process changes:

- `init/species_table.f90` -- `n_mion=27`, `n_melem=10`, `n_mphot=17`; Fe I/II/III
  rows (`f_sp` cols 31-33); `iel_Fe=10`; `melem_top` appends `2`.
- `functions/cross_sec.f90` -- `sigma_FeI`, `sigma_FeII` (single-shell VFKY96).
- `init/set_energy_vectors.f90` -- `sigma_tab(:,16)=sigma_FeI`, `(:,17)=sigma_FeII`.
- `radiation/Cool_coeff.f90` -- Huang recombination `alpha_rec_FeI/FeII_Huang`
  with `rec_FeII`/`rec_FeIII` array wrappers; Voronov `ion_coeff_FeI/FeII`;
  dispatcher cases for `FeII`/`FeIII` recombination and `FeI`/`FeII` ionization.
  No `cool_coeff_metal` case for Fe (Phase 2).
- `init/parameters.f90`, `files_IO/input_read.f90`,
  `files_IO/metals_input_read.f90` -- `n_species=33`; abundance `X_Fe`;
  `metals.inp` keys `FeI`/`Fe`; `melem_ab(iel_Fe)` row; `thereis_metals` trigger.
- `init/set_IC.f90`, `files_IO/load_IC.f90` -- neutral-iron IC for `f_sp` cols 31-33.

### Iron abundance

Asplund et al. (2009) solar Fe/H = 3.16e-5 (log eps = 7.50). The Huang+2023
WASP-121b model adopts 4x solar metallicity, so the validation run uses
`X_Fe = 1.26e-4` by number relative to H. C/N/O/Mg are held at solar so the run
isolates the iron change.

### Validation on WASP-121b

Two converged runs (identical planet / SED / grid; power-law SED, He/H=0.0851,
Roche domain, "Rate/2 + Mdot/2") plus the frozen Phase-1b C/N/O+Mg baseline as
the regression reference:

| key | metals | iron | log10 Mdot | T_peak [K] | r(T_peak) [R_p] | steps |
|---|---|---|---:|---:|---:|---:|
| `ref` | C,N,O,Mg | -- (frozen Phase-1b run A) | 12.85 | 12197 | 1.899 | 15382 |
| `off` | C,N,O,Mg | off (`X_Fe=0`) | 12.85 | 12197 | 1.899 | 15382 |
| `on`  | C,N,O,Mg + **Fe** | on, 4x solar | 12.87 | 12223 | 1.915 | 16566 |

Mdot in log10 g/s; Huang+2023 Case A reference log10 Mdot = 12.57.

Observations (for evaluation, not pass/fail):

- **Regression, force-zeroed iron (`off` vs `ref`).** The Phase-1c binary carries
  three extra species and two extra MINPACK unknowns, yet run `off` reproduces the
  frozen baseline on every shared column to the round-off floor: worst relative
  difference **2.0e-13** across the seven hydro columns and **3.7e-13** (N I)
  across the 30 shared species. The three new Fe columns are **identically zero,
  no NaNs**. When `X_Fe=0` the iron block reduces to identity rows in the residual
  (`fvec=x`, rooted at zero) and contributes exactly zero electrons, so the
  23-unknown system returns the 21-unknown answer -- adding iron is
  non-perturbative when off.
- **Iron ionization structure (`on`).** Fe II dominates the base (fraction 0.89 at
  $r=1$, peaking ~0.96 near $r=1.02$); Fe III dominates the escaping flow (0.97 at
  the outer edge), with the Fe II/Fe III crossover near $r\simeq1.25$; Fe I is
  minor everywhere (max ~0.11 at the base, a small bump ~0.08 near $r=1.2$). This
  is the qualitative arrangement of Huang+2023 Fig. 12 -- Fe I confined to the
  base, Fe II at the base, Fe III through the outflow.
- **Iron as an electron donor.** At 4x solar, iron supplies up to **~26%** of the
  free-electron density, so it feeds back on the whole coupled solution, not just
  its own columns.
- **Switching iron on (`off` -> `on`).** $\Delta\log_{10}\dot M = +0.02$ dex,
  $\Delta T_{\rm peak} = +26$ K (12197 -> 12223). The change is small and in the
  *heating* direction: with no Fe line cooling yet, iron acts only as an added
  electron donor / photoelectric-heating channel, so the thermosphere runs
  slightly hotter rather than cooler.
- **Numerics.** Run `on` converged with zero NaNs; per-element X/H is conserved to
  **8.3e-13** across all five active elements including iron.

### Caveats and gaps

- **No iron line cooling** (`mion_iscool=.false.`) -- Phase 2. Because Fe II is the
  dominant low-altitude coolant in Huang+2023, the run-`on` **temperature profile
  is not yet comparable to Huang** and is expected to sit above it. The physics
  quantity validated at Phase 1c is the iron *ionization structure*, not $T(r)$.
- **No iron charge exchange with H** -- Phase 1d.
- **Single-shell VFKY96 Fe I/II cross sections** under-represent inner-shell XUV
  vs Huang's Zatsarinny+2019 Fe I (flagged in `cross_sec.f90`).
- **Recombination is the Huang fit**, not Badnell (deliberate; Badnell-swap is
  Phase 6).
- **Mdot is still ~2x Huang Case A** (12.87 vs 12.57); the missing metal cooling
  and the Phase-0 geometry confounds still apply.
- **Huang Fig. 12 is not digitized/overlaid** -- only the qualitative Fe II-base /
  Fe III-outflow ordering is checked.

### Regression reference for Phase 1d / Phase 2

The run-`on` output (C/N/O/Mg + 4x-solar Fe) is frozen at
`phase1c_validation/.gate_ref_cno_mg_fe/` (input, `metals.inp`, `run.log`, and
`output/`) as the post-Phase-1c iron-on baseline. Phase 1d (iron charge exchange)
and Phase 2 (iron line cooling) will both *change* the iron ionization and the
thermal structure deliberately, so this baseline is the "before" snapshot for
those comparisons rather than a bit-reproduction target.

### Status and relation to the Huang plan

Phase 1c **implementation is complete** and the iron ionization structure is
presented for validation against Huang+2023 Fig. 12 (numbers and figures in the
notebook; no pass/fail asserted here -- that judgement is the reader's). The
ionization-network half of Phase 1 now carries the full Huang element set with
verified iron atomic data. The thermal-balance half remains deferred: until Fe II
line cooling lands in Phase 2 the metals-on temperature profile will not match
Huang. Next: **Phase 1d** -- charge exchange with H for all metals (Huang Table 4),
including Fe$^+$+H$^+\!\leftrightarrow$Fe$^{2+}$+H.

Scratch runs: `phase1c_validation/{regression_cno_mg (off), cno_mg_fe (on),
.gate_ref_cno_mg_fe (frozen)}/`; comparison notebook
`phase1c_validation/Phase1c_validation.ipynb`.

## 6. Phase 1d: charge exchange with hydrogen (Huang Table 4)

*Added 2026-06-05.*

Phase 1d adds the charge-exchange (charge-transfer, CT) reactions of Huang+2023
Table 4 as a single table-driven module. The default active set is the 23
**Group A** metal + H / H$^+$ reactions, which couple each metal's ionization
balance to the H$^+$ fraction. A runtime switch `cx_full 1` in `metals.inp`
enables the rest of the table -- **Group B** (He + H, 2 rows), **Group C**
(metal + He, 6 rows), **Group D** (metal + metal, 32 rows). This **replaces** the
previous hard-coded Kingdon & Ferland (1996) C/N/O + H charge-exchange terms in
the solver with a generic implementation covering all ten metals plus He.

Re-sourcing C/N/O to Huang Table 4 has one deliberate consequence: Table 4
tabulates only the neutral$\leftrightarrow$singly-ionized CT for C, N, O
(rows A13--A18), so the switch **drops** the doubly-ionized
C$^{2+}$/N$^{2+}$/O$^{2+}$ + H recombination terms the old KF96 implementation
carried. This is an intentional behavior change to match the paper; it raises the
C/N/O double-ionization fractions and reshapes the thermal profile (documented
below). No solver structure, cooling loop, heating, output, or post-processing
code changed beyond the CT wiring.

### Atomic data and rates (sources, not fabricated)

All 63 rate rows are transcribed verbatim from `docs/charge_exchange_table4.md`,
the image-verified transcription of Huang+2023 Table 4 (23 Group A + 2 Group B +
6 Group C + 32 Group D). The paper text says "65 reactions"; the 2-row
discrepancy is immaterial and noted in the reference doc.

- **Functional forms** (`radiation/charge_exchange.f90`): KF96 bracket
  `a*T4^b*(1+c*exp(d*T4))*exp(-Ek/T4)` (`kf`); power law `a*(T/300)^b*exp(-Ek/T4)`
  (`gp`); 4th-order $\ln T$ polynomial (`p4`) for the radiative-CT rows of
  N / S / Na / K (Lin/Zhao/Dutta/Watanabe); and the Satta (2013) $f_{\rm Si}(T)$
  and $D_{\rm CSi}(T)$ helpers for the Si + He and C + Si rows.
- **Spot checks at $T=10^4$ K** (notebook step 3, against the table):
  Fe + H$^+$ $= 4.000\times10^{-9}$ (RV72), Fe$^+$ + H $= 2.011\times10^{-12}$,
  Fe$^{2+}$ + H $= 1.260\times10^{-9}$ (ND87), O + H$^+$ $= 2.091\times10^{-9}$
  (near-resonant, S99) -- all match.
- **Rate clamp** $[0,\,10^{-6}]$ cm$^3$ s$^{-1}$: the upper bound guards the $\ln T$
  polynomials, which diverge as $T\to 1$ K (never reached); the lower bound zeroes
  a few KF96 bracket fits that dip slightly negative below their validity floor
  (e.g. Mg$^+$ + H$^+$ at $T\lesssim1500$ K), where the true rate is negligible.
- **Reference-doc notes:** Ca has **no** direct H charge exchange (only Group D);
  K has only the forward K + H$^+$ row. Reverse rates for pairs lacking direct
  data are Huang's detailed-balance estimates, carried verbatim.

### Touched files

A new module plus minimal wiring; the solver's old per-element C/N/O CT block is
deleted, not extended:

- `radiation/charge_exchange.f90` **(new)** -- the full Table-4 descriptor arrays,
  the `cx_rate` dispatcher, and the generic residual assembly. Public entry points:
  `cx_init` (build the active reaction list from `cx_full`), `cx_set_cell(T)`
  (evaluate per-cell rates), `cx_add_to_fvec` (add CT source terms to the MINPACK
  residual via the `cx_fvidx` element/stage $\to$ row map), and the `cx_full`
  switch. The zero-abundance guard is automatic: each rate is
  $k\,n_{\rm donor}\,n_{\rm acceptor}$, which vanishes when either reactant is absent.
- `nonlinear_system_solver/System_HeH_metals.f90` -- the hard-coded C/N/O + H
  charge-exchange terms are **removed**; a single `call cx_add_to_fvec(...)` now
  feeds CT contributions into the H, He, and every metal-stage residual row.
- `radiation/ionization_equilibrium.f90` -- `call cx_set_cell(T_K(j))` before each
  cell's `hybrd1` solve, storing the per-cell rates used inside the residual.
- `files_IO/metals_input_read.f90` -- parses `cx_full <0|1>` from `metals.inp`.
- `files_IO/input_read.f90` -- `if (thereis_metals) call cx_init` after the input
  read, before the ionization sweep.

### Validation on WASP-121b

Two new converged runs (identical planet / SED / grid as Phase 1c) compared
against the frozen Phase-1c iron-on output as the regression reference:

| key | charge exchange | log10 Mdot | T_peak [K] | r(T_peak) [R_p] | steps |
|---|---|---:|---:|---:|---:|
| `ref`  | KF96 C/N/O + H (frozen Phase-1c run `on`) | 12.87 | 12223 | 1.915 | 16566 |
| `cx`   | Huang Table 4 Group A (all metals + H), default | 12.89 | 12955 | 1.778 | 13942 |
| `full` | Huang Table 4 full (A+B+C+D), `cx_full 1` | 12.88 | 12951 | 1.765 | 11661 |

Mdot in log10 g/s; Huang+2023 Case A reference log10 Mdot = 12.57. All three runs
converged with **zero NaNs**.

Observations (for evaluation, not pass/fail):

- **H$^+$ sink magnitude.** In the ionized bulk ($x_{\rm HII}>0.5$) the net metal
  charge exchange is at most $3.7\times10^{-3}$ of radiative recombination (median
  $1.0\times10^{-3}$) -- negligible as an H$^+$ sink, consistent with Huang
  \S4.2 ("Fe--H CX lower than H photoionization and recombination"). At the
  near-neutral base ($x_{\rm HII}\sim7.8\times10^{-8}$) the *ratio* is large
  ($\sim$118) only because H$^+$ is itself a trace species there.
- **Near-resonant balance.** O + H$^+$ $\leftrightarrow$ O$^+$ + H is near-resonant
  and sits in near-equilibrium at the base (gross destroy $3.66\times10^5$ vs
  produce $3.69\times10^5$, net $+3.0\times10^3$). The net base H$^+$ sink is
  instead driven by the fast, one-sided Fe + H$^+\!\rightarrow$Fe$^+$ + H channel
  (Fe gross destroy $3.31\times10^4$, produce $4.06\times10^3$, net
  $-2.9\times10^4$).
- **C/N/O double ionization rises** (the dropped 2+ CT). Column-max fractions
  `ref`$\to$`cx`: O\,III $0.871\to0.946$, N\,III $0.897\to0.944$,
  C\,III $0.851\to0.856$. Removing the C$^{2+}$/N$^{2+}$/O$^{2+}$ + H recombination
  channels lets the double-ions build up.
- **Iron base shift.** Fe\,II base fraction $0.889\to0.896$, Fe\,III essentially
  vanishes at the base ($0.002\to0.000$): Fe + H$^+\!\rightarrow$Fe$^+$ + H pushes
  iron toward Fe$^+$ at the base.
- **Thermal structure.** $T_{\rm peak}$ **rises** $12223\to12955$ K and moves
  inward ($1.915\to1.778\,R_p$). This is the **combined** effect of (i) the C/N/O
  re-sourcing that drops the 2+ recombination terms and (ii) the new metal + H CT.
  Charge exchange adds **no radiative-loss channel**, so it cannot cool the gas;
  Fe\,II / metal line cooling remains deferred to Phase 2, and the metals-on
  profile therefore still runs *hotter* than Huang, now more so.
- **Group A vs full Table 4 (`cx` vs `full`).** Turning on Groups B/C/D changes the
  bulk temperature by at most 1.5% (worst relative difference), with
  $\max|\Delta T| = 556$ K at $r=1.008$; at $r=1.02$ the full run is 357 K cooler
  with more base H$^+$ and less He\,II. $\dot M$ moves $-0.01$ dex and the outflow
  is nearly unchanged -- the extra He--H / metal--He / metal--metal rows mostly
  redistribute the base He\,II / H$^+$ balance.
- **Dynamics.** $v_{\rm out}$ and sonic radius: `ref` 18.83 km/s @ $1.884\,R_p$;
  `cx` 20.28 km/s @ $1.826\,R_p$; `full` 20.40 km/s @ $1.826\,R_p$. The hotter,
  more inward thermal peak drives a slightly faster wind with the sonic point
  $\sim0.06\,R_p$ further in.

### Caveats and gaps

- **No metal line cooling yet** (Phase 2). Charge exchange adds no loss term, so
  $T(r)$ is not comparable to Huang; the C/N/O re-sourcing makes the base hotter,
  not cooler.
- **C/N/O 2+ + H terms intentionally dropped** to match Table 4 -- a real change
  from the Phase-1c (KF96) reference, and the dominant driver of the O\,III/N\,III
  rise rather than the new metal CT alone.
- **Mdot still $\sim$2x Huang Case A** (12.88--12.89 vs 12.57); the missing metal
  cooling and the Phase-0 geometry confounds still apply.
- **Paper lists 65 reactions, 63 transcribed** -- immaterial, documented in
  `charge_exchange_table4.md`.
- **No Huang figure digitized/overlaid** for the H$^+$ budget; only the
  relative-magnitude statement (net CT $\ll$ recombination in the bulk) is checked.

### Regression reference for Phase 2

The default Group-A run is frozen at
`phase1d_validation/.gate_ref_cno_mg_fe_cx/` (input, `metals.inp`, `run.log`,
`output/`) as the post-Phase-1d baseline. Because the C/N/O re-sourcing
deliberately shifted both the ionization and the thermal structure away from the
Phase-1c frozen reference, **this** snapshot -- not the Phase-1c one -- is the
"before" reference for Phase 2 (metal line cooling, including Fe\,II).

### Status and relation to the Huang plan

Phase 1d **implementation is complete**: the 23 Group-A metal + H reactions are
active by default and the full Table 4 sits behind `cx_full`. With this, the
**ionization-network half of Phase 1 is physically complete** relative to
Huang+2023 \S2.4 + Table 4 -- the full element set, verified atomic data, and the
metal--H charge-exchange coupling are all in place. The thermal-balance half
remains deferred: until metal line cooling (Fe\,II in particular) lands in
**Phase 2**, the metals-on temperature profile stays above Huang. Next: **Phase 2**
-- metal line cooling.

Scratch runs: `phase1d_validation/{cno_mg_fe_cx (Group A),
cno_mg_fe_full (full Table 4)}/`; comparison notebook
`phase1d_validation/Phase1d_validation.ipynb`.

## 7. Phase 2: metal line cooling — validation gate closed

Phase 2 brings the CHIANTI-based metal line cooling (Mg I/II, Ca II, Na I, Fe II,
on top of the C/N/O coolants) into the coupled energy balance and **closes the
Phase 2 validation gate** against Huang et al. (2023) Fig. 10. The per-coolant
atomic data and the CHIANTI effective-cooling tables are documented in
`Update_EXHALE_early_phase` (Parts II–III); this section records
the integrated WASP-121b gate and the diagnostic used to close it.

### Active line coolants

The line-cooling ions are flagged by `mion_iscool` in `species_table.f90`
(canonical order): C I, C II, O I, O II, N I, N II, Mg I, Mg II, Ca II, Na I,
Fe II. Each contributes `beta_esc * ne * n_ion * Lambda_ion(T)` with the
effective cooling coefficient `Lambda_ion` interpolated from the hardcoded
`cool_logL_*` tables in `Cool_coeff.f90` (NCOOLT = 41, log10 T = 3.0–5.0,
dlogT = 0.05). In the optically thin limit `beta_esc = 1` (see
`Update_EXHALE_early_phase`, Part II, "Adopted assumption: 100% escape"). N I/II carry the
`iscool` flag but have no tabulated line cooling (`Lambda = 0`), so they
contribute exactly zero.

### Per-channel cooling diagnostic (this phase)

To close the gate quantitatively, `eval_cool` (`util_ion_eq.f90`) gained an
optional `cool_chan` output — an **exact** decomposition of the total cooling
rate into H/He recombination, collisional ionization, collisional excitation,
bremsstrahlung, and one column per metal ion. It is read straight from the arrays
`eval_cool` already builds, so the channel sum reproduces the `cool` total to
machine precision (no offline re-derivation, no divergence risk). A companion
writer `write_cool_breakdown_eq` (same module, called once from `EXHALE_main.f90`
after the final equilibrium write) dumps `output/Cooling_breakdown.txt`: columns
`r/Rp, T, ne, cool_total`, the four H/He channels, then 27 metal-ion columns. The
run prints `max |sum(channels)/cool - 1|` as an internal consistency check.

### Validation on WASP-121b (the gate)

Converged equilibrium run (full metals, the curated `metals.inp`; same SED/grid
as Phase 1d):

- Internal consistency: `max |sum(channels)/cool - 1| = 3.7e-16` (machine eps) —
  the breakdown is a bit-exact decomposition of the `Hydro_ioniz.txt` cool column.
- In the **1.15 < r/Rp < 1.4** band the dominant radiative coolants are
  **Fe II (~30%)** and **Mg II (~25%)**, jointly ~55%; metals are **57–85%** of
  all radiative cooling through the band. Fe II leads the inner edge (top coolant
  at 1.1–1.2 R_p, 53% at 1.1 R_p), Mg II takes over by 1.3–1.4 R_p, and H
  recombination only overtakes above ~1.5 R_p — the Fe II → Mg II handoff Huang
  report in their Fig. 10.
- Mdot unchanged (log10 Mdot = 12.51 g/s); zero NaNs.

This matches the Phase 2 gate ("a metal-driven cooling feature appears at
1.15–1.4 R_p, Fig. 10") and `Huang_update_expected_results.md` (Mg II + Fe II
become the most important coolant in 1.15 < r/Rp < 1.4).

### Touched files

- `src/modules/radiation/util_ion_eq.f90` — optional `cool_chan` out-arg on
  `eval_cool` (exact per-channel split); new `write_cool_breakdown_eq` writer;
  `mion_fsp` added to the `species_table` use-list.
- `EXHALE_main.f90` — `use utils_ion_eq, only: write_cool_breakdown_eq`; call it
  after the final equilibrium `write_output`.
- Cooling coefficients themselves live in `Cool_coeff.f90` (provenance in
  `Update_EXHALE_early_phase`, Parts II–III); unchanged here.

### Caveats and deferred

- The per-ion metal decomposition is exact in the default branch
  (`use_2lev_cool = .false.`, the WASP-121b setting). In the optional two-level
  fine-structure branch the per-ion terms are the resonance-line approximation
  and need not sum to the two-level `cool_M`.
- **Fe II** uses Boltzmann/coronal lower-level populations (E_cut = 38459 cm⁻¹),
  validated vs Huang Fig. 6 coronal curve within ~1.5× — sufficient for the gate.
  The Hollenbach & McKee IR fine-structure + density-dependent two-level
  refinement (Huang "Sum" curve) is deferred.
- **Fe I** line cooling (NIST f-values + Van Regemorter, Huang Fig. 5) is not yet
  added; Fe I is a minor coolant in the gate band. *(Now added in §9; confirmed
  negligible — peak 3e-5 of total cooling.)*
- The **CHIANTI Lyα swap** is intentionally skipped in Phase 2; the 0.35× Lyα is
  a Phase-3 radiative-transfer suppression, not a cooling-rate change.
- Carrying metal cooling into the post-processed `_adv` profile is handled
  separately (runtime `pp_metals`, default frozen); see the post-process notes.

### Status and relation to the Huang plan

**Phase 2 is complete and the gate is closed** (user-validated 2026-06-05). The
metals-on thermal balance now reproduces the Huang Fig. 10 Fe II/Mg II cooling
feature, closing the thermal-balance half of the network that Phase 1 left open.
Next: **Phase 3** — in-code excited H(n=2) + Lyα radiative transfer.

## 8. Phase 3a: in-code excited H(n=2) + Balmer proton source and heating

*Added 2026-06-05.*

Phase 3a is the first half of Phase 3: an **in-code H(n=2) population model** whose
output feeds back into both the hydrogen ionization balance (an extra proton
source) and the energy equation (Balmer photoelectric heating), fully coupled.
The second half (a Lyα radiative-transfer solve to supply a self-consistent mean
intensity J̄_Lyα(r)) is deferred; this phase uses a **parameterized** J̄_Lyα with a
switch to read an external RT profile instead. Everything is **default-off**: with
no stellar parameters in `input.inp` the build is byte-identical to the Phase-2
result (`use_excited_H = .false.`).

### Physics and atomic data (sources, not fabricated)

- **H(2s)/H(2p) populations:** a 2×2 statistical-equilibrium system following
  Christie, Arras & Li (2013) — recombination cascade into n=2, collisional
  2s↔2p mixing, two-photon 2s→1s decay, Lyα 2p→1s with a 2p radiative-pumping
  term, and 2s/2p photoionization. Solved per cell (`n2_populations`).
- **Lyα mean intensity** (the 2p pump): parameterized as
  J̄_Lyα ≈ 0.1 · F_LyC / Δν_D (Huang et al. 2017, Eq. 6) with a **top-down
  1/(1+τ_Lyα) line-center attenuation** so the optically thick base
  (τ_Lyα ~ 1.1e8 here) is not over-pumped. Oscillator strength f_Lyα = 0.4162,
  line cross-section constant C_Lyα = 1.49736e-2 · f_Lyα.
- **n=2 Balmer-continuum photoionization rate Γ₂:** from a **diluted stellar
  blackbody** above the n=2 threshold E₂ = 3.40 eV (`gamma_n2_balmer`), with
  σ_LyC = 6.3e-18 cm² and threshold cross-section s₂ = 1.4e-17 cm². Γ₂ carries the
  **same ξ dayside-dilution factor as ground-state EUV ionization** (Rate/2 → 0.5,
  Rate/4 → 0.25); for WASP-121b's "Rate/2 + Mdot/2" (ξ = 0.5) this gives
  Γ₂ ≈ 613 s⁻¹, comparable to Huang's ≈576 s⁻¹.
- **Heating:** Balmer photoelectric heating H_pe = n₂ · (photoelectron energy ×
  Γ₂) per cell (`heat_n2_balmer`); optional collisional de-excitation heating
  behind `incl_deexc_heat` (off by default — it overlaps the existing H I
  collisional-excitation cooling).

### Coupling into the solver

The model runs **decoupled (lagged one outer step)** to keep it out of the
per-cell Newton solve: `excited_H_update` fills two global per-cell arrays from
the previous converged state, and `ioniz_eq` injects them:

- `gph_balmer_HI` [s⁻¹] is added to the H I photoionization rate `P_HI` — an
  extra **proton source** (S_proton = Γ₂ · n₂).
- `heat_balmer` [erg cm⁻³ s⁻¹] is added to the photoheating rate.

`excited_H_update` is called once per outer iteration in `EXHALE_main.f90` before
`ioniz_eq`, and once more after convergence to refresh the diagnostics. The
post-loop refresh deliberately does **not** re-solve `ioniz_eq` (the Balmer
heating is ~1e-4 of the total, not worth desyncing T/f_sp from the converged
hydro state).

### Runtime interface

Append to `input.inp` (omit for the default Phase-2-identical behavior):

```
Stellar Teff [K]: 6459.0
Stellar radius [R_sun]: 1.458
```

These set `use_excited_H` (both must be > 0). Optional further tail lines:
`Jlya RT file: <path>` switches J̄_Lyα to **mode 1** (read an external RT profile
and use it directly, instead of the mode-0 parameterized field); `Deexc heat:
True` enables collisional de-excitation heating. As with all trailing lines the
scan is by keyword, so order and blank lines do not matter.

### Diagnostics

`write_excited_H` dumps `output/Excited_H.txt` (14 columns) from the converged
state: `r/Rp, T, nHI, ne, J̄_Lyα, n2s, n2p, S_proton(n=2), H_pe, H_dx, τ_Lyα`,
then the three **proton-budget** columns added in this phase —
**col12 ground-state photoionization rate** (P_HI · nHI), **col13 collisional
ionization** (C · n_e · nHI), **col14 recombination** (α_B · n_e · nHII), all
volumetric [cm⁻³ s⁻¹]. The three rate coefficients are captured inside `ioniz_eq`
on the converged pass (globals `gph_ground_HI`, `cion_HI`, `arec_HII`), so the
budget is self-consistent with the solver rather than re-derived offline.

### Touched files

- `init/parameters.f90` — Phase-3a globals: `use_excited_H`, `T_star_eff`,
  `R_star`, `gamma2_bal`, `hpe2_bal`, `jlya_mode`, `jlya_rt_file`; per-cell
  feedback arrays `gph_balmer_HI`, `heat_balmer`; diagnostics `Jlya_arr`,
  `n2s_arr`, `n2p_arr`, `Sproton_arr`, `Hpe_arr`, `Hdx_arr`; and the proton-budget
  capture arrays `gph_ground_HI`, `cion_HI`, `arec_HII`.
- `radiation/excited_hydrogen.f90` **(new)** — `excited_H_update` (per-cell n=2
  solve, ξ factor, Γ₂, J̄_Lyα mode 0/1, source and heating), `n2_populations`,
  `gamma_n2_balmer`, `heat_n2_balmer`, `load_jlya_rt`, and `write_excited_H`.
- `radiation/ionization_equilibrium.f90` — inject `gph_balmer_HI`/`heat_balmer`
  under `use_excited_H`; capture the ground-state / collisional / recombination
  coefficients for the proton-budget dump.
- `files_IO/input_read.f90` — scan the stellar / `Jlya RT file` / `Deexc heat`
  tail lines; set `use_excited_H`, `jlya_mode`.
- `EXHALE_main.f90` — in-loop `excited_H_update` before `ioniz_eq`; post-loop
  refresh + `write_excited_H`.
- `Makefile` — `excited_hydrogen.f90` in `SRC`.

### Validation on WASP-121b (coupling ON vs OFF)

Two converged runs, same build, differing only in the `use_excited_H` toggle
(stellar tail lines present/absent). Both reach steady state via `du < du_th`
with **log₁₀ Mdot = 12.51** (unchanged from Phase 2 — the coupling does not
destabilize the energy solver). Diagnostics in
`WASP-121b/{analyze_excited_H.py, analyze_proton_budget.py, compare_on_off.py}`.

**Proton budget (cf. Huang Figs. 11/27), volumetric rates [cm⁻³ s⁻¹]:**

| r/Rp | S_n2 | S_ground | S_recomb | S_n2/S_ground | S_n2/S_source |
|---|---:|---:|---:|---:|---:|
| 1.10 | 5.1e4 | 6.6e6 | 6.6e6 | 0.008 | 0.008 |
| 1.30 | 1.2e4 | 4.3e5 | 4.4e5 | 0.029 | 0.028 |
| 1.50 | 4.9e3 | 7.6e4 | 8.1e4 | 0.064 | 0.060 |
| 1.80 | 1.9e3 | 1.5e4 | 1.7e4 | 0.125 | 0.111 |
| 2.00 | 2.0e3 | 5.9e3 | 7.9e3 | 0.34 | 0.25 |
| 2.01 | 2.9e3 | 5.1e3 | 7.9e3 | 0.57 | 0.36 |

Observations (for evaluation, not pass/fail):

- **Self-consistency:** at the base S_ground ≈ S_recomb (6.6e6 vs 6.6e6 at
  1.1 R_p), i.e. photoionization equilibrium holds; collisional ionization is
  negligible (T too low).
- **n=2 share:** the n=2 Balmer channel is a minor proton source below ~1.5 R_p
  (~1–6% of ground-state photoionization) but climbs to **~34–57% of the
  ground-state rate (~25–36% of the total proton source) near 2 R_p**, where the
  gas is already ionized and EUV is attenuated (ground-state photoionization
  depleted) while the thin-Lyα-fed n=2 channel is not. This is the "significant
  proton source approaching 2 R_p" behavior Huang describe, here concentrated near
  2 R_p rather than flat across 1–2 R_p.

**Heating budget (cf. Huang Figs. 10/26):** Balmer photoelectric heating is
present in the budget, peaking at **~1.8% of total heating** at the Lyα
photosphere (τ_Lyα ~ 1 near r ~ 2 R_p) and ~1e-4 globally — a small but nonzero
entry, dominating nowhere, as expected.

**ON vs OFF temperature and ionization:** the coupling makes the gas slightly
**more ionized** (up to ~24% less neutral H near 2 R_p, Δx_HII up to ~1.1e-2) and
slightly **cooler** (up to ~179 K at r ~ 2 R_p, mean ~42 K). The cooling is
**physical, not numerical**: the extra n=2 photoionization converts HI→HII,
lowering nHI; the dominant EUV photoelectric heating scales with nHI, so ~24% less
nHI at r ~ 2 R_p means correspondingly less local EUV heating there, which the
small Balmer H_pe does not offset.

### Caveats and gaps

- **Parameterized J̄_Lyα (mode 0).** The 0.1·F_LyC/Δν_D field with 1/(1+τ_Lyα)
  attenuation is a placeholder for the Phase-3b Lyα RT solve. Because our n=2
  contribution is peaked near 2 R_p rather than flat across 1–2 R_p, if Huang's
  Figs. 11/27 are flatter/larger across the band the parameterized field is
  under-pumping n=2 in 1.1–1.8 R_p; mode 1 (external RT J̄_Lyα) is the path to test
  this once a Lyα RT result exists.
- **Convergence stop is momentum-uniformity** (`du < du_th`), not a thermal/
  ionization steady-state criterion, so a small residual convergence-path
  sensitivity rides on top of the physical ON–OFF difference.
- **De-excitation heating off by default** (overlaps the existing H I coex
  cooling); left behind `incl_deexc_heat`.

### Two analysis gotchas (cost false alarms; recorded so they are not repeated)

- `Hydro_ioniz.txt` columns 6/7 (heat/cool) are **already physical**
  (`write_output` multiplies by q0); do **not** multiply by q0 again in analysis.
- `Ion_species.txt` column 1 is `r/Rp`, so **nHI = col2, nHII = col3** — the
  recombination sink uses col3, not col2.

### Status and relation to the Huang plan

**Phase 3a is complete and the proton-source comparison is user-accepted**
(2026-06-05). The H(n=2) model is fully coupled into the ionization balance and
the energy equation; the n=2 photoionization proton source reproduces the
qualitative Huang Figs. 11/27 behavior and the Balmer photoelectric heating
appears in the budget without destabilizing the solver. Next: **Phase 3b** — a
Lyα radiative-transfer solve to replace the parameterized J̄_Lyα with a
self-consistent J̄_Lyα(r) (feeding mode 1).

## 9. Phase 2 refinement: Fe I line cooling (NIST f-values + Van Regemorter)

*Added 2026-06-05.*

Fe I line cooling was one of the two coolants left deferred when the Phase 2 gate
closed (§7, "Caveats and deferred"). This update adds it. CHIANTI v11 has **no
Fe I model atom**, so — following Huang et al. (2023) §2.5 — Fe I cooling is built
from NIST permitted-line oscillator strengths plus the Van Regemorter
collision-strength approximation, with the lower levels populated in Boltzmann
equilibrium over the metastable manifold. The new coolant flows through the same
`beta_esc * ne * n_FeI * Lambda_FeI(T)` machinery as the other Phase-2 metals.

### Atomic data and method (sources, not fabricated)

The cooling coefficient per `(n_e * n_FeI_total)` in the optically thin limit
[erg cm³ s⁻¹] is

```
Lambda(T) = sum_lines  f_l(T) * q_lu(T) * dE_lu ,
  f_l(T) = g_l exp(-E_l/kT) / U(T)                 Boltzmann lower-level fraction
  U(T)   = sum_{even levels, E < e_cut} g exp(-E/kT)
  q_lu   = 2.16 a^-1.68 exp(-a) T^-3/2 f_lu         Van Regemorter (Huang Eq. 10),
           a = dE/kT,  f_lu from NIST Aki + g's + lambda.
```

- **Raw data:** the auditable NIST dumps `fe1_nist_lines.tsv` (permitted E1 lines:
  Aki, g_i, g_k, E_i, E_k) and `fe1_nist_levels.tsv` (every level g, E for the
  partition function), fetched by `fetch_fe1_nist.py`.
- **Metastable manifold = EVEN-parity levels below the lowest odd level**
  (z⁷D°, 19351 cm⁻¹ = 2.40 eV). Below that energy there is no lower level for an
  E1 decay, so the even levels are genuinely metastable and hold the Fe I
  population. This is **not** a plain energy cut: the Fe II treatment can use an
  energy cut because its cut conveniently lands at the first odd level, but for
  Fe I a pure energy cut would wrongly include the radiatively-depopulated z⁷D°
  odd levels. The default cut yields **17 even levels and 828 permitted cooling
  lines**.
- Upper levels of the cooling lines are odd (E1) and assumed to decay
  radiatively. The builder is `cooling_data/fe1_cooling.py`; it emits the
  41-point `cool_logL_FeI(NCOOLT)` Fortran block on the ATES grid
  (log10 T = 3.0–5.0, dlogT = 0.05).

Coefficient values (per `n_e n_FeI`, erg cm³ s⁻¹): `Lambda(5e3 K) = 1.27e-22`,
`Lambda(1e4 K) = 1.23e-20`, `Lambda(2e4 K) = 1.31e-19` — the rising, then
flattening shape of Huang Fig. 5.

### Touched files

- `cooling_data/fetch_fe1_nist.py`, `cooling_data/fe1_cooling.py` **(new)** — raw
  NIST fetch and the Boltzmann-manifold + Van Regemorter builder.
- `radiation/Cool_coeff.f90` — `cool_logL_FeI(NCOOLT)` table constant; new
  `cool_FeI(T, out)` (interpolates the table via `interp_cool_table`); dispatch
  cases `FeI` in `cool_coeff_metal` and `cool_coeff_metal_scalar`. The table-header
  comment block now documents Fe I alongside the other coolants.
- `init/species_table.f90` — `mion_iscool(FeI) = .true.` (canonical ion index 25),
  so the `eval_cool` metal-cooling loop and the per-channel breakdown both pick up
  Fe I automatically.

### Validation on WASP-121b

Re-ran the converged full-metals WASP-121b equilibrium (same SED/grid as the
Phase 2 gate) with Fe I cooling on; the per-channel `Cooling_breakdown.txt` then
gives the Fe I share directly. Diagnostics:
`WASP-121b/analyze_FeI_cooling.py` (left panel = the ported `Lambda_FeI(T)`
coefficient read straight from `Cool_coeff.f90` vs the Fe II coronal curve,
Huang Fig. 5 style; right panel = per-channel cooling vs radius, Huang Fig. 10
style), figure `WASP-121b/FeI_cooling_validation.png`.

Observations (for evaluation, not pass/fail):

- **Fe I is a negligible coolant** along the whole profile: its peak share of the
  total cooling is **3.0e-5 at r = 1.06 R_p (T = 8303 K)**, and it is far smaller
  elsewhere (~3e-8 of total at the base, where Fe II dominates). Adding Fe I does
  **not** perturb the gate-band Fe II → Mg II handoff.
- **Mdot unchanged** at log10 Mdot = 12.51 (Fe I adds no measurable cooling), zero
  NaNs.
- This is the expected outcome: at WASP-121b base temperatures iron is almost
  entirely Fe II (Fe I fraction < 0.11, §5), and Fe I's permitted lines have high
  excitation energies, so its Boltzmann-suppressed cooling cannot compete with the
  Fe II / Mg II forbidden + resonance lines. Fe I is carried for completeness and
  to mirror Huang's coolant list, not because it shapes the WASP-121b energy
  budget.

### Status and relation to the Huang plan

Fe I line cooling is **implemented and its coefficient matches the Huang Fig. 5
shape/magnitude**; the WASP-121b run confirms it is radiatively negligible there,
surfaced as a figure for the reader's judgment (no pass/fail asserted). One of the
two Phase-2 deferred cooling items is now closed. The remaining one — the Fe II
density-dependent (IR fine-structure) low-temperature refinement — is §10.

## 10. Phase 2 refinement: density-dependent Fe II line cooling (SE table)

*Added 2026-06-05. (Implemented step-by-step; this section is updated as each
step lands — currently: the 2-D table generator is built and verified.)*

This is the second deferred Phase-2 cooling item (§7): replacing the **coronal**
Fe II line-cooling coefficient with a **density-dependent** one that follows the
Fe II level populations from the coronal (low-`n_e`) limit into the collisionally
saturated (high-`n_e`, LTE) limit. It directly affects the cool, dense base of the
WASP-121b thermosphere, where the coronal table overestimates Fe II cooling by
~10⁴×.

### Why the coronal table is wrong at the base

The Phase-2 Fe II table is the optically-thin **coronal** limit: only the ground
level is populated and every collisional excitation is assumed to be followed by a
radiative decay, so the cooling is strictly ∝ `n_e` (per `n_e n_FeII` it is a
function of `T` alone). That assumption fails at the cool, **dense** base of a hot
Jupiter thermosphere, because the low-`T` Fe II cooling there is carried almost
entirely by **forbidden** metastable / a⁶D fine-structure lines (A ≲ 5e-3 s⁻¹,
far-IR 26–87 µm), whose critical densities

```
n_crit = A_ul / sum_c n_c k_ul^c  ~  1e4 - 1e7 cm^-3   (electron collider)
```

sit far **below** the base `n_e ~ 1e8 cm⁻³`. Those lines are therefore
collisionally **saturated** (thermalized to their LTE populations); the coronal
coefficient overestimates their cooling by ~`n_e / n_crit`.

### The density-dependent coefficient (multilevel statistical equilibrium)

`cooling_data/fe2_cooling.py` solves the full CHIANTI multilevel
statistical-equilibrium (SE) problem for the Fe II level populations at each
`(T, n_e)` — 175 levels, electron-collisional excitation/de-excitation
(Burgess–Tully descaled Upsilon from the `.scups` file) plus radiative decay
(`.wgfa` A-values), normalized to Σ_i n_i = 1 — so the cooling smoothly
interpolates between the two limits:

```
Lambda_vol / n_FeII = sum_{u>l} n_u A_ul dE_ul          [erg s^-1 per ion]
Lambda_eff(T, n_e)  = (Lambda_vol / n_FeII) / n_e        [erg cm^3 s^-1]
```

At low `n_e`, `Lambda_eff -> Lambda_coronal(T)`; at high `n_e` it falls off ~`1/n_e`
as the numerator saturates. This `Lambda_eff` is exactly the per-`(n_e n_FeII)`
coefficient ATES already tabulates, so the existing assembly
`cool_M = beta_esc * n_e * sum_i n_ion(i) * Lambda(i)` carries it unchanged — the
`n_e` in the prefactor times `Lambda_eff` gives `beta_esc * n_FeII *
(cooling per ion)`, i.e. the correct saturated base cooling.

**Why electron-only SE is sufficient at the (mostly neutral) base.** In the
saturated limit the upper levels reach their **Boltzmann (LTE) populations
regardless of which collider thermalizes them**, and the per-ion cooling becomes
collider- and `n_e`-independent. At the base `n_e ~ 1e8` already exceeds `n_crit`
for the dominant forbidden lines, so the electron-only SE already reaches that LTE
plateau; adding the (more abundant) neutral-H collider would only push the few
still-unsaturated higher-`n_crit` lines further toward the same LTE limit — a
second-order correction. An explicit H I / proton collider channel is flagged as a
future refinement in the module.

### Step 1 — the 2-D table generator (done)

`fe2_cooling.py` gained an **optimized** table builder `lambda_eff_table(Tgrid,
negrid)` and Fortran emitters (`fortran_grid_1d`, `fortran_block_2d`,
`emit_fortran_table`; run `python3 fe2_cooling.py --table`). The optimization is
required because the naive per-cell solve rebuilds a cubic spline for every one of
the ~4300 transitions at every grid cell: instead the descaled Upsilon is
evaluated **once per transition over the whole T-grid** (`upsilon()` accepts an
array T), and the rate matrix is split `M = n_e · Q_coll(T) + A_rad`, so the
ne-loop only reassembles and solves. The full 41×29 table builds in ~3 s and
matches the naive per-cell solve to **1.5e-15** (machine precision).

The emitted table is on the ATES temperature grid (`NCOOLT = 41`, log10 T =
3.0–5.0, dlogT = 0.05) × an electron-density axis `cool_logne(NCOOLNE = 29)`,
log10 `n_e` = 0–14, dlog = 0.5 (brackets the WASP-121b run's log10 `n_e` ≈
7.9–9.7 with the coronal plateau below and LTE saturation above). It is stored as
`cool_logL_FeII_ne(NCOOLT, NCOOLNE) = log10 Lambda_eff`.

Two built-in sanity checks (printed to stderr by `--table`):

- **Low-`n_e` edge reproduces the coronal table.** At `n_e = 1 cm⁻³` the SE
  coefficient equals `cooling_effective(pop='coronal')` to **0.1%** across
  T = 1e3–1e5 K — so swapping the 1-D coronal table for the `n_e = 1` column of the
  2-D table is a no-op, and the Phase-2 gate's coronal Fe II validation (vs Huang
  Fig. 6, within ~1.5×) carries over unchanged.
- **Base suppression.** At the WASP-121b base (`T ≈ 2239 K`, `n_e = 1e8`) the
  coronal coefficient `2.80e-20` drops to the SE value `1.79e-24` — a
  **15,663× suppression**. The suppression is strongest at low T (forbidden lines)
  and weak at high T (T ≳ 3e4 K, where the resonance lines with large A dominate
  and stay near-coronal): e.g. only ~1.4× at T = 3e4, n_e = 1e8.

### Step 2 — the 2-D table and bilinear interpolator in `Cool_coeff.f90` (done)

The generated table and its interpolators are now in
`src/modules/radiation/Cool_coeff.f90`:

- **Table.** `integer, parameter :: NCOOLNE = 29`, the `cool_logne(NCOOLNE)`
  axis (log10 `n_e` = 0–14, dlog = 0.5), and the data block
  `cool_logL_FeII_ne(NCOOLT, NCOOLNE)` (a `reshape([...], [NCOOLT, NCOOLNE])`
  with T on the inner index, `n_e` on the outer), pasted straight from
  `fe2_cooling.py --table`. The `n_e = 1` column matches the existing 1-D
  `cool_logL_FeII` to the fifth decimal (first node −19.78002 vs −19.77975).
- **Bilinear interpolator.** `interp_cool_table_2d(logL2d, T, ne, out)` —
  log-log bilinear in (log10 T, log10 `n_e`). The T axis reuses the existing
  uniform `cool_dlogT` scheme; the `n_e` axis uses the uniform `dlogne = 0.5`.
  Both axes clamp to the nearest edge outside the table (so `n_e` below the
  grid gives the coronal plateau, above it the LTE branch). A scalar twin
  `interp_cool_table_2d_scalar(logL2d, Ts, nes)` mirrors it bit-for-bit for the
  per-cell post-process solve.
- **Coolant wrapper.** `cool_FeII_ne(T, ne, out)` calls the 2-D interpolator on
  `cool_logL_FeII_ne`. The 1-D `cool_FeII` is kept as the low-`n_e` edge / the
  scalar fallback; its header `!To Be Checked` note was rewritten to record that
  the density-dependent table now supersedes the coronal-only rate in the
  cooling assembly.

Builds clean (gfortran, full `make`); the new routines are additive and not yet
wired, so behavior is unchanged until Step 3.

### Step 3 — wiring the density-dependent Fe II into the cooling assembly (done)

The 2-D coefficient is now wired into both cooling paths so the local `n_e`
selects the coronal → LTE-saturated branch automatically:

- **Main solver (`util_ion_eq.f90`, `eval_cool`).** `n_e` is already computed at
  the top of `eval_cool` (`call calc_ne`). Immediately after the canonical
  `c_metal(:,i)` fill loop, a short loop locates the Fe II ion by name and
  overrides its column with `call cool_FeII_ne(T_K, ne, ...)`. Nothing else
  changes: the existing assembly `cool_M = beta_esc * n_e * Σ_i nm(:,i)·c_metal(:,i)`
  then multiplies by `n_e`, which **cancels the 1/`n_e`** baked into
  `Lambda_eff = (Σ_u n_u A_ul dE_ul)/n_e`, leaving the physically correct
  per-Fe II-ion cooling. Because the saturated (LTE) populations are
  collider-independent, the electron-only SE solve is exact in that limit; at
  low `n_e` `Lambda_eff` reduces to the coronal rate, so optically-thin
  upper-atmosphere cells are unaffected. The override flows automatically into
  the `use_2lev_cool` branch (which already reads `c_metal(:,26)`) and into the
  per-channel `cool_chan` diagnostic (which reads `c_metal(:,i)`), so the
  Cooling_breakdown dump stays an exact decomposition.
- **Post-process per-cell solve (`T_equation.f90`).** This path computes its own
  local `n_e` (from H, He, and the frozen metal stages) and sums the metal
  coolants with the scalar dispatcher. The Fe II term is special-cased to call
  the new `cool_FeII_ne_scalar(TT, n_e)` (a thin wrapper over
  `interp_cool_table_2d_scalar`), so the temperature the post-process converges
  to balances the same density-dependent Fe II cooling `eval_cool` applies.
- **`species_table.f90`.** The `!To Be Checked` Fe II comment was rewritten to
  state that `eval_cool` now applies the density-dependent SE coefficient (no
  longer "deferred").

Builds clean (gfortran, full `make`, all four touched modules relink).

### Step 4 — WASP-121b run and the cooling budget (numbers for the reader)

A fresh WASP-121b run (solar metals, `pp_metals 1`, otherwise the same
`input.inp` as the Fe I-cooling run) was compared against the snapshot of the
previous (coronal Fe II) run kept in `output_coronalFeII/`. Both share the same
grid, so `Cooling_breakdown.txt` columns are differenced directly
(`analyze_FeII_density_correction.py` → `FeII_density_correction_compare.png`).

**At the dense base (`r = 1.0`, `T = 2358 K` fixed by the lower BC).** The Fe II
line cooling drops from the coronal `2.23e-4` to the SE `1.71e-8`
erg cm⁻³ s⁻¹ — a **≈1.3e4× suppression**, consistent with the offline
fixed-`(T, n_e)` check in Step 1 (15,663× at `T = 2239 K`, `n_e = 1e8`). Because
coronal Fe II carried **99.1%** of the base cooling, the *total* radiative
cooling at the base falls **≈122×** (`2.25e-4 → 1.84e-6`); Fe II's share there is
now **0.9%**. The previous run's post-process printed
*"35 of 502 cells fell back to eq T (stiff base band)"*; with the spurious base
overcooling removed that message is **gone** — the per-cell T solve no longer
stiffens at the base.

**Transition region 1.15–1.4 R_p (Huang Fig. 10).** The band-integrated
metal-line cooling share in the new run is

| coolant | share |
|---|---|
| Mg II | 44.0% |
| O II  | 18.9% |
| Fe II | 13.4% |
| O I   |  7.3% |
| C II  |  6.3% |
| Ca II |  5.1% |

so **Mg II is now unambiguously the dominant metal coolant** across the
transition region, with Fe II demoted to a ~13% contributor — the qualitative
ordering Huang's Fig. 10 shows. (In the old coronal run Fe II still topped the
budget out to ~1.2 R_p, the bottom-right panel's red curve.)

**Global response.** Removing the base overcooling raises the heating efficiency,
so the steady-state mass-loss rate rises from `log10 Mdot = 12.51` to
**`12.80` g/s** (≈2×). The self-consistent wind is correspondingly denser and,
above the base, slightly **cooler** (e.g. `T(1.05)`: `7810 → 5744 K`) — the freed
energy goes into mechanical escape and adiabatic expansion rather than being
radiated, the expected energy-limited behavior. Whether this Mdot and the
transition-region budget match Huang's WASP-121b quantitatively is left to the
reader's judgment against the paper; the figure and table above are surfaced for
that comparison, not asserted as a pass.

**Caveat on the figure's panel (b).** The "suppression factor" there is the
*ratio of Fe II cooling between the two runs*, not a fixed-condition coefficient
ratio. At the base (BC-fixed T) it equals the coefficient suppression (~1.3e4×);
above ~1.2 R_p the coronal and SE coefficients converge (the hot resonance lines
dominate and stay near-coronal), so the residual ratio there reflects the
changed — denser — wind profile, not the coefficient.

## 11. Post-process refinement: re-solve metal ionization in the advected wind (`pp_metals 2`)

*Added 2026-06-05.*

The `_adv` profiles that `TPM.py` reads are built by `post_process_adv.f90`, which
re-solves the **per-cell temperature** at the advection-corrected structure (the
wind is advected, then its `T(r)` is re-converged against heating/cooling). Until
now the metal *ion* densities carried into those profiles were handled in one of
two ways, selected at runtime by `pp_metals` in `metals.inp`
(→ `pp_metal_mode`; see §8 caveats):

- `pp_metals 0` — **metal-free legacy** `_adv` (the pre-metals behavior).
- `pp_metals 1` — **frozen equilibrium** metals (the default): the converged
  equilibrium metal stages are carried into `_adv` unchanged. They therefore
  track the *equilibrium* structure, not the (slightly hotter, more ionized)
  advected wind.

This step implements the deferred third option:

- `pp_metals 2` — **re-solve**: the metal ionization balance is recomputed in
  each cell at the advection-corrected H/He and the post-process temperature, so
  the metal stage split tracks the advected wind self-consistently.

### Implementation — reuse the coupled residual with H/He pinned

The metals are pure equilibrium species (no advection equations of their own): in
the main loop they are solved *coupled* with H/He through `ion_system_HeH_metals`
+ `hybrd1` (photo + collisional ionization + recombination + charge exchange,
no time derivative). "Re-solve" therefore means re-evaluating that same
equilibrium balance at the post-process local `T` and the advected H/He — **not**
a new physics path.

Rather than duplicate the residual, mode 2 reuses it verbatim through a thin
wrapper, `ion_system_metals_pp` (new, in `post_process_adv.f90`):

```fortran
subroutine ion_system_metals_pp(N_in,x,fvec,iflag,params)
  ...
  call ion_system_HeH_metals(N_in,x,fvec,iflag,params)  ! full coupled residual
  fvec(1) = x(1) - pp_xHII_fix     ! pin HII  fraction
  fvec(2) = x(2) - pp_xHeII_fix    ! pin HeII fraction
  fvec(3) = x(3) - pp_xHeIII_fix   ! pin HeIII fraction
end subroutine
```

The wrapper overwrites only the three H/He rows with identity equations that pin
H/He to their advection-corrected fractions (stored in module-level
`pp_xHII_fix` / `pp_xHeII_fix` / `pp_xHeIII_fix`), while leaving the
`2·n_melem = 20` metal-stage rows untouched. So `hybrd1` moves only the metal
unknowns; H/He stay fixed at the advected wind. Two properties fall out for free:

- **Self-consistent `n_e`.** The reused residual computes the electron density
  *internally* from the pinned H/He plus the metals it is currently solving —
  exactly as the equilibrium solver does — so the re-solve sees the same
  `n_e(metals)` feedback, not a frozen `n_e`.
- **Charge exchange with the advected H.** The Kingdon & Ferland H↔metal
  couplings inside `ion_system_HeH_metals` see the *pinned* H⁺/He⁺ densities, so
  the metal balance reacts to the advected hydrogen ionization, not the
  equilibrium one.

### The per-cell sweep in `post_process_adv`

Mode 2 adds one sweep over cells, placed after the post-advection
`calc_ne` and before the photoheating refresh, inside the same `k`-loop that
re-converges `T`, so metals, H/He and `T` are driven to a *joint* fixed point.
Per cell `j` (skipping `nh(j) ≤ 0`):

1. `cx_set_cell(T_K(j))` — charge-exchange rate coefficients at the local `T`.
2. Build the per-element coefficients (canonical order) from the existing
   post-process rate arrays: `meg_g0 = P_m` (photoionization), `meg_b0 =
   aion_m_pp` (collisional ionization), `meg_a1 = rec_m_pp` (recombination), and
   the second-stage `g1/b1/a2` where `melem_top ≥ 2`; hand them to
   `set_metal_coeffs`.
3. Pin `pp_xHII_fix = nhii/nh`, `pp_xHeII_fix = nheii/nhe`,
   `pp_xHeIII_fix = nheiii/nhe`; set `params(7)=nh`, `params(8)=nhe`.
4. Seed `sys_x` with the pinned H/He fractions plus the current metal split and
   call `hybrd1(ion_system_metals_pp, …)`.
5. **Conserve element totals**: the pre-solve totals `nm_tot_pp(:,im)` (summed
   over each element's stages once, before the `k`-loop) are multiplied by the
   re-solved fractions to write back `nm_w`. Only the *stage split* changes; the
   element abundance is untouched.

Supporting edits: the `use` block pulls in `ion_system_HeH_metals`,
`set_metal_coeffs`, and `cx_set_cell`; `params` was widened from `dimension(25)`
to `dimension(60)` to match the residual's dummy length (only entries 1–11 are
consumed by the metal solve); module-level `pp_x*_fix` and locals
`nm_tot_pp` / `meg_*` / `i0,top,im` were added. Builds clean (gfortran, full
`make`).

### Clean A/B and what changes

Because the **equilibrium** solve in the main RK loop is mode-independent and
deterministic, a full `pp_metals 2` run reproduces `Hydro_ioniz.txt`,
`Ion_species.txt`, and `Cooling_breakdown.txt` **byte-identically** to the
`pp_metals 1` run (verified with `diff -q`); *only* the `*_adv.txt` profiles
differ. So the A/B isolates exactly the post-process metal treatment. The mode-1
default is preserved in `output_ppm1_dd/`, the mode-2 run in
`output_ppm2_resolve/`, and the canonical `output/` is left at the mode-1 default
(`metals.inp` reset to `pp_metals 1`); `analyze_metal_resolve.py` differences the
two and writes `metal_resolve_compare.png`.

**Singly-ionized fraction X⁺/X_tot in the `_adv` wind (eq frozen → re-solve):**

| r/R_p | T_eq | T_adv | C | O | N | Mg | Ca | Na | Fe |
|---|---|---|---|---|---|---|---|---|---|
| 1.00 | 2358 | 2358 | .966→.959 | .000→.000 | .000→.000 | .866→.849 | .809→.826 | .706→.662 | .934→.992 |
| 1.05 | 5744 | 5798 | .981→.984 | .002→.002 | .014→.015 | .957→.957 | .636→.588 | .945→.953 | .997→.996 |
| 1.10 | 7779 | 7984 | .961→.964 | .016→.012 | .119→.132 | .924→.913 | .653→.600 | .995→.996 | .993→.989 |
| 1.15 | 8729 | 9044 | .912→.916 | .096→.058 | .378→.403 | .918→.905 | .661→.588 | .998→.999 | .992→.991 |
| 1.20 | 9651 | 10035 | .784→.774 | .366→.255 | .505→.466 | .872→.864 | .609→.538 | .999→.999 | .925→.956 |
| 1.30 | 10180 | 10274 | .600→.600 | .340→.339 | .345→.344 | .657→.658 | .380→.380 | .999→.999 | .387→.401 |
| 1.40 | 10624 | 10680 | .462→.462 | .232→.232 | .234→.234 | .500→.499 | .245→.245 | 1.000→1.000 | .182→.186 |
| 1.48 | 10980 | 11049 | .374→.374 | .173→.173 | .175→.175 | .404→.403 | .179→.179 | 1.000→1.000 | .114→.116 |

The corrections are **modest and concentrated in the 1.1–1.25 R_p ionization
transition**, where `T_adv` runs a few hundred K hotter than `T_eq`:

- **Fe II** is *higher* across the 1.18–1.27 R_p shoulder (the blue curve sits
  above red in the figure; e.g. 0.925→0.956 at 1.20 R_p), the largest meaningful
  density change being **+9.5%** in n(Fe II) at r = 1.232 R_p.
- **Ca II** is *lower* by ~5–7 points through 1.05–1.25 R_p (0.661→0.588 at
  1.15 R_p); max **−11.9%** in n(Ca II) at r = 1.191 R_p.
- **Mg** and **Na** are essentially unchanged (≤2% in n(Mg II); the Mg/Na panels
  overlap), and **every element converges back to the equilibrium split in the
  outer wind** (≳1.27 R_p) where `T_adv` rejoins `T_eq`.
- The post-process printed **no** stiff-base T-fallback message in mode 2 (the
  density-dependent Fe II cooling of §10 already removed the base overcooling).

**Tiny-denominator artifacts (not physical).** The per-element "max rel change"
scan reports **O +6538%** at r = 1.000 and **N +17.9%** at r = 1.013. These are
at the dense base where O and N are essentially *fully neutral* (O II fraction
≈ 0.000 in the table), so n(O II) is a vanishingly small number; a tiny absolute
shift (n⁺ ≈ 9.7e2 → 6.4e4 cm⁻³, against an O total ~1e8) reads as a huge
*relative* change but is negligible in absolute and in any observable. They are
flagged here as denominator artifacts, not a physical re-solve effect.

**Reader's judgment, not a pass.** The re-solve makes the metal stage split in the
`_adv` profiles consistent with the advected wind rather than the equilibrium
structure, with the physically expected sign (more ionization where the advected
wind is hotter) and magnitude (single-digit-percent, transition-localized). The
default remains `pp_metals 1` (frozen): the deliberate, documented choice (the
equilibrium split is the converged one, and `TPM.py`'s lines are dominated by H
and He). Mode 2 is provided for studies that need the metal ionization to track
the advected wind. The table and `metal_resolve_compare.png` are surfaced for
that comparison, not asserted as a validation pass.

---

## 12. Phase 3b: in-code Lyα radiative transfer (J̄_Lyα via escape probability)

Phase 3a drove the H(n=2) population with a **parameterized** Lyα mean intensity
`J̄_Lyα = 0.1 F_LyC/Δν_D / (1 + τ_Lyα)` (mode 0). That estimate peaks the n=2
population near 2 R_p and **under-pumps** it across the transmission-relevant
1.1–1.8 R_p band, leaving a `jlya_mode = 1` hook (`load_jlya_rt`) to ingest a real
RT field. Phase 3b supplies the real `J̄_Lyα(r)`.

Three findings reshaped the original plan (which was a Monte Carlo with an
offline two-pass and a β-on-cooling term):

1. **Method → escape probability, not Monte Carlo.** The WASP-121b column is
   extremely thick at line center: `τ₀(Lyα) ≈ 5.5×10⁷` at the base (Voigt
   `a ≈ 5×10⁻⁴`, so `a·τ₀ ≈ 10⁴ ≫ 1`). A faithful acceleration-free MC needs
   `~τ₀ ≈ 10⁷–10⁸` scatterings/photon (intractable), and core-skipping is both
   forbidden and incompatible with the core-dominated pumping field. In the
   `a·τ₀ ≫ 1` regime the **Neufeld (1990)/Harrington analytic wing-escape**
   solution is accurate, so we use a local escape-probability source function.
2. **Computed in-code, not offline.** Because the formula is cheap, J̄ is evaluated
   **every timestep** inside the run (`jlya_mode = 2`, the same lagged-explicit
   decoupling as the Phase-3a Balmer feedback), so a single run is
   self-consistent. The `jlya_mode = 1` option to read an external (e.g. real MC)
   J̄ profile is kept.
3. **No β-on-cooling term.** Reading Huang §2.5/Fig 7 directly: the "Lyα cooling
   10× lower" is the **Black (1981) vs CHIANTI collisional-excitation
   rate-coefficient** difference, **not** photon trapping. ATES already uses the
   CHIANTI-consistent Cen `coex_rate_HI` (§7), so there is no 10× to apply. The
   genuine trapping feedback on the energy balance is the **Phase-3a H(n=2)
   heating** (photoelectric `Hpe` + collisional de-excitation `Hdx`), which
   strengthens automatically once J̄ is correct. The scaffolding β-on-`coex` hook
   was therefore reverted.

### The escape-probability field (`lya_rt.f90 : jlya_escape_prob`, mode 2)

Per cell, from the converged `T`, `n_HI = n_1s`, `n_HII`, `n_e`:

- Doppler width `Δν_D = ν_Lyα √(2kT/m_H)/c`; Voigt `a = (A_2p1s/4π)/Δν_D`.
- Line-center optical depth `τ(r) = ∫_r^∞ n_1s σ₀ dr'`, `σ₀ = C_Lyα/Δν_D` (same
  trapezoid as mode 0).
- **Wing escape probability** `β(r) = π^(-1/4) √(a/τ)` (≤1): a photon escapes once
  it diffuses to `x* = (aτ/√π)^{1/2}` where the wing optical depth is 1
  (Neufeld/Harrington).
- **Internal sources** (Huang's two): `P = α_B n_e n_HII + q_{1s2p} n_e n_HI`
  (recombination cascade + electron-impact 1s→2p). The escape-probability closure
  `J̄ = S(1−β)` with the trapped pile-up `n_2p = P/(A_2p1s β)` gives
  `J_int = (2hν³/c²)(g_1s/g_2p) P (1−β)/(A_2p1s β n_1s)`. The **(1−β)** factor is
  essential: `J_int → S` when thick and **→ 0 when thin** (β→1), instead of the
  source function `S ∝ 1/n_1s` diverging as the gas ionizes in the outer wind (which
  spuriously raised J̄ outward).
- **Stellar beam (trapped).** Lyα is resonantly *scattered*, not destroyed, so the
  beam is **not** attenuated by `e^{-τ}`. The broad stellar line (half-width
  `Δv_s ≈ ±70 km/s`) penetrates through its wings: a stellar photon at Doppler
  offset `x` reaches depth `r` if `|x| > x₁ = √(aτ/√π)`, so for a Gaussian stellar
  profile of width `X_s = Δv_s/v_th` the penetrating fraction is
  `T_*(r) = erfc(x₁/(√2 X_s))` (smooth; a flat profile gives a kink). The
  profile-averaged mean intensity uses the **stellar line width** `Δν_⋆ = ν₀ Δv_s/c`
  (≈ constant), **not** the planetary `Δν_D` — for a stellar field flat across the
  narrow planetary line `J̄ = ∫ J_ν φ dν ≈ J_ν(ν₀)`; dividing by `Δν_D` would make
  J_star rise spuriously as the wind cools. A bounded trapping buildup
  `E = 1 + (boost−1)(1−β)T_*` enhances the penetrating beam in the thick wind
  (→ `boost` where `T_*~1, β≪1`; → 1 in the thin outer wind and at the deep
  penetration edge), avoiding both the `1/β` over-count (which blew the n=2 source
  to 6300× ground-state) and a flat cap (which over-enhances the dense base). Thus
  `J_star = ξ F_⋆ T_* E /(4π Δν_⋆)`, with `ξ` the dayside dilution.
- `J̄_Lyα = J_int + J_star`.

New `input.inp` lines (all default-off → byte-identical to Phase 3a):
```
Jlya escape-prob: True               # jlya_mode = 2 (in-line escape-probability RT)
Stellar Lya flux [erg/cm2/s]: 1.0e5  # incident stellar Lyα at the planet (Case D ×0.35)
Lya stellar halfwidth [km/s]: 70     # broad stellar line Δv_s (penetration + Δν_⋆)
Lya stellar boost [-]: 5             # bounded trapping-buildup cap (left un-tuned)
Deexc heat: True                     # Phase-3a collisional de-excitation heating
```

A 3b-0 scaffolding step first verified the plumbing was inert (WASP-121b
byte-identical to the pre-change binary; frozen reference
`WASP-121b/output_pre3b_baseline/`).

### WASP-121b validation setup (spherical Case A, matched to Huang)

Huang's Fig 11 is **Case A** (spherical, no RLOF), so the comparison must match his
geometry **and** radial frame — which exposed a setup error. The run had been
Roche-truncated at ~2 R_p, and ATES's base (r=1) sits at **1 μbar**
(`n_0=3.1×10¹²`, T=2358 K); but Huang defines `R_p` at the **4 mbar transit radius**
and places his 1 μbar hydro bottom at **~1.25 R_p** (Fig 9 top-pressure axis /
Fig 12 green "lower atm." band; the earlier "1.46 R_p" note was wrong). The setup
was therefore corrected to **spherical Case A** with `R_0` set to Huang's 1 μbar
physical radius (`R_0 = 1.25 × 1.766 = 2.2075 R_J`, planet mass unchanged → correct
base gravity) and outer radius 15.2 (= 19 R_p,Huang); ATES `r` maps to Huang's frame
as `R_Huang = 1.25 r`. The Case-A `Ṁ = 0.077 M_p/Gyr` is within ~1.5× of Huang's
0.052 (Table 3).

### Trend vs Huang Fig 11 (qualitative — exact match is not expected)

Huang's J̄_Lyα is a **plane-parallel Monte Carlo** RT; ours is an analytic
escape-probability formula, so the two cannot agree quantitatively — a meaningful
quantitative comparison is only sensible once **all** Phase-3b features (incl. the
velocity term) are in. Qualitatively (`WASP-121b/lya_caseA_huang11.py`), both show a
**significant inner peak** (ATES ~0.6 @1.35, Huang ~0.5 @1.6 R_p) and comparable
magnitude (0.05–0.6). ATES has a **dip at ~1.8 R_p** and **rises gently outward**
where Huang gently declines; the outward rise traces to the **n2s recombination
cascade** in the Phase-3a H(2s/2p) model (n=2 stays recombination-populated while
ground-state H vanishes in the ionized wind) — **not** J̄, which after the two fixes
above is physically sound (J_int peaks then declines; J_star ≈ const outward).
H(n=2) heating stays ~1% of total with `Hdx ≳ Hpe` near the n_2 peak (cf. Huang
Fig 10).

> **Two J̄ bugs fixed during the trend check (2026-06-06):** (1) `J_star` used the
> planetary `Δν_D` instead of the stellar line width `Δν_⋆`, making J̄ rise
> spuriously as the wind cooled; (2) `J_int` was missing the `(1−β)` escape closure,
> so the source function diverged (∝1/n_1s) in the ionized outer wind. Both are now
> corrected.

### Wind-velocity (Sobolev) escape — implemented, negligible here

`β_tot = 1 − (1−β_static)(1−β_Sobolev)` with `τ_S = C_sob n_1s/|dv/dr|`
(`C_sob = (πe²/m_e c) f λ`) is now in `jlya_escape_prob` (the outflow `v(r)` is
threaded through `excited_H_update`). For WASP-121b it changes J̄ by **<1%**: where
it could help (the thick inner wind) `τ_S ~ τ_static ~ 10⁷` so `β_Sob ~ 10⁻⁷`, and
the thin outer wind already has `β_static ≈ 1`. The velocity gradient and the Lyα
opacity scale together (both ∝ n_1s), so the wind velocity opens **no new escape
route** for this extremely thick line — a physical result, not a tuning. **Phase 3b
is thus feature-complete** (escape-probability J̄ = internal + trapped stellar beam;
de-excitation heating; Sobolev escape).

**Deferred (per the user's trend-first guidance):** the post-all-features
quantitative comparison and the n2s outward-rise revisit (the residual rise is the
Phase-3a 2s recombination cascade, not the RT). The `boost`, `Δv_s`, and `π^(-1/4)`
prefactor are left un-tuned.

---

## 13. Phase 4: RLOF / tidal potential (Case B validated; Case D blocker = escape radius, transonic IC added, geometry-limited)

Phase 4 adds the Roche-lobe-overflow (RLOF) tidal physics that drives Huang's
enhanced mass loss — **Case B Ṁ = 0.32, Case D = 1.03 M_p/Gyr** — the regime that
matches the observed transit depths (spherical Case A underestimates them).

**ATES already has the 1D Roche potential** (`grav_field.f90`, the
`spherical_domain=.false.` branch: planet + stellar-tidal + centrifugal), the
L1-truncated domain (`r_max` = Roche lobe), and the tidal momentum source
(`Source.f90`) — Caldiroli's ATES-v2 Roche mode. So the substellar tidal
hydrodynamics is already in place; no new hydro was needed.

### The per-case base gravity (Huang Table 3)

The key to each case is the 1 μbar **base radius** R₀. Huang's Table 3 lists
**log g @ 1 μbar = GM_p/r₁μbar²** per case — exactly the R₀-setting quantity
(R₀ = √(GM_p/g)):

| Case | tidal | M_p (M_J) | dT | XUV | Lyα | Fe | log g@1μbar | → R₀ (R_J) | Ṁ (Huang) |
| :-- | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| A | N | 1.1824 | 0 | 1 | 1 | 1 | 2.84 | 2.058 | 0.052 |
| B | Y | 1.1824 | 0 | 1 | 1 | 1 | 2.70 | 2.418 | 0.32 |
| D | Y | 1.1204 | +350 | 0.74 | 0.35 | 4 | 2.62 | 2.581 | 1.03 |

R₀ grows A→B→D as the tidal force + the dT+350 lower-atmosphere heating puff the
1 μbar level outward (weaker base gravity), driving the RLOF.

### Case B — validated

WASP-121b in Roche mode (tidal on) + the corrected base + metals + Lyα mode-2 gives
**Ṁ = 0.337 M_p/Gyr ≈ Huang Case B (0.32), within ~5%**. With the base-radius fix
both cases agree (Case A 0.077 ≈ 1.5× Huang's 0.052; Case B 0.337 ≈ 1.05× Huang's
0.32). ATES's existing Roche potential thus captures the RLOF enhancement — no
geometric-reconstruction code was needed for the Case B Ṁ gate.

### Case D — real blocker found (escape radius), transonic IC added

Case D's base (R₀ = 2.581 R_J, log g = 2.62) sits at **~76% of the L1 distance**:
the L1-truncated Roche domain reaches only **r_max = 1.353 R_p**. Re-investigating the
earlier "won't launch" report uncovered the **actual** blocker, and it is *not* the
initial condition.

**Root cause — the escape radius lies outside the L1-truncated domain.** `define_grid`
sets the constant-momentum convergence window `[j_min:N]` from the **escape radius**
(`Escape radius [R_p]:`, default **1.50**). It scans for the first cell with
`r(j) ≥ r_esc`; if *no* cell qualifies — which happens precisely when
`r_esc > r_max` — the `DO` loop falls through leaving `j = N+Ng+1`, so **j_min > N**
and the window `[j_min:N]` is **empty**. Fortran `maxval`/`minval` over a zero-size
array return **∓HUGE**, so `mom_max = −HUGE`, `mom_min = +HUGE` → `du = Inf`, and the
time-derivative norm `dtu = −HUGE`. Since `dtu = −HUGE < dtu_th`, the loop reports
**"steady state reached" at step 0** — a *false* convergence that has nothing to do
with the wind. The earlier bracket ("works R₀ ≤ 2.30, fails R₀ ≥ 2.40") is exactly
**r_max crossing the 1.50 escape radius at R₀ ≈ 2.37** — confirming the escape radius,
not the IC, as the discriminant. (So the earlier "isothermal IC cannot launch the
deep-RLOF wind" diagnosis was wrong.)

**Fixes.** (1) A guard in `define_grid.f90`: if `j_min > N`, re-anchor the window to
mid-domain (`r ≥ 1 + ½(r_max−1)`) and print a loud WARNING telling the user to set a
smaller escape radius — so this can never again silently masquerade as convergence.
(2) For Case D, set `Escape radius [R_p]: 1.20` (inside the 1.353 domain). With a valid
window the du/dtu diagnostics become finite and meaningful, and **Case D launches and
runs** — from *both* the hydrostatic and the transonic IC. (The earlier
`du = |(mom_max−mom_min)/max(mom_min, 1e-30)|` floor at `EXHALE_main.f90:308` is retained
as defense-in-depth against a genuine transient zero-momentum cell, but it was *not* what
unblocked Case D — with an empty window `mom_min = +HUGE`, not zero, so the floor never
engaged.)

**Transonic-wind IC (implemented, `Transonic IC: True`).** `set_IC.f90` gained an
optional transonic isothermal-wind initial profile, solved from the ATES potential via
the algebraic (Bernoulli) integral of the steady isothermal-wind equation,

$$\tfrac12 v^2 - c^2\ln v - 2c^2\ln r + \phi(r) = B,\qquad c^2=(1+\delta p_{bc})/\rho_{bc},$$

with `B` fixed by the sonic point (`v=c_s` where `Dφ(r_c)=2c²/r_c`, located by
bisection). At each radius the cubic-like relation has a subsonic (`v<c_s`) and a
supersonic (`v>c_s`) root; the wind takes the subsonic root for `r<r_c` and the
supersonic root for `r>r_c` (robust 1-D bisection on each monotone branch), and ρ then
follows from mass conservation `ρv r² = Ṁ` pinned to the base BC. Because that IC
already satisfies `ρv r² = const`, `du ≈ 0` at step 0 would itself trip the
momentum-constant exit before the cold wind heats up; the flag therefore **auto-enables
`force_start`** so the heating develops first. If no interior sonic point exists it
**falls back to the hydrostatic IC**. Default **off** ⇒ all existing cases are
byte-identical (Case B re-verified at **Ṁ = 0.337**).

**Case D outcome — launches, but the L1-truncated 1-D domain undershoots.** With the
escape-radius fix the wind launches and relaxes, but it **never reaches a sonic point
inside the domain**: the heated metal wind is **deeply subsonic everywhere**
(Mach ≲ 0.05 — rising to ~0.050 at r ≈ 1.25, then *falling* to ~0.040 at L1). With a
subsonic outer boundary the zero-gradient outflow BC is only marginally well-posed, and
`ρv r²` stays non-constant (peaks ~1.7× the L1 value), so the run settles into a
**breathing limit cycle**: after the ~10² heating transient damps, du keeps oscillating
in ~0.04–0.8 (touching but never crossing the strict `du_th = 0.02`) for as long as it
runs. The mass flux is correspondingly time-variable — the quasi-steady **Ṁ oscillates
over ~0.3–0.6 M_p/Gyr**, well *below* Huang's 1.03, the deficit set by the L1 truncation
cutting off the still-accelerating subsonic wind before it reaches the sonic point. The
transonic IC behaves identically to the hydrostatic IC here: both launch the same
subsonic flow. (Warm-starting via `load_IC` to accelerate
convergence reintroduces a base NaN — the documented grid re-read inconsistency — so it
is not a workaround.)

**Conclusion.** Case D's deep RLOF (Ṁ = 1.03) is **not reachable in ATES's 1-D
L1-truncated Roche mode**: the wind is subsonic at L1, so the truncated domain
structurally undershoots. Capturing it needs the 3-D geometric treatment (Phase 5) or a
domain extended past L1 — *not* a different IC. The transonic IC and the escape-radius
guard are correct, necessary infrastructure (they remove the false convergence and let
the deep-RLOF wind launch and run with meaningful diagnostics), but they are not
sufficient on their own for Case D's full mass-loss rate.

### Geometric reconstruction → Phase 5

The substellar → terminator equipotential mapping, the triaxial radii
(R_px/R_py/R_pz), and the **4/9 mass-loss conversion** feed the transmission spectrum
(Phase 5, §14); they also supply the geometry that the subsonic L1-truncated Case D wind
needs.

## 14. Phase 5: velocity-broadened transmission spectrum (in progress)

Phase 5 is the headline deliverable — the NUV/optical **transit-depth spectrum** the whole
upgrade exists to match (Huang Table 3 / Figs 17–24). It extends the existing spherical
post-processor `TPM.py` in two stages: **5a** adds the metal resonance lines (spherical,
gated on Case A); **5b** replaces the spherical geometry with the 3-D Roche-equipotential
reconstruction + velocity broadening (gated on Case D). Pure Python post-processing — no
Fortran change.

**Validation targets** — these are the **effective transit radius `R_eff/R_star`**, not an
absorption percent (the expected-results doc header states *"단위: R_p/R_*"*; the "~30%" in
its prose is a loose gloss of 0.30). Our code returns the line-center absorption fraction
`h`, which maps to the same quantity via `R_eff/R_star = √((R_p/R_star)² + h)`:

| Line | Case A | Case D | Observed |
| :-- | :-: | :-: | :-: |
| Mg II λ2796 (4 Å) | 0.182 | 0.302 | 0.309 ± 0.036 |
| Ca II K λ3934 | 0.199 | 0.278 | 0.281 ± 0.009 |
| Hα λ6563 | 0.201 | 0.185 | 0.186 ± 0.003 |
| Hβ λ4861 | 0.174 | 0.135 | 0.143 ± 0.005 |
| Na D2 λ5890 | 0.152 | 0.147 | 0.147 ± 0.002 |

### 5a — metal resonance lines (Mg II, Ca II, Na I D), spherical — done

`TPM.py` gained Mg II λ2796, Ca II K λ3934, Na I D2 λ5890 as an **isolated, appended block**
(the validated He 10830 / Lyα / Hα / Hβ pipeline above it is untouched). They are resonance
lines whose lower level is the ion ground state; at ~10⁴ K the excited fine-structure levels
are Boltzmann-negligible, so the lower-level density is the ion density itself
(`n_lower ≈ n_ion`). The ion densities live past the H/He block of `Ion_species_adv.txt`, in
`mion_fsp` order (`species_table.f90`): numpy columns **17 = Mg II, 23 = Ca II, 25 = Na I**.
The block reuses the existing spherical chord / Voigt (`wofz`) / trapezoid-τ / disk-average /
instrument-convolution machinery (`resonance_depth()`), so it is mechanically identical to
the He/Lyα lines.

**Case A result (line-center % / 4 Å-band %), and diagnosis:**

| Line | line-center | 4 Å-band | R_eff reached | Huang Case A |
| :-- | :-: | :-: | :-: | :-: |
| Mg II 2796 | 26.1% | 2.3% | **4.2 R_p** | 0.182 |
| Ca II K | 5.5% | 0.33% | 1.9 R_p | 0.199 |
| Hα | 3.9% | — | 1.6 R_p | 0.201 |
| Na D2 | 0.8% | 0.06% | 0.7 R_p | 0.152 |

Two real gaps surface, both informative: (1) **Mg II already reaches Huang's extent**
(line-center 26% → R_eff 4.2 R_p, vs Huang's clustered ~3.6 R_p) — the metal machinery is
sound — but **Ca II / Na I / Hα are too optically thin** (R_eff 0.7–1.9 R_p), a
line-specific ionization / abundance / n=2-population question, not a geometry one. (2) **The
lines are too narrow**: the Case A wind reaches only **29 km/s**, so a 26%-deep core washes
out to ~2% over a 4 Å (±215 km/s) bin, whereas Huang's lines stay deep *across* the bin —
they are velocity-broadened to hundreds of km/s by the fast Roche outflow + rotation, i.e.
**Phase-5b physics**. (The existing Hα tool shows the same shortfall, so it is systematic,
not introduced by the metals.)

**Units — resolved.** The mismatch above is one of *units*, not physics: Huang's table is
`R_eff/R_star` (radius ratio), while the numbers quoted here are absorption percent. The
Case A spherical run over-extends (its 15 R_p domain lets Mg II reach 4.2 R_p →
`R_eff/R_star = 0.51`, vs Huang's 0.182); the proper test is the L1-confined Roche run in
`R_eff/R_star` units (see 5b).

### 5b — Roche-equipotential reconstruction + velocity broadening — implemented

New module **`roche_recon.py`**: the 3-D reconstruction. Assuming a uniform thermodynamic
state on Roche-equipotential surfaces, the field at any (x,y,z) equals the 1-D substellar
state at the equivalent radius `r_eff` where `phi_sub(r_eff) = phi_3D(x,y,z)`. The
dimensionless Roche potential (R_p units, `GM_p/R_p` energy units) mirrors `grav_field.f90`:

  φ(x,y,z) = −1/r − q/r_s − (1+q)/(2a³)·[(x − x_cm)² + y²],  q = M_*/M_p, a = orbit [R_p].

Setting the tidal terms → 0 recovers φ = −1/r (spherical), so the triaxial reconstruction
reduces to the 5a spherical case. The module provides `roche_phi`, the L1 root of
`dφ_sub/dr = 0`, a monotone `phi_sub → r_eff` inverse (`ReconMap`), the triaxial radii, and
the Eggleton Roche-lobe radius. `TPM.py` gained a `geometry='triaxial'` switch
(`triaxial_depth()`): the transit LOS runs along +x; the state at each 3-D point is the
substellar state at `r_eff`; and the **per-sector LOS velocity** `v_sub(r_eff)·x/r − Ω·y`
(wind + tidally-locked rotation, Huang Eq. 16) is integrated over **20 angular sectors**
(vectorized over wavelength).

**Demonstrated on the converged Case B (Roche) run.** Geometry self-consistent: L1 =
1.574 R_p, R_RL/L1 = 0.714 ≈ 2/3, and the equipotential is correctly tidally triaxial
(R_px = 1.05 > R_py = 0.919 > R_pz = 0.893: substellar bulge, polar compression).

**Result — the model reproduces Huang's transit radii.** Converting the Case B triaxial
depths to `R_eff/R_star` (the table's unit):

| Line | Huang A | Huang D | **Case B (this work)** |
| :-- | :-: | :-: | :-: |
| Mg II | 0.182 | 0.302 | **0.234** |
| Ca II | 0.199 | 0.278 | **0.222** |
| Na D2 | 0.152 | 0.147 | **0.163** |
| Hα | 0.201 | 0.185 | 0.239 |
| Hβ | 0.174 | 0.135 | 0.232 |

Case B is the *intermediate* model (between Case A and Case D), and the **metals (Mg II,
Ca II, Na) land right between Huang's Case A and Case D** — the lobe-confined reconstruction
correctly reproduces his transit radii in both magnitude and trend. The Balmer lines come out
~25% high (Hα 0.239 vs the expected ~0.19), consistent with the 5a finding that the n=2 / T is
slightly over-predicted.

**Note (correction):** an earlier draft of this section claimed the lobe-confined model
"undershoots Huang because it omits the beyond-L1 escaping tail." That was wrong on two
counts: (1) Huang's reconstruction is *also* lobe-confined (equipotential mapping to the
Roche-lobe boundary, §1.8) — neither model includes a beyond-L1 tail, so it cannot be the
difference; and (2) the apparent shortfall was a **units error** — comparing our absorption
*percent* against Huang's `R_eff/R_star`. In the correct units the two agree (table above).
The small triaxial-vs-spherical reduction (terminator compression `R_py < R_px`) is a minor
geometric effect; Huang's area-preserving normalization `R_py·R_pz = R_p²` would offset it.

### 5b — direct Case A comparison: hydro validated, but the metals are over-extended

A clean **spherical Case A at the correct 1 μbar base R0 = 2.058 R_J** (Table 3 log g = 2.84;
the old `output_caseA` used R0 = 1.766, the 4 mbar transit radius) was run to convergence.

- **Hydro / base validated:** Ṁ = **0.056 M_p/Gyr ≈ Huang's Case A 0.052** (the old wrong
  base gave a different Ṁ). So the base radius and the wind are right.
- **Metal transit radii — consistent once compared like-for-like.** Comparing in Huang's own
  definition (Mg II is quoted in a **4 Å NUV bin**; Ca II / Na are optical line-center):

  | Line | Huang A | this work | |
  | :-- | :-: | :-: | :-- |
  | Mg II (4 Å) | 0.182 | **0.220** | ✓ 1.2× |
  | Na D2 (LC) | 0.152 | **0.175** | ✓ 1.15× |
  | Ca II K (LC) | 0.199 | **0.306** | 1.5× high |

**The "over-extension" was a comparison-definition artifact, not a physics problem.** An
earlier draft reported Mg II R_eff/R_star = 0.60 (3.3× Huang) — that compared the *line-center*
R_eff against Huang's *4 Å NUV-bin* value. In the same 4 Å bin, Mg II = 0.22 ≈ Huang's 0.182.
Na matches too; only Ca II is genuinely ~1.5× high.

**The metal ionization is correct and is *not* the lever.** All channels are wired and
physical (Verner+1996 σ, Badnell RR+DR, Voronov collisional, Huang+2023 charge exchange;
`a_tau = 0`, no spurious attenuation). The Mg II/III balance is a close competition between a
weak photoionization (Γ ≈ 2e-3 s⁻¹, weak because σ_MgII ≈ 0.24 Mb at threshold) and RR·n_e;
DR (resonances ~5×10⁵ K) and collisional ionization (`exp(−15eV/kT)` at 10⁴ K) are
negligible. Crucially, a test showed **×30 more Mg II photoionization moves R_eff only
0.57→0.32** — the resonance line is so optically thick that ionizing it harder cannot confine
it, so the ionization is a dead end as a lever (and is not the cause once the bin is matched).

**Status:** the full Phase-5 machinery (metal + Balmer lines, spherical & triaxial geometry,
wind + rotation broadening, triaxial-radii prediction) is implemented and runs end-to-end;
the hydro/Ṁ reproduces Huang (Case A 0.056, Case B 0.337); and the metal transit radii match
Huang once compared in his own bins — Mg II (4 Å) 0.22 vs 0.182, Na 0.175 vs 0.152. Remaining
items are minor: **Ca II ~1.5× high** (line-center; plausibly Ca II→Ca III or the Ca
abundance), the slightly-high Balmer n=2 (TPM-internal J̄_Lyα estimate), the `R_py·R_pz = R_p²`
normalization (methodological), a terminating Case D run, and the `phase5_transmission.ipynb`
packaging. Note: comparisons must use Huang's per-line definition (4 Å bin for the NUV metals,
line-center for the optical lines) — line-center vs 4 Å differs by ~3× for the optically-thick
Mg II.

## 15. Post-process advection guard at a "breathing" base (2026-06-07)

**Symptom.** The spherical-vs-Roche example figure (user manual Fig. 4, WASP-121b)
showed a sharp temperature spike at r ≈ 1.1 R_p in the advection-corrected (`_adv`)
Roche profile, with a sawtooth band below it. The converged equilibrium (`eq`)
profile is smooth there — so the spike is a post-processing artifact, not hydro.

**Diagnosis — two compounding causes (both confined to the base):**
1. **Breathing, inflowing base (1-D limitation).** WASP-121b nearly fills its Roche
   lobe → the 1-D wind is subsonic at L1 and the dense lower atmosphere recirculates
   with small *negative* (inflow) velocities and a stagnation point (v = 0) at
   r ≈ 1.1. Intrinsic to spherical symmetry, not a bug.
2. **Spurious hot root in the metal-cooled post-processor.** The `_adv` solve
   re-solves the per-cell energy/ionization balance assuming an outflow (upwinds
   from the next-inner cell). In the dense base the non-monotonic metal line-cooling
   curve gives the energy equation a second, spurious *hot* root that MINPACK can
   land on; the upwind coupling cascades one bad cell outward → sawtooth + spike.
   The pre-existing 2× band guard did not catch the ~1.4× spike.

**Fix (option c — eq fallback where the flow is not a clean outflow).** In
`post_process_adv.f90`, wherever `v ≤ 0` (the breathing inflow base) the
post-processor now *skips* the advection correction and keeps the converged `eq`
ionization and temperature. Three guards, all gated on `pp_metal_on` (metals-off is
byte-identical): the no-He and He H/He advection loops, and the per-cell T solve.
This removes the spurious root **and** breaks the upwind cascade at its source,
while the genuine `_adv` correction is retained in the outflow above the stagnation
point.

**Validation (PP-only re-run on the converged Case B IC).**
- `eq` reproduced: max |ΔT|/T = 2.6e-3, max |Δv| = 0.015 km/s; log₁₀ Ṁ = 13.38
  unchanged (≈ 0.34 M_p/Gyr); Mg II reproduced to 1.6% in the line-forming region.
- `_adv` base oscillation eliminated: T sign-changes in [1.02, 1.10] R_p dropped
  **36 → 0**; spike gone; the eq→`_adv` handoff at v = 0 is continuous (0.7% step)
  and the few-percent wind correction is preserved.
- Spherical Case A (v > 0 everywhere) is untouched — the guard never triggers.

Cause 1 (the breathing base) remains a 1-D limitation but is now represented there
by the smooth `eq` solution; resolving the inflow itself would need 2-D/3-D.

**Documentation updates (this session).** User manual gained §1.3 *The gravitational
potential*, §3.5 *Spectral coverage for metal photoionization* (neutral-metal
thresholds, lowest K I at 4.34 eV/2856 Å), §8 *Atomic and molecular data* (by-process
source table; code is purely atomic/ionic — no molecular data), and the §4
breathing-base caveat. Fig. 4 was regenerated from the corrected `_adv` output.


## 16. He 2³S triplet + metals simultaneous solve, Lyα-flux guard, and convergence workflow (2026-06-08)

**Motivation.** Until now the He metastable triplet (He 2³S, "HeITR") and the
trace-metal ionization network could not be solved in the same run: both placed an extra
unknown at `x(4)` of the coupled ionization residual, and `input_read.f90` allocated the
metal unknowns only `if (thereis_metals .and. .not. thereis_HeITR)`. Modelling HD189733b
with the triplet *and* solar metals requires both at once.

**Merged solver (factual).** New module
`nonlinear_system_solver/System_HeH_TR_metals.f90` (`ion_system_HeH_TR_metals`, 84
SLOC). The H/He/triplet rows (`x1=HII, x2=HeII, x3=HeIII, x4=He 2³S`) are taken verbatim
from `System_HeH_TR`; the metal rows are taken verbatim from `System_HeH_metals` and
shifted to `x(5+2*(e-1))`. The two blocks share one electron density
`ne = nHII + nHeII + 2 nHeIII + Σ(nX+ + 2 nX++)` (the neutral triplet is excluded), and
charge exchange is added through `cx_add_to_fvec`. Solved with `hybrd1`. Supporting
edits: `charge_exchange.f90` gained a settable base index `cx_metal_base` (default 4; the
merged dispatch sets 5); `input_read.f90` sets `N_eq = 4 + 2*n_melem` when both are
active; `ionization_equilibrium.f90` uses a generalized metal base index (`mbase = 4`, or
5 with HeITR) plus a both-on dispatch branch; the new source was added to the `Makefile`.
Only the equilibrium *solver* needed merging — the radiation grid, cross sections,
photo-heating, IC setup and post-processing already handled both species.

**Regression & physical check.** Metals on / triplet off is byte-identical to the
pre-merge binary over 8048 steps (new branch not taken). With both on, the coupled solve
is clean (no NaN) and was exercised end-to-end on WASP-121b with the full stack (triplet
+ metals + excited-H + Lyα), reaching a flat converged solution (peak He 2³S fraction
~5×10⁻⁶).

**Lyα escape-probability input guard (factual).** `input_read.f90`: if
`Jlya escape-prob: True` (`jlya_mode=2`) but `Stellar Lya flux [erg/cm2/s]:` is missing
or ≤ 0, the run now stops with an explanatory message rather than silently dropping the
stellar beam (in mode 2 `J̄_Lyα = J̄_int + J̄_star` with `J̄_star ∝ F_Lya_star`, so
`F_Lya_star = 0` would make the halfwidth/boost settings silent no-ops).

**Convergence threshold and automated two-stage reconstruction (factual).**
- `du_th` is now a runtime variable (was a compile-time `parameter`) with default
  `1.0e-3` — the original ATES-Code-main value. EXHALE had carried `2.0e-2`, at which
  a run flagged "converged" can still show ~2% spatial spread in the supersonic mass flux
  ρvr².
- A new runtime variable `du_th_plm` plus the input line `du_th [PLM,WENO3]: <du_plm>
  <du_final>` enables an *automatic* two-stage run: `EXHALE_main.f90` starts in PLM and
  switches `rec_method` to WENO3 once `du < du_plm` (or PLM stalls), then converges at
  `du < du_final`. This automates the manual PLM → `Load IC` + WENO3 workflow recommended
  in the ATES README. If `du_plm ≤ du_final` (or the line is absent) the run is
  single-stage.
- `CFL` is now runtime and settable via the input line `CFL: <value>`.

**Other source-term changes (kept, with caveats).** `energy_semi_implicit.f90` now damps
the implicit energy update with `dF/dT = 1 + c·|dC/dT|` (was `max(0, dC/dT)`), adding
damping on the falling cooling branch; it removed one oscillation in testing but did
*not* by itself fix the convergence behaviour below. `Cool_coeff.f90` now interpolates
the cooling tables with monotone PCHIP (C¹) rather than linear (C⁰); this was implemented
to test an interpolation-kink hypothesis (not supported, below) and is retained as a
smoothness improvement.

**Convergence — provisional observations (not conclusions).** *The following is a
research log, deliberately tentative; several intermediate "non-convergence" readings in
earlier sessions were later traced to methodology, not physics.* Two artifacts are
recorded so they are not repeated: an over-short step cap (weak-wind cases judged
"floored" near 25k steps while `du` was still descending — HD209458b later reached
`du < 2e-2` near 50k steps), and a run-script bug that killed the backgrounding subshell
rather than the `EXHALE.x` child, so runs continued orphaned and the reported `du` was a
premature snapshot (use `pkill -x EXHALE.x`). Because of these, the earlier framing that
"metal line cooling destabilizes the wind" should be treated as **not established**. With
the two-stage strict-`du_th` pipeline we currently see, across planets:
- **WASP-121b** (inflated, strongly irradiated): converges quickly and flat (~0.16%
  supersonic mass-flux spread, vs ~2% for single-stage PLM). This case looks solid.
- **HD209458b**: reaches a small `du ~ 0.02–0.025` that then oscillates rather than
  dropping to 1e-3. *Tentatively* usable as a quasi-steady state: a He 2³S + 7-metal
  snapshot from the two-stage pipeline gives a smooth transonic wind with
  Tmax ≈ 6.8×10³ K and log₁₀ Ṁ ≈ 10.1 (median over the supersonic region), consistent
  with the literature ~10¹⁰ g/s; see `heitr_metals_dev/HD209458b_quasisteady.ipynb`.
- **HD189733b** (strongly bound, β₀ ≈ 192): `du` oscillates around ~2 in a sustained
  limit cycle and does not reach the threshold — the problematic case.

A plausible (but unproven) reading is that the difficulty tracks how weakly driven /
strongly bound the wind is, rather than the metal content per se.

**Base-breathing hypothesis and the fixes that did *not* work.** The HD189733b
oscillation resembles the known base-breathing behaviour, and the lower boundary
(`Apply_BC.f90`, `BC_component_constrho`) hard-pins the ghost density and pressure to
cold reservoir values with a one-way velocity valve `max(v1,0)` — a configuration that
can reflect acoustic waves. Three variants were tried and reverted:

| Test | Outcome |
| :-- | :-- |
| Lower CFL (0.6 → 0.2) | `du` floor barely moved (~2.8 → ~1.9) |
| Zero-gradient base velocity (drop the valve) | `du` still ~2.3 |
| Riemann-invariant base velocity `v = v1 + 2(c_b−c1)/(γ−1)` | did not help HD189733b **and** broke WASP-121b |

So the base-breathing reading is *consistent with* the observations but is **not
demonstrated**, and these simple BC variants are not solutions. Candidate routes if
pursued (both substantial, uncertain payoff): a characteristic / non-reflecting
(NSCBC-style) base BC; a steady-state Newton / BVP solver (which would not orbit a
time-marching limit cycle if a steady solution exists); or accepting HD189733b as a
marginal, time-averaged quasi-steady case. A separate, cheaper idea for *speed* (not the
oscillation) is local per-cell time-stepping, since the global `dt` is currently set by
the smallest base cell (`eval_dt.f90`). Full working notes:
`docs/heitr_metals_and_convergence_notes.md`; solver-option discussion:
`docs/numerical_methods.md`.

**Code-size note.** The merge adds one Fortran module (`System_HeH_TR_metals.f90`, 84
SLOC), so the 2026-06-07 file tally of 57 becomes 58; the remaining changes are small
edits to `charge_exchange.f90`, `input_read.f90`, `ionization_equilibrium.f90`,
`parameters.f90`, `EXHALE_main.f90`, `energy_semi_implicit.f90` and `Cool_coeff.f90`.

---

## 17. Metals in the bulk-gas mass/charge/particle budget (`eos_metals`)

*Added 2026-06-13.*

**Motivation.** Until this update the trace metals were strictly *passive*
species: they contributed to heating/cooling and were solved in the coupled
ionization equilibrium, but they were excluded from the mass density, the base
density `rho_bc`, and the electron/particle densities used everywhere
*outside* the ionization solver (EOS temperature, opacities, cooling calls,
energy update). The only self-consistent place was the charge balance *inside*
the MINPACK residuals (`System_HeH_metals.f90`). At solar abundances the
missing metal mass is ~1.3% of the H/He mass — small but systematic, and
inconsistent with solving the metal ionization to high precision. Recorded in
`TO_BE_DONE.md` on 2026-06-12; implemented 2026-06-13.

**Runtime interface.** New key `eos_metals 0|1` in `metals.inp` (documented in
`inputdata/metals.inp.example`). `1` (the new default) includes the metals in
the bulk-gas budget; `0` restores the legacy trace approximation exactly.
Metals-off runs (no `metals.inp`) are identical either way.

### What changed (factual)

- **Composition constants** (`input_read.f90`): the mass per hydrogen nucleus
  is now `mass_per_H = 1 + 4*HeH + sum(melem_ab*melem_A)` and the base density
  becomes `rho_bc = mass_per_H/(1 + HeH)` (in units of `m_H*n0`). Atomic
  weights live in `species_table.f90:melem_A`. The normalization `n0` keeps
  its H/He-nuclei meaning, so `n_H = n0/(1+HeH)` is unchanged. (Correction to
  the original `TO_BE_DONE` entry: `v0`, `p0`, `b0` are m_H-based
  normalizations and never contained a mean molecular weight — they were never
  affected.)
- **Mass density**: `calc_rho` takes an optional `nm` argument and adds
  `melem_A(elem)*nm` over all stages; passed from `ionization_equilibrium`.
- **Electron density**: `calc_ne` takes an optional `nm` and adds `stage*nm`;
  passed at every site where `nm` is available — `composition` (main-loop
  EOS), `ionization_equilibrium` (opacity/cooling), `util_ion_eq` (`eval_cool`
  + cooling dump), `energy_semi_implicit`, `excited_hydrogen`,
  `post_process_adv`. (The scalar `T_equation` already carried metal electrons
  via `pp_nm_cell`.)
- **Particle density**: `calc_ntot` likewise adds the metal nuclei, so the EOS
  `T = p/(n_tot + n_e)` sees them.
- **Ghost pressure / `dp_bc`**: the lower-BC ghost pressure is now
  `ntot_bc + dp_bc` with `ntot_bc = (1 + HeH + sum(melem_ab))/(1 + HeH)` — so
  T = T0 holds *exactly* at the base — and `dp_bc` includes the metal
  electrons at the first cell. `set_IC`/`load_IC` use `mass_per_H`/`ntot_bc`
  consistently. Full formulas: `EXHALE_BC_and_IC.tex`.

### Validation (WASP-121b regression matrix)

- `eos_metals 0` is **byte-identical** to the pre-change goldens for both
  regression cases (`wasp_full` 7296 steps, `wasp_he23off` 7300) — the legacy
  path is untouched.
- `eos_metals 1` (new default) converges cleanly (7145/7142 steps, ~2% fewer);
  profile shifts are at the expected ~1% mass scale: median |dn| ~ 2.4%,
  |dp| ~ 2.3%, |dT| ~ 0.9% (max |dT| = 132 K in the breathing-base region
  r ~ 1.09). The *number* flux n·v·r² drops 1.1% (log10 Mdot 13.31 → 13.30),
  which offsets the +1.2% mass per particle — the *mass*-loss rate is
  essentially unchanged.
- Goldens re-snapshotted under the new default; a repeat run reproduces them
  byte-identically (single-thread determinism intact).

**Still open (small).** The `use_2lev_cool` legacy cooling branch remains
AIOLOS-hard-coded (default off), and the per-cell scalar `T_equation` uses
`pp_nm_cell` only in the post-process path (the time-marching per-cell solve
receives n_e via the now-metal-aware `energy_semi_implicit`).

---

## 18. Automatic initial-condition selection (`IC mode: auto`)

*Added 2026-06-13.*

**Motivation.** The code now carries three IC families — cold hydrostatic
(default), transonic isothermal wind (added for the Phase-4 deep-RLOF Case D,
where the cold IC false-converges), and the hot-Parker warm seed (benchmarked:
no speedup) — and the choice was manual. Kubyshkina+2018 build a per-planet IC
automatically; the design study is `docs/auto_ic_design.md` and the literature
comparison `docs/code_comparison.tex`.

### The selector (one exact probe, no tunable threshold)

`select_IC_auto` (`set_IC.f90`) evaluates the *cold* base sound speed in code
units, `c2_cold = (ntot_bc + dp_bc)/rho_bc`, and calls the existing
`find_sonic` bisection for a critical point `phi'(r_c) = 2c²/r_c` in the
*actual* potential (Roche or spherical). If an interior sonic point exists,
the cold isothermal atmosphere already wants to blow a transonic wind, so the
transonic IC is selected; otherwise the cold hydrostatic IC is kept. One
criterion catches both regimes that need the transonic start: deep RLOF
(through the Roche/L1 topology) and low-gravity boil-off
(r_c ≈ b0/(2c0²) in a point-mass potential). The Jeans parameter b0
(Kubyshkina's Λ) is logged as a diagnostic but is *not* a decision variable;
the hot-Parker seed is never auto-selected.

**Runtime interface.** New input line
`IC mode: <cold|transonic|hot_parker|auto>`. The named modes simply set the
legacy flags; `auto` defers the decision to `set_IC`. The legacy keys
(`Transonic IC:`, `Hot Parker IC:`) take precedence over `auto`, and the setup
report (`ATES.out`) records the family actually in effect, with a
"(chosen automatically: IC mode = auto)" tag. Touched files: `parameters.f90`
(`ic_mode`), `input_read.f90` (key parsing), `set_IC.f90` (`select_IC_auto`),
`write_setup_report.f90`.

### Validation gates (all run 2026-06-13)

| Gate | Result |
| :-- | :-- |
| A-1 no `IC mode` key | byte-identical regression (defaults untouched) |
| A-2 `IC mode: cold` | byte-identical to A-1 |
| A-3 spherical tutorial, `auto` | b0 = 83.2, no interior r_c → cold; `IC_dump` byte-identical |
| B-1 Case D (deep RLOF), `auto` | r_c = 1.297 (just inside L1 = 1.353) → transonic, same as the manual key |
| B-2 boil-off planet (0.0189 MJ, 0.446 RJ, 1100 K) | b0 = 8.5, r_c = 5.2 (as predicted) → transonic; wind launches (initial residual 0.19 vs 14 cold) |
| B-3 `wasp_full` (Roche), `auto` | → transonic; du-stop in 5876 vs 7145 steps (−18%) |

**Honest speed verdict.** The −18% step count is *Roche/deep-RLOF-specific*
(the steady wind there is close to the cold isothermal transonic solution).
For the mild spherical boil-off case B-2 the cold control *also* launched and
converged first (log10 Mdot = 11.20), while the auto run relaxed monotonically
toward the same flux from ~2× above — the cold-c² Parker profile
*overestimates* the flux of a weakly heated wind. So auto is not harmful (same
attractor) but not faster there.

**Key byproduct (a trap for all Mdot work).** While adjudicating gate B-3 it
was found that the `du < 1e-3` stopping criterion is **path-dependent at the
~5% level**: the cold-start and transonic-start `wasp_full` runs both
"converge" by `du` but sit 4.8% apart in mass flux, each only 1–2% flux-flat.
Newton-finishing *both* states (`Load IC` + `Solver: Newton`) collapses them
onto the *same* fixed point (flux ratio 0.9996, flatness 0.005%). Practical
rule, now standing: **for quantitative Mdot always finish with
`Solver: Newton`**; never compare bare du-stops at the percent level.

**Deferred.** A Phase-C retry ladder (escalate cold → warm seed on NaN/stall)
is designed in `auto_ic_design.md` §6 but *not* implemented — the B-2 result
suggests failure-triggered escalation is the right shape if a real case ever
needs it. IC formulas and the selection logic are also documented in
`EXHALE_BC_and_IC.tex` §5 and the user manual.

## 19. Base boundary conditions and the flux-based convergence decision (2026-06-15)

**Context.** A systematic study of how the lower boundary and the stabilizers
interact with convergence on HD 209458b is collected in
`docs/EXHALE_BC_and_IC.tex` (test matrix, two stabilizer tables, and a
convergence-criterion history). This entry records the resulting default changes
to the code; the companion document holds the full evidence.

**Shapiro filter: off by default.** The optional Shapiro (1970) spatial filter —
ported from CETIMB (Koskinen et al. 2013a), whose use of it was found through
Huang et al. (2023), to damp the same gravity-unbalanced base sound-wave
instability — damps the transient "breathing" of the base cell, but on a clean
cold-start it also drives the wind into a slow infall and prevents convergence.
Isolated in a 2×2 (filter × base-valve) sweep, the filter — not the valve — was
the cause. It is therefore *off* by default (`shapiro_eps < 0`) and opt-in only
via `Shapiro filter: <eps> [every]`.

**Base boundary condition: density vs pressure.** The base can now be anchored
either on density (the legacy fixed n₀) or on a fixed base pressure, selected by
`Base BC: density|pressure [p_ubar]` (default density; `pressure` fixes p_base in
μbar and derives n₀). In practice the converged wind is nearly insensitive to
which anchor is used — the real lever on Mdot is the base *density* itself, not
the BC form.

**The flux-based convergence decision (history).** EXHALE measures two distinct
quantities at the base: the flux metric `du` (the radial spread of ρvr², which is
exactly the constant-mass-flux criterion of ATES, ΔMdot/Mdot < 1e-3, and of
CETIMB, F_c = ρvr² altitude-flat), and the steady residual ‖R‖ = ‖∂_t u‖ over
mass/momentum/energy, a *stricter* EXHALE-specific check. The convergence rule
evolved: before the Newton finish it was purely flux-based; the Newton work then
added the stronger ‖R‖ criterion; but the reference codes converge on flux alone.
The standing decision, adopted here, is therefore **flux-based** — `du < du_th` is
the acceptance test, and ‖R‖ is computed and reported *for reference only* (it
gates the stop only if `Resid tol:` is set).

**Volume-weighted residual norm.** When ‖R‖ *is* reported, the default is now a
cell-volume-weighted relative norm, Σⱼ|Rⱼ|·rⱼ²Δrⱼ / Σⱼ|uⱼ|·rⱼ²Δrⱼ (`resid_vol`,
settable `Resid norm: vol|Linf`), rather than the bare L∞ maximum. On the
non-uniform radial grid the L∞ max is dominated by the few tiny near-base cells
and overstates the residual; the volume weighting gives a physically
representative number. The infall diagnostic used in the test matrix was likewise
switched to a *radius-weighted* negative-velocity fraction (weighted by Δr over
the linear radius range), so that sub-10 m/s near-base noise is not counted as
bulk infall.

**du-keyed Newton/JFNK hand-off.** The production Newton finish (`Solver: Newton`)
now hands off from the marching warm-up to the JFNK steady solve when the *flux*
metric drops below `newton_du_switch` (default 1e-2), consistent with the
flux-based decision and adjustable via an optional third token,
`Solver: Newton [du_switch]`. This replaces the earlier residual-keyed hand-off
(‖R‖ < 5e-2, variable `newton_R_switch`, now removed); the two are equivalent in
practice, but keying on `du` avoids tying the hand-off to the stricter residual.
The JFNK solve still *targets* ‖R‖ < `resid_th` (1e-3); on failure it keeps its
best iterate and falls back to marching with the du-based stops re-armed.

**The full-physics residual floor (base j=1).** With He 2³S + metals, the
du-keyed hand-off plus the volume-weighted norm drove the JFNK residual from
~0.17 down to ~2e-3 — about two orders of magnitude better — but it then *stalls*
(`info ≠ 0`) on a localized momentum imbalance at the base cell j=1 (the
"breathing base"), short of 1e-3. The wind is nonetheless flux-converged (e.g.
HD 209458b full physics: log10 Mdot steady, ρvr² spread < 1%), which is a
converged run by the reference standard. Closing the last factor of ~2 in ‖R‖
requires physically damping the base momentum — the explicit-viscosity task at
the top of `TO_BE_DONE.md`.

**Touched files.** `parameters.f90` (`shapiro_eps`, `base_bc_mode`/`base_p_ubar`,
`base_v_massflux`, `resid_vol`, `newton_du_switch`; `newton_R_switch` removed);
`input_read.f90` (`Base BC:`, `Resid norm:`, `Solver: Newton [du_switch]`
parsing); `steady_residual.f90` (`residual_norms_vol`); `steady_newton.f90`
(`resid_relnorm` vol/L∞ branch); `EXHALE_main.f90` (du-keyed JFNK trigger,
always-on flux/residual diagnostic). Full detail and the test matrix:
`docs/EXHALE_BC_and_IC.{tex,pdf}`.

## 20. Tidal (Roche) potential with an extended outer boundary (2026-06-16)

**Motivation.** The two domain modes were tied to two different potentials.
`Domain mode: Spherical` extends the grid to a user `Outer radius` but *drops* the
stellar-tidal and centrifugal terms (pure `-b0/r`); the default Roche mode keeps
the full tidal potential but truncates at the inner Lagrange point L1. Neither
matched Yan et al. (2022) / Huang et al. (2023), which retain the tidal force
*and* integrate out to ~10 R_p (well past L1). For WASP-52b, L1 ≈ 2.53 R_p, so
Roche mode cut the domain far short of Yan's 10 R_p, truncating the outer
He 10830 / Hα absorbing layers.

**Change.** The outer-boundary extent is now decoupled from the potential choice.
In Roche mode (no `Domain mode` line → full tidal potential of `grav_field.f90`
retained), an explicit `Outer radius [R_p]: >1` sets `r_max` to that value instead
of the L1 radius:

```fortran
   else  ! Roche mode: full tidal potential retained
      r_max = (3.0*Mrapp)**(-1.0/3.0)*atilde      ! L1 (default)
      if (r_out_user .gt. 1.0d0) r_max = r_out_user  ! NEW: extend past L1
   endif
```

This reproduces the Yan/Huang "tidal force + extended domain" setup with a
two-line change and no new keyword (the existing `Outer radius` is simply honored
in Roche mode too).

**Caveats and usage.** Beyond L1 the Roche potential is past its maximum (net
outward force), so the 1-D radial flow there is a spherical approximation to an L1
funnel — the same approximation Yan/Huang make. The cold hydrostatic IC is invalid
past L1, so pair the extended domain with `IC mode: auto` (validated: the auto
selector finds an interior cold sonic point and seeds a transonic wind). The tidal
singularity of the `-b0*q/(atilde-r)` term sits at the star (r = atilde ≫ 10), so
a 10 R_p domain is numerically safe. `Domain mode: Spherical` (no tidal) remains
available for Huang Case-A-like comparisons.

**Touched files.** `input_read.f90` (honor `r_out_user` in the Roche branch of the
`r_max` assignment). Documented in `docs/EXHALE_user_manual.tex` (potential section
and the input-options table).

## 21. Photoionization cross-section benchmarks and a TOPbase high-energy extension of He 2³S (2026-07-01)

**Motivation.** The He 2³S (metastable triplet) photoionization cross section
`sigma_HeI23S` — the drain term that sets the He 10830 lower-level population — was
a four-segment broken power-law fit to Norcross (1971), hard-cut to zero above
59.2 eV. Two questions: (i) are EXHALE's H/He/He 2³S photoionization cross sections
consistent with modern atomic data, and (ii) does the unphysical 59 eV cutoff
matter? Both are settled in the new memo `docs/photoion_cross_sections.{tex,pdf}`
(reproducible via `docs/compare_photoion_cross_sections.py`).

**Benchmarks.** H I and He II use the exact hydrogenic form; He I uses the 2-term
ATES fit, which agrees with the Verner+1996 (VFKY96) single-shell fit to within 14 %
over 24.6–500 eV. The He 2³S broken PL reproduces its parent Norcross/p-winds table
to ~10 % and — checked against the **TOPbase / Opacity Project** R-matrix cross
section (extracted from the CDS server; raw 922-point table saved to
`data/topbase_HeI_2_3S_photo.dat`) — agrees to ~8 % in the smooth near-threshold
region (5–20 eV). TOPbase additionally resolves a 40–55 eV autoionizing-resonance
forest (peaks to 3223 Mb); their **net effect on the 2³S rate is only ~1–2 %** once
the few numerically under-resolved spikes are excluded (the raw table must *not* be
fed into a broadband rate integral — a quadrature artifact). Extending the cross
section above 59 eV raises the 2³S rate by only 0.1–1.2 %.

**Change.** `sigma_HeI23S` is extended with a minimal, fully-connected edit
(`src/modules/functions/cross_sec.f90`): the broken PL is kept unchanged up to the
45.6 eV bump; the **last power-law segment is re-aimed** so it runs continuously
from the bump peak to a power-law fit of the TOPbase smooth background at the
junction E_J = 70 eV,

```
sigma_tail(E) = 10^5.4623 * E^-2.9964  Mb      (E > 70 eV)
```

(slope of the last segment changes only −3.039 → −2.658; intercept derived in-code
for exact continuity). The old 59 eV cutoff is dropped, so 59–1240 eV now carries
the TOPbase tail instead of zero. The result is C⁰-continuous at *both* joins
(45.6 eV and 70 eV); the original pre-modification formula is kept in comments.
A unified full-range PCHIP node table and a tail-only power law are documented as
alternatives, but this is the recommended drop-in (smallest rate change, retains the
validated threshold-to-bump shape exactly).

**Validation.** A controlled before/after on WASP-52b (`fxuv0p25_he98_L1`,
`Include He23S? True`, ε Eri SED with X-rays), same IC, only `cross_sec.f90`
differing: the bulk wind is unchanged (Mdot log = 11.49 both), and the He 10830 line
changes by < 0.001 % in every metric (line-center depth 15.2675 % → 15.2675 %;
EW 0.0877 Å → 0.0877 Å) — below the Newton convergence noise. The extension is more
physically correct yet observationally harmless for He 10830, as predicted by the
+0.65 % rate estimate (even smaller here because the ε Eri SED is soft).

**Recombination cross-check.** A companion memo
`docs/recombination_coefficients.{tex,pdf}` (`compare_recombination.py`) documents
the H/He/He 2³S recombination: H II, He II, He III use the Hui & Gnedin (1997)
case-B fits (≈ Draine 2011 / Verner & Ferland to ~2 %); He II→He I splits into
Benjamin+1999 singlet/triplet (Oklopčić & Hirata 2018) when the 2³S network is on.
The He 2³S triplet rate is identical to the latest escape models (Benjamin+1999).
Newer nebular-community data (Porter+2012/2013, Del Zanna & Storey 2022) differ
~5–15 % but are not yet adopted by escape codes; H/He case-B needs no update.

**Touched files.** `src/modules/functions/cross_sec.f90` (`sigma_HeI23S`: re-aimed
last segment + TOPbase tail, 59 eV cutoff removed, original kept in comments). New
docs: `photoion_cross_sections.{tex,pdf}`, `compare_photoion_cross_sections.py`,
`data/topbase_HeI_2_3S_photo.dat`, `recombination_coefficients.{tex,pdf}`,
`compare_recombination.py`, and figures under `docs/figures/`.

## 22. He I (1¹S) photoionization → Verner+1996 by default, with ATES-rate switches (2026-07-01)

**Motivation.** The He I *ground-state* (1¹S, singlet) photoionization used the
legacy ATES two-term fit, which runs ~12–14 % below the Verner+1996 (VFKY96)
single-shell fit in the 50–100 eV band (`docs/photoion_cross_sections.tex`, §4).
Update it to Verner by default, keep the old behavior behind a switch, and add a
parallel switch for the recombination data.

**Singlet/triplet separation (clarification).** He photoionization in EXHALE is
already split by spin: `sigma_HeI` is the **1¹S ground (singlet)** cross section
(updated here to Verner — Verner's He I *is* the 1¹S ground state), while
`sigma_HeI23S` is the **2³S metastable (triplet)** cross section (Norcross/TOPbase,
§21). The 2¹S singlet metastable is not tracked (only 2³S, the 10830 lower level),
so the Verner update correctly touches only the singlet ground state.

**Change.** A single backward-compatible switch was added in `input.inp`:
`ATES_photoionization_rate: True` reverts He I (1¹S) photoionization to the legacy
two-term fit. **Default = Verner+1996** (`sigma_VFKY96` with He I params
`[E_th,E_0,σ_0,y_a,P,y_w,y_0,y_1] = 24.59, 13.61, 949.2, 1.469, 3.188, 2.039, 0.4434,
2.136`). No recombination switch is needed — see below.

**Comparison (WASP-52b, `fxuv0p25_he98_L1`, He23S active; same converged IC).**
Switching He I photoionization legacy→Verner changes the He ionization balance and
hence the metastable He population:

| metric | legacy (ATES 2-term) | Verner (new default) | change |
|---|---|---|---|
| Ṁ (log g/s) | 11.49 | 11.48 | −2 % |
| He 10830 line-center depth | 14.516 % | 14.933 % | **+2.9 %** |
| He 10830 equivalent width | 0.0877 Å | 0.0910 Å | **+3.8 %** |

So unlike the >59 eV cross-section tail (§21, <0.001 % effect), the ground-state
photoionization choice is a **real ~3 % effect on He 10830** — it matters for the
absolute depth and the inferred He abundance. `ATES_photoionization_rate: True`
reproduces the old result exactly (same formula).

**Recombination status — no update needed (EXHALE already current).** The intended
companion update (He recombination → Porter+2012/2013 / Del Zanna & Storey 2022) was
investigated and found to be a **no-op**: those are *emissivity* papers — they update
the He I line cascade / line ratios, **not** the total (effective) recombination to
2³S, and provide no drop-in α(2³S)/α(1¹S). The escape-relevant quantity — total
recombination into the triplet ladder, which all funnels to the metastable 2³S — is
set fundamentally by (case-B total recombination) × (triplet spin fraction 3/4) and
is robust across modern treatments:
- EXHALE (Benjamin+1999): α₃ = 2.10×10⁻¹³ T₄⁻⁰·⁷⁷⁸, α₁ = 1.54×10⁻¹³ T₄⁻⁰·⁴⁸⁶
- 2025 thermosphere paper (arXiv:2509.14499): α₃ = (3/4)(2.72×10⁻¹³ T₄⁻⁰·⁷⁸⁹) =
  2.04×10⁻¹³ T₄⁻⁰·⁷⁸⁹; α₁ = 1.54×10⁻¹³ T₄⁻⁰·⁴⁸⁶ (**identical** singlet)

α₃ agrees to ~2–3 %, α₁ is identical. So EXHALE's He recombination already matches the
latest escape-modeling practice; no change is warranted. (Separately, a trial direct
sum of Cloudy's `he_iso_recomb.dat` per-level radiative recombination gave a total
~2× the case-A value owing to its unmapped 1642-level indexing — confirming that the
raw file is not a usable shortcut.) No recombination switch was added.

**Touched files.** `src/modules/init/parameters.f90` (`ates_photoion_rate` switch),
`src/modules/files_IO/input_read.f90` (parse the optional key),
`src/modules/functions/cross_sec.f90` (`sigma_HeI` branches Verner default / legacy
two-term). Verified: Verner is the default (He I σ at 30/50/100 eV =
5.36/2.02/0.394 Mb), `ATES_photoionization_rate: True` recovers the legacy values
(5.24/1.78/0.347 Mb).

## 23. He 2³S photoionization → two VFKY96 (Verner) forms + power-law bridge (2026-07-01)

**Motivation.** `sigma_HeI23S` was a bespoke 5-segment broken power law
(Norcross-1971 fit + TOPbase-tail extension, §21). Since EXHALE already uses the
Verner+1996 (VFKY96) `sigma_VFKY96` routine for He I and all metals, the metastable
triplet is now expressed in the same form, so all photoionization cross sections share
one analytic kernel.

**Change** (`src/modules/functions/cross_sec.f90`). `sigma_HeI23S` is rebuilt from two
VFKY96 single-shell forms joined by a log-linear (power-law) bridge, used piecewise:
- **wing A** (E ≤ x3 ≈ 34.70 eV, threshold → Cooper minimum):
  `[E_th,E0,σ0,ya,P,yw,y0,y1] = 4.78, 2.645, 20.8, 1.0e12, 3.42, 2.681, 1.956, 2.603`
- **bridge** (34.70 < E < 45.59 eV): straight log–log line joining wing A(x3) to
  wing B(x4) — this is the "power-law fitting part in the middle".
- **wing B** (E ≥ x4 ≈ 45.59 eV, bump → tail):
  `[E_th,E0,σ0,ya,P,yw,y0,y1] = 45.59, 49.68, 1052, 0.04393, 2.941, 1.717, 5.488e-5, 1.118`

Each wing is a ~4% fit to the resonance-averaged cross section; the construction is
**continuous at both x3 and x4**, threshold (E<4.78 eV) and the smooth high-E tail
(~E^−7/2) are handled by the wings themselves (no hard cutoff). The former broken power
law is retained, commented out, in the routine.

**Verification.** Rebuilt EXHALE.x; `sigma_HeI23S` reproduces the former broken power law
to ~4% (12 eV: 1.466 vs 1.465; 30: 0.289 vs 0.289; 45.6: 2.57 vs 2.68 — bump 4% low; 70:
0.855 vs 0.858; 100: 0.294 vs 0.295; 300: 0.0110 vs 0.0110), continuous across the bridge
(34.70 eV: 0.222 → 38: 0.502 → 40: 0.795 → 45.59: 2.570 → 46: 2.521).

**Touched files.** `src/modules/functions/cross_sec.f90` (`sigma_HeI23S`: two VFKY96 wings
+ power-law bridge; former broken PL kept commented out).

## 24. Temperature-dependent Penning ionization of He 2³S (2026-07-01)

**Motivation.** The metastable He(2³S) destruction by neutral hydrogen, He(2³S)+H, was
carried as a *temperature-independent* constant Q31 = 5×10⁻¹⁰ cm³ s⁻¹ — the ATES /
Oklopčić & Hirata (2018) value, itself the Roberge & Dalgarno (1982) *sum* of Penning and
associative ionization (also used by Lampón et al. 2020). Taylor et al. (2025, ApJ 989:68)
replaced it with a Maxwell–Boltzmann-averaged, **temperature-dependent** fit to the cross
sections of Morgner & Niehaus (1979) and Cohen & Lane (1971). Penning ionization is often
the *dominant* 2³S loss near the base (abundant neutral H), so it directly affects the
He 10830 prediction.

**New rate (Taylor et al. 2025, Table 2), cm³ s⁻¹, T in K:**

- T ≤ 4000 K:  `Q31 = 1.9e-9 * (300/T)**0.07`
- T > 4000 K:  `Q31 = 9.1e-9 * (300/T)**0.50`

This is ~3× the old constant across ~2000–9000 K. The new rate is the Penning channel only
(He(2³S)+H → He(1¹S)+H⁺+e⁻); the associative-ionization channel (→ HeH⁺+e⁻) is not
represented separately in EXHALE (no HeH⁺), so the single Q31 keeps carrying the dominant
2³S+H destruction term `−n_HeITR·n_HI·Q31`.

**Unconditional default.** The temperature-dependent rate is used *unconditionally*; the old
5×10⁻¹⁰ constant is no longer selectable. (`ATES_photoionization_rate` still reverts only
the He I (1¹S) photoionization cross section to the legacy ATES 2-term fit; it does not
affect the Penning rate.)

**Implementation.** Q31 became a per-cell array (was a scalar). New `penning_HeI_23S(T,coeff)`
in `Cool_coeff.f90` evaluates the piecewise fit over the whole grid; `HeITR_coeffs` calls it
unconditionally. The two callers (`ionization_equilibrium.f90`, `post_process_adv.f90`) pass
`Q31(j)` per cell into `params`; the `System_HeH_TR*` residuals are unchanged (they already
read Q31 from `params`).

**Verification.** Rebuilt EXHALE.x; ran HD189733b (metals + He23S, step-capped) with the new
rate and, for comparison, with the former constant. Relative to the old constant, He(2³S) is
reduced near the base (ratio ≈ 0.28 at 1.0 Rp, 0.92 at 1.1 Rp; column-integrated over
1–2 Rp drops to 0.59×), and the profiles converge to unity above ~1.3 Rp where H is ionized
and Penning switches off — the expected neutral-H-gated signature.

**Touched files.** `src/modules/radiation/Cool_coeff.f90` (new `penning_HeI_23S`);
`src/modules/radiation/util_ion_eq.f90` (`HeITR_coeffs`: Q31 array + unconditional call);
`src/modules/radiation/ionization_equilibrium.f90` and
`src/modules/post_process/post_process_adv.f90` (Q31 array, per-cell `params` assignment).

## 25. He/H diffusive separation — Phase 1 (2026-07-01, experimental, default OFF)

**Motivation.** EXHALE is single-fluid and re-solves composition by *local* ionization
equilibrium (which conserves the element ratio), so He/H was frozen at the input `HeH` at all
radii — EXHALE could not represent the diffusive He/H separation seen in Taylor et al. (2025)
(8%→2.5%) and Xing et al. (2023). This is the top recommendation of
`docs/methodology_aiolos_taylor_xing.md`.

**Method (Phase 1).** New module `species_diffusion.f90` transports the He element ratio
`f = n_He/n_H` with bulk advection **and** a molecular-diffusion drift of He relative to H,
`Φ = −D n_H(∂f/∂r + fG)`, `G = (m_He−m_H)g/(kT)`, `D = 1.52e18(1/m_H+1/m_He)^½ T^½/n`
(Banks & Kockarts 1973), plus eddy diffusion `K_zz` on the gradient term. Solved as a **fully
conservative** implicit (backward-Euler, tridiagonal) update of `n_He`; the new element split
is written back into `f_sp` conserving mass (`Σ m_s f_sp = 1`, metals frozen to H). Gated on
`He_diffusion` (default `.false.`; `He_Kzz` sets K_zz). Design/derivation:
`docs/design_hehe_diffusion.md`.

**Debugging highlights (see design doc §7).** (a) Explicit schemes were unstable at the cold,
dense base (tiny He scale height) → implicit tridiagonal. (b) A non-conservative
material-advection form drained He → rewritten conservative. (c) **Root-cause bug:** a local
time scale named `t0` (`=R0/v0`) shadowed the *global temperature* `T0` (Fortran is
case-insensitive) in `TK = Tcode*T0`, zeroing `TK` (floored to 1 K) and inflating the
settling coefficient `G ∝ 1/TK` by ~5800×, causing total over-settling. The two constants are
numerically near-equal (`R0/v0 ≈ 2.83e4 s`, `T0 ≈ 2.83e4 K`), which hid it. **Fix:** rename
the local to `tscale`.

**Result (HD 209458b, He23S+metals).** Physical mild separation — (He/H)/HeH ≈ 0.6–0.9 in the
inner thermosphere, falling to ~0.17–0.24 aloft — consistent in magnitude with Taylor (→0.3×)
and Xing. Stable (no NaN); `He_diffusion` OFF is byte-identical.

**Known Phase-1 limitations.** A near-base pile-up from the lagged-`n_H` nonlinearity is held
by a physical limiter `fHe ≤ HeH` (aloft observable unaffected); proper fix is a
self-consistent `n_H` inner iteration. Ambipolar field, thermal diffusion, metal diffusion,
and a Newton-finished quantitative `Ṁ` are Phase-2. Flag stays default OFF (opt-in).

**Touched files.** new `src/modules/functions/species_diffusion.f90`;
`src/modules/init/parameters.f90` (`he_diffusion`, `he_kzz`);
`src/modules/files_IO/input_read.f90` (`He_diffusion`, `He_Kzz` keys);
`EXHALE_main.f90` (call after temperature, before ionization); `Makefile`.

## 26. He/H diffusion — Phase 2 (2026-07-02): ambipolar, thermal diffusion, metals

Extends §25 (design doc §7e). The Phase-1 kernel was factored into a shared subroutine
`solve_1elem` (one conservative implicit advection–diffusion–settling solve for any element
vs the H background), reused by He and each metal.

- **P2b — ambipolar-corrected settling (default ON, `He_ambipolar`).** The polarization
  field lifts ions, so the He-vs-H effective settling mass is `Δm_eff = 3 − 0.5(Z̄_He − Z̄_H)`
  (mean charges from the local ionization state): 3 at the neutral base, 2.5 in the fully
  ionized wind. Modest effect (~2% less depletion aloft), physically correct direction.
- **P2c — thermal diffusion (default `He_alphaT = 0`).** Adds `α_T ∂lnT/∂r` to the settling
  coefficient; off by default, activatable via `He_alphaT`.
- **P2d — per-element metal diffusion (default OFF, `He_metal_diffusion`).** Each trace metal
  element diffuses independently vs n_H with its own mass (`melem_A`), binary D, and ambipolar
  correction; ion stages rescaled to the diffused total. Metals trace → no n_H feedback.
  **Validated (HD 209458b C/N/O):** heavier elements deplete more — at 3 R_p, element/base ≈
  He 0.171, C 0.158, N 0.136, O 0.105 (monotonic in mass; matches Xing et al. 2023).
  Regression: OFF ⇒ He/H unchanged, C/H frozen at 1.0.
- **P2a — self-consistent n_H: attempted, NOT adopted.** A Picard n_He↔n_H iteration diverges
  (the near-base settling genuinely concentrates He, driving n_H→0); the physical cap
  `fHe ≤ HeH` remains the correct limiter. Documented in design doc §7e.

**Touched files.** `species_diffusion.f90` (`solve_1elem` kernel; ambipolar `dmeff`; thermal
term; metal-diffusion loop); `parameters.f90` (`he_ambipolar`, `he_alphaT`,
`he_metal_diffusion`); `input_read.f90` (`He_ambipolar`, `He_alphaT`, `He_metal_diffusion`).

## 27. Metal-diffusion rescale: ratchet bug found in review, fixed (2026-07-02)

A code review of the diffusion changes found two bugs in the §26 metal-stage rescale that
together made metal depletion *irreversible*: (i) an `rX ≤ 1` clamp — redundant for the
metal/H ≤ reservoir cap (already enforced by `min(nX, fXbase·nHl)`) but forbidding any
*replenishment* of a previously depleted cell (a one-way ratchet); and (ii) a skip of cells
where the element was negligible, which made exhausted cells permanent holes. During early
relaxation (wind undeveloped) settling transiently depletes metals; the ratchet locked that
in. **Fix:** rX is no longer clamped above (the cap alone bounds it), and an exhausted cell
that the solve replenishes is re-seeded through the neutral stage (ionization equilibrium
re-partitions next step).

**Retraction:** the previously reported WASP-121b "metal homopause" (Fe → 0 by ~1.25 Rp, with
a ~2500 K hotter thermosphere) was an artifact of this ratchet — with the fix, WASP-121b
metals track He (even Fe is advection-dominated there, w_s/v ~ 1e-3) and the v1/v2
temperature structures agree. HD 209458b results return to the pre-hybrid values (He I 10830
70.5% → 27.1%, −2.6×; mass ordering He > C > N > O aloft unaffected).
`docs/version_compare.{md,tex,pdf}` updated accordingly.

**Touched files.** `species_diffusion.f90` (metal-stage rescale: cap-only bound + neutral-stage
re-seed of exhausted cells).

## 28. Lower-atmosphere connection: Tier 1-3 machinery (2026-07-02, all opt-in)

Survey + proposal + implementation: `docs/lower_atmosphere_coupling.{md,tex,pdf}`.
New folder `src/modules/lower_atmosphere/`; everything default-off (regression: standard
HD 209458 b 3000-step run reproduces v1.0's Mdot 8.97 exactly).

- **Tier 1** `lower_column.f90` + key `Lower column: <R_1bar RJ>`: Koskinen+2022 analytic
  isothermal-Teq hypsometric column with Visscher chemical-equilibrium H2/H/He; reports the
  derived 1-ubar base radius (equilibrium + fully-atomic bracket), base q_H2/q_H/q_He and mu
  next to the input "Planet radius". Gate passed: Model A r0/R1bar=1.343 (paper 1.34),
  q_H2=0.838 (0.84), q_H=0.027 (0.026).
- **Tier 2 foundation** `h3p_cooling.f90` (Miller+2013 Table-5 LTE emission fits + Table-6
  non-LTE factor; anchors reproduced <0.5%; paper archived as
  `references/Miller_2013_JPCA_117_9770.pdf`) and `mol_rates.f90` (Koskinen+2022 Table-1
  rates R1-R23, verified against the PDF). Remaining: coupled System_HeH_mol solver + EOS.
- **Tier 2a** key `Molecular base: True`: EOS-only base correction removing the H2-bound
  particles from ntot_bc via the equilibrium fit (chemistry stays atomic; crude, documented).
- **Tier 3** optional `base.inp` (keys T_base/r_base/HeH_base/Kzz_base; echoed, no-op when
  absent) consumed by input_read + driver `src/utils/run_lower.py` (isothermal or Guillot
  2010 semi-grey T(p)); VULCAN photochemistry documented as the upgrade path.

Physics note from the Tier-1 gate: the chemical-equilibrium base stays strongly molecular
up to T~2000 K, so the atomic base of Teq~1000-2000 K hot Jupiters rests on photochemical
dissociation (Moses 2011; Koskinen 2013a) - quantified motivation for the VULCAN tier.

**Touched files.** new `src/modules/lower_atmosphere/{lower_column,h3p_cooling,mol_rates}.f90`,
`src/utils/run_lower.py`; `parameters.f90` (`lower_col_r1bar`, `molecular_base`);
`input_read.f90` (keys + `read_base_inp` + ntot_bc adjustment); `EXHALE_main.f90` (Tier-1
report); `Makefile`.

## 29. Tier-2 molecular chemistry core: H2/H2+/H3+/HeH+ (2026-07-03, opt-in)

New coupled equilibrium system `System_HeH_mol.f90` (7-8 unknowns: H+, He+, He++,
[He 2^3S], H2, H2+, H3+, HeH+ as element-conserving fractions; hybrd1): atomic rows use
EXHALE's own rate arrays so the molecule-free limit reproduces the atomic systems; the
molecular channels are the Koskinen+2022 Table-1 network (`mol_rates`). H2 photoionization
cross section sigma_H2 (Yan+1998 Eqs. 17-19; verified vs their Table 7; PDF archived) added
to `cross_sec.f90` and wired into the opacity/photoionization/heating integrals
(`PH_heat_HHe`, optional args). H3+ IR cooling (Miller+2013 + non-LTE factor) enters the
`cool` array. EOS is molecule-aware (`calc_ne/calc_ntot/calc_rho`, `get_species_densities`:
each molecule = 1 particle, molecular ions carry electrons, correct masses). f_sp grows
33->37 (H2 H2p H3p HeHp; output schema + Load-IC names extended; columns zero when off).
Key: `Molecular chemistry: True` (requires He; v1 excludes trace metals; pair with
`Molecular base: True` for a consistent ghost-cell base pressure).

Gates: (0) mol-off regression Mdot=8.97 unchanged; (1) HD 209458 b mol-on: fully molecular
base with a sharp H2->H front at r=1.019 Rp and an essentially atomic wind above it
(Mdot 8.91 vs 8.97, He 2^3S peak +1.4%) -- the hot-Jupiter atomic assumption is now a
RESULT, not an input; (2) hot-Uranus-like (0.0457 MJ, Teq=1140 K): front rises to 1.148 Rp,
H3+ active below it. Documented caveats: local equilibrium only (no molecular advection --
Koskinen's high-altitude H2 replenishment is not reproduced), no Lyman-Werner
photodissociation (<=1.4x on Mdot per Koskinen), P4/P5 photo channels and the 4.48 eV
dissociation sink omitted, gamma=5/3 retained, metals excluded (v1).

## 30. Newton-diffusion co-convergence, Riemann p_min floor, OMP CRITICAL removal, fcheck regression, lya_rt lower-BC audit (2026-07-03)

- **Stall-based Newton hand-off.** He_diffusion runs plateau the flux metric just above
  the hand-off threshold (observed: du frozen at 1.084e-2 vs the 1.00e-2 switch for 1e6
  steps, so the JFNK finish never fired). The Newton warm-up now tracks the du plateau
  (stall counter) and hands over to JFNK when du has stalled within 5x the switch.
  **Validated** (HD 209458 b He_diffusion run): the hand-off engaged at the plateau
  (du = 2.2e-2). The subsequent JFNK then hit the KNOWN base-momentum holdout (worst
  cell j=1, momentum row; ||R|| stalling at ~7.7e-3 while the volume-weighted reference
  residual was already 2.2e-4) and fell back to marching as designed -- so the
  co-convergence loop is validated up to its JFNK-success precondition; the base-cell
  viscosity work remains the enabler (see HD209-BC recommendations memo).
- **Newton + He_diffusion co-convergence.** `Solver: Newton` previously froze the diffused
  He/H field at the hand-off state (the JFNK residual has no operator-split diffusion).
  Now an outer iteration (excited-H pattern) alternates the JFNK solve with 500 diffusion
  relaxation steps at the converged wind, until the He/H field drift < 1e-3 (max 5 passes;
  drift printed per pass). Molecular chemistry + He_diffusion combination is refused (v1).
- **Riemann speed-estimate floor** (`speed_estimate_HLLC/ROE`): `p_min` floored at
  1e-30 p_max before `Q = p_max/p_min` -- no-op for healthy states, prevents a
  pathological division (code-review deferred item, now applied).
- **`PH_heat_HHe` OMP CRITICAL removed**: each thread writes only its own j elements of
  the shared arrays, so the critical section serialized the loop for no correctness
  benefit; values unchanged.
- **`regression/run_fcheck.sh`**: the periodic runtime-checked regression recommended by
  the 2026-07-02 review (rebuild with -fcheck=bounds,do,mem, bounded HD 209458 b run,
  fail on any runtime trap, restore the production build).
- **`lya_rt` lower-boundary audit** (vs Huang et al. 2017): our escape-probability
  closure Jbar = S(1-beta) is local and has NO absorbing bottom boundary, whereas Huang's
  Monte Carlo applies a pure-absorber bottom (H2 accidental-resonance true absorption at
  N_H2 ~ 1e14 cm^-2 makes the molecular layer a photon sink). Consequence: in the bottom
  few scale heights our Jbar (hence the n=2 population there) is likely overestimated --
  the downward-loss channel is missing. The Halpha-forming region (1e-4 to ~1 ubar in
  Huang) lies above the base, so the impact on the TPM Balmer spectra is expected to be
  limited; a quantitative check (add a bottom-loss beta channel, or compare against an
  absorbing-bottom MC) is left as the follow-up. AUDIT note recorded, no code change.

## 31. Tier-3 first light: VULCAN -> base.inp (2026-07-03)

VULCAN (public photochemical kinetics; cloned to `../VULCAN`, FastChem compiled, HD 189733 b
SNCHO-2025 network with `use_photo=True`) + the new converter `src/utils/vulcan_to_base.py`
(.vul pickle -> `base.inp`: photochemical q_H2/q_H/q_He at 1 ubar, hypsometric r_base using
VULCAN's own mu(p) and T(p), molecular mixing ratios as comments; NO metal release -- outside
VULCAN's scope). First light (intermediate state of the converging HD 189733 b run):
q_H2 = 0.63, q_H = 0.23 at 1 ubar -- versus the chemical-equilibrium column's q_H = 0.020,
i.e. **photochemistry dissociates ~11x more H than equilibrium**, directly quantifying the
Tier-1 caveat (Moses 2011; Koskinen 2013a). r_base 1.168 vs analytic 1.174 R_J; T_base 863 K
(the Moses11 T(p), cooler than isothermal Teq=1183 K). The run converged
("successfully run to steady-state", 2356 steps, ~6 h wall): the base state is unchanged
from the quoted numbers (the intermediate save was already converged). End-to-end chain
demonstrated: VULCAN -> vulcan_to_base.py -> base.inp -> EXHALE startup override
(T0=863.4 K, R0=1.1685 R_J, He/H=0.0969, K_zz=1e9 echoed).

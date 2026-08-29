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
wind is undefined there and the sonic-point treatment breaks down. Dropping the
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
echoed in `EXHALE_setup.out`.

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
by several arguments for each ion, and threaded the MINPACK coefficients in each cell
through a hand-indexed flat `params(60)`. This step converts that plumbing to a
**species-metadata table + array-indexed (2D) data flow** so that adding an
element becomes "add rows to a table" instead of "thread N new arguments through
M routines." **No physics changes** -- the network is still C/N/O + Mg, and the
result is gated bit-for-bit against §2 (below).

### What changed (Steps A-F)

- **A -- Metadata module** (`init/species_table.f90`, new). One canonical table
  holds metadata for each ion (length `n_mion = 12`): f_sp column (`mion_fsp`), parent
  element (`mion_elem`), stage / charge^2 (`mion_stage`, `mion_z2`),
  photo-ionizable flag and photo-table column (`mion_isphot`, `mion_iphot`),
  threshold (`mion_ethr`), cooling flag (`mion_iscool`); and metadata for each element
  (length `n_melem = 4`): nuclear charge (`melem_Z`), neutral-ion index
  (`melem_i0`), top stage (`melem_top`). The canonical ion order
  (CI..MgIII = `f_sp` columns 7-18) is defined here once. Element indices
  `iel_C/O/N/Mg`.
- **B -- `write_output` 2D.** Metal ion densities pass as `nm(:, 1:n_mion)` and
  are written with an implied-do over the table (`(nm(j,i)*n0, i = 1,n_mion)`),
  not named writes for each ion.
- **C -- `PH_heat_HHe` + column density 2D.** Optical depth, photoheating, and
  the photoionization rates loop over photo-ionizable ions through
  `sigma_tab(:,mion_iphot(i))` and `Nm_col`, returning a 2D rate array
  `P_m(:, 1:n_mion)` (inert top stages stay zero) instead of `P_CI..P_MgII`
  arguments, one for each ion.
- **D -- `eval_cool` metal channels 2D.** Recombination and collisional-ionization
  coefficients return as 2D arrays (`rec_m`, `aion_m`); bremsstrahlung sums
  `mion_z2(i) * GF_elem(:,mion_elem(i)) * nm(:,i)` over ions. No coefficient
  arguments for each ion.
- **E -- Generalized MINPACK system.** `ion_system_HeH_metals` keeps hybrd1's fixed
  `(N_eq, x, fvec, iflag, params)` signature but (i) receives each cell's metal
  coefficients through a **module-level block** in `System_HeH_metals`
  (`set_metal_coeffs` stores `met_ntot/g0/g1/b0/b1/a1/a2`), set by the driver before
  each `hybrd1` call -- safe because the ionization-equilibrium cell loop is
  serial -- and (ii) assembles `fvec` by **looping over elements** (two balance
  equations each, force-zeroing absent elements). `params` now carries only H/He
  (1-11), HeITR (12-18), and the C/N/O charge-transfer rates (40-48); the old
  slots 49-55 for each element are gone. The driver builds metal densities, initial
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
- **Atomic data for each ion and abundance for each element remain manual entries** -- the
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
- `nonlinear_system_solver/System_HeH_metals.f90` -- module-level metal block for each cell
  + element-loop residual (E).
- `radiation/ionization_equilibrium.f90` -- metadata-driven density build, params
  packing (H/He + charge transfer only), initial guesses, solution unpack (E).
- `init/set_energy_vectors.f90` -- min-over-active-metals grid floor (F).

### Status and relation to the Huang plan

Phase 1a is complete: the plumbing for each element in §2 is now table-driven. Adding
the Phase-1b metals should require only metadata rows plus atomic/abundance data,
with **no changes to `System_HeH_metals.f90` or `ionization_equilibrium.f90`**. Next:
**Phase 1b** -- Si, Ca, Na, K, S added together as a uniform-template batch.

Gate scratch run: `/tmp/ates_gate/cno_mg/` vs frozen reference
`mg_validation/cno_mg/output/`.

## 4. Phase 1b: correct C/N/O cross sections and add Si, Ca, Na, K, S

*Added 2026-06-05.*

Two changes are made together in this update, both enabled by the Phase 1a (§3)
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
first. The solver keeps a uniform packing of **two unknowns for each element**
(`N_eq = 3 + 2*n_melem = 21` with all nine elements in the table) and a
`melem_top` flag for each element:

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
| `set_hco_metals` (setter for each cell) | `set_metal_coeffs` |
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
  across the full T range; the X/H ratio of each element is conserved to ~1e-13 for **all nine
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
  rates. The Badnell-swap for iron is deferred to Phase 6. *(2026-08-19: Phase 6
  is closed and the swap will not happen -- the Badnell DR project does not cover
  either iron stage, section 63.2. The choice recorded here is therefore the
  permanent one, not an interim step.)*
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
- **Numerics.** Run `on` converged with zero NaNs; the X/H ratio of each element is conserved to
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
Huang. *(It landed -- Fe II now cools through a CHIANTI multilevel
statistical-equilibrium table and Fe I through NIST f-values plus Van Regemorter;
see the Phase 2 entries below.)* Next: **Phase 1d** -- charge exchange with H for all metals (Huang Table 4),
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

A new module plus minimal wiring; the solver's old C/N/O CT block for each element is
deleted, not extended:

- `radiation/charge_exchange.f90` **(new)** -- the full Table-4 descriptor arrays,
  the `cx_rate` dispatcher, and the generic residual assembly. Public entry points:
  `cx_init` (build the active reaction list from `cx_full`), `cx_set_cell(T)`
  (evaluate rates in each cell), `cx_add_to_fvec` (add CT source terms to the MINPACK
  residual via the `cx_fvidx` element/stage $\to$ row map), and the `cx_full`
  switch. The zero-abundance guard is automatic: each rate is
  $k\,n_{\rm donor}\,n_{\rm acceptor}$, which vanishes when either reactant is absent.
- `nonlinear_system_solver/System_HeH_metals.f90` -- the hard-coded C/N/O + H
  charge-exchange terms are **removed**; a single `call cx_add_to_fvec(...)` now
  feeds CT contributions into the H, He, and every metal-stage residual row.
- `radiation/ionization_equilibrium.f90` -- `call cx_set_cell(T_K(j))` before each
  cell's `hybrd1` solve, storing each cell's rates used inside the residual.
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
Phase 2 validation gate** against Huang et al. (2023) Fig. 10. The atomic data
for each coolant and the CHIANTI effective-cooling tables are documented in
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

### Channel-by-channel cooling diagnostic (this phase)

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
  `eval_cool` (exact split into each channel); new `write_cool_breakdown_eq` writer;
  `mion_fsp` added to the `species_table` use-list.
- `EXHALE_main.f90` — `use utils_ion_eq, only: write_cool_breakdown_eq`; call it
  after the final equilibrium `write_output`.
- Cooling coefficients themselves live in `Cool_coeff.f90` (provenance in
  `Update_EXHALE_early_phase`, Parts II–III); unchanged here.

### Caveats and deferred

- The metal decomposition for each ion is exact in the default branch
  (`use_2lev_cool = .false.`, the WASP-121b setting). In the optional two-level
  fine-structure branch the terms for each ion are the resonance-line approximation
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
cell-by-cell Newton solve: `excited_H_update` fills two global arrays (one entry for each cell) from
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
  `R_star`, `gamma2_bal`, `hpe2_bal`, `jlya_mode`, `jlya_rt_file`; feedback
  arrays over the cells `gph_balmer_HI`, `heat_balmer`; diagnostics `Jlya_arr`,
  `n2s_arr`, `n2p_arr`, `Sproton_arr`, `Hpe_arr`, `Hdx_arr`; and the proton-budget
  capture arrays `gph_ground_HI`, `cion_HI`, `arec_HII`.
- `radiation/excited_hydrogen.f90` **(new)** — `excited_H_update` (n=2
  solve in each cell, ξ factor, Γ₂, J̄_Lyα mode 0/1, source and heating), `n2_populations`,
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
`beta_esc * ne * n_FeI * Lambda_FeI(T)` channel as the other Phase-2 metals.

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
  so the `eval_cool` metal-cooling loop and the breakdown by channel both pick up
  Fe I automatically.

### Validation on WASP-121b

Re-ran the converged full-metals WASP-121b equilibrium (same SED/grid as the
Phase 2 gate) with Fe I cooling on; the channel-by-channel `Cooling_breakdown.txt` then
gives the Fe I share directly. Diagnostics:
`WASP-121b/analyze_FeI_cooling.py` (left panel = the ported `Lambda_FeI(T)`
coefficient read straight from `Cool_coeff.f90` vs the Fe II coronal curve,
Huang Fig. 5 style; right panel = cooling in each channel vs radius, Huang Fig. 10
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
coefficient EXHALE already tabulates, so the existing assembly
`cool_M = beta_esc * n_e * sum_i n_ion(i) * Lambda(i)` carries it unchanged — the
`n_e` in the prefactor times `Lambda_eff` gives `beta_esc * n_FeII *
(cooling per ion)`, i.e. the correct saturated base cooling.

**Why electron-only SE is sufficient at the (mostly neutral) base.** In the
saturated limit the upper levels reach their **Boltzmann (LTE) populations
regardless of which collider thermalizes them**, and the cooling of each ion becomes
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
required because the naive cell-by-cell solve rebuilds a cubic spline for every one of
the ~4300 transitions at every grid cell: instead the descaled Upsilon is
evaluated **once per transition over the whole T-grid** (`upsilon()` accepts an
array T), and the rate matrix is split `M = n_e · Q_coll(T) + A_rad`, so the
ne-loop only reassembles and solves. The full 41×29 table builds in ~3 s and
matches the naive cell-by-cell solve to **1.5e-15** (machine precision).

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
  cell-by-cell post-process solve.
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
  cooling for each Fe II ion. Because the saturated (LTE) populations are
  collider-independent, the electron-only SE solve is exact in that limit; at
  low `n_e` `Lambda_eff` reduces to the coronal rate, so optically-thin
  upper-atmosphere cells are unaffected. The override flows automatically into
  the `use_2lev_cool` branch (which already reads `c_metal(:,26)`) and into the
  channel-by-channel `cool_chan` diagnostic (which reads `c_metal(:,i)`), so the
  Cooling_breakdown dump stays an exact decomposition.
- **Post-process cell-by-cell solve (`T_equation.f90`).** This path computes its own
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
overcooling removed that message is **gone** — the cell-by-cell T solve no longer
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

The `_adv` profiles that `EXHALE_transit.py` reads are built by `post_process_adv.f90`, which
re-solves the **temperature in each cell** at the advection-corrected structure (the
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

### The cell-by-cell sweep in `post_process_adv`

Mode 2 adds one sweep over cells, placed after the post-advection
`calc_ne` and before the photoheating refresh, inside the same `k`-loop that
re-converges `T`, so metals, H/He and `T` are driven to a *joint* fixed point.
Per cell `j` (skipping `nh(j) ≤ 0`):

1. `cx_set_cell(T_K(j))` — charge-exchange rate coefficients at the local `T`.
2. Build the coefficients for each element (canonical order) from the existing
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

**Tiny-denominator artifacts (not physical).** The "max rel change" scan for each element
reports **O +6538%** at r = 1.000 and **N +17.9%** at r = 1.013. These are
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
equilibrium split is the converged one, and `EXHALE_transit.py`'s lines are dominated by H
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
   rate-coefficient** difference, **not** photon trapping. EXHALE already uses the
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
Roche-truncated at ~2 R_p, and EXHALE's base (r=1) sits at **1 μbar**
(`n_0=3.1×10¹²`, T=2358 K); but Huang defines `R_p` at the **4 mbar transit radius**
and places his 1 μbar hydro bottom at **~1.25 R_p** (Fig 9 top-pressure axis /
Fig 12 green "lower atm." band; the earlier "1.46 R_p" note was wrong). The setup
was therefore corrected to **spherical Case A** with `R_0` set to Huang's 1 μbar
physical radius (`R_0 = 1.25 × 1.766 = 2.2075 R_J`, planet mass unchanged → correct
base gravity) and outer radius 15.2 (= 19 R_p,Huang); EXHALE `r` maps to Huang's frame
as `R_Huang = 1.25 r`. The Case-A `Ṁ = 0.077 M_p/Gyr` is within ~1.5× of Huang's
0.052 (Table 3).

### Trend vs Huang Fig 11 (qualitative — exact match is not expected)

Huang's J̄_Lyα is a **plane-parallel Monte Carlo** RT; ours is an analytic
escape-probability formula, so the two cannot agree quantitatively — a meaningful
quantitative comparison is only sensible once **all** Phase-3b features (incl. the
velocity term) are in. Qualitatively (`WASP-121b/lya_caseA_huang11.py`, not in the tree;
the same Fig. 11 case-A comparison is now `WASP-121b/lya_vs_huang11.py`), both show a
**significant inner peak** (EXHALE ~0.6 @1.35, Huang ~0.5 @1.6 R_p) and comparable
magnitude (0.05–0.6). EXHALE has a **dip at ~1.8 R_p** and **rises gently outward**
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

**EXHALE already has the 1D Roche potential** (`grav_field.f90`, the
`spherical_domain=.false.` branch: planet + stellar-tidal + centrifugal), the
L1-truncated domain (`r_max` = Roche lobe), and the tidal momentum source
(`Source.f90`) — Caldiroli's ATES-v2 Roche mode. So the substellar tidal
hydrodynamics is already in place; no new hydro was needed.

### The base gravity for each case (Huang Table 3)

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
0.32). EXHALE's existing Roche potential thus captures the RLOF enhancement — no
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
optional transonic isothermal-wind initial profile, solved from the EXHALE potential via
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

**Conclusion.** Case D's deep RLOF (Ṁ = 1.03) is **not reachable in EXHALE's 1-D
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
post-processor `EXHALE_transit.py` in two stages: **5a** adds the metal resonance lines (spherical,
gated on Case A); **5b** replaces the spherical geometry with the 3-D Roche-equipotential
reconstruction + velocity broadening (gated on Case D). Pure Python post-processing — no
Fortran change.

**Validation targets** — these are the **effective transit radius `R_eff/R_star`**, not an
absorption percent (the expected-results doc labels the table *"Units: R_p/R_star"*; the
"~30%" in its prose is a loose gloss of 0.30). Our code returns the line-center absorption fraction
`h`, which maps to the same quantity via `R_eff/R_star = √((R_p/R_star)² + h)`:

| Line | Case A | Case D | Observed |
| :-- | :-: | :-: | :-: |
| Mg II λ2796 (4 Å) | 0.182 | 0.302 | 0.309 ± 0.036 |
| Ca II K λ3934 | 0.199 | 0.278 | 0.281 ± 0.009 |
| Hα λ6563 | 0.201 | 0.185 | 0.186 ± 0.003 |
| Hβ λ4861 | 0.174 | 0.135 | 0.143 ± 0.005 |
| Na D2 λ5890 | 0.152 | 0.147 | 0.147 ± 0.002 |

### 5a — metal resonance lines (Mg II, Ca II, Na I D), spherical — done

`EXHALE_transit.py` gained Mg II λ2796, Ca II K λ3934, Na I D2 λ5890 as an **isolated, appended block**
(the validated He 10830 / Lyα / Hα / Hβ pipeline above it is untouched). They are resonance
lines whose lower level is the ion ground state; at ~10⁴ K the excited fine-structure levels
are Boltzmann-negligible, so the lower-level density is the ion density itself
(`n_lower ≈ n_ion`). The ion densities live past the H/He block of `Ion_species_adv.txt`, in
`mion_fsp` order (`species_table.f90`): numpy columns **17 = Mg II, 23 = Ca II, 25 = Na I**.
The block reuses the existing spherical chord / Voigt (`wofz`) / trapezoid-τ / disk-average /
instrument-convolution path (`resonance_depth()`), so it is mechanically identical to
the He/Lyα lines.

**Case A result (line-center % / 4 Å-band %), and diagnosis:**

| Line | line-center | 4 Å-band | R_eff reached | Huang Case A |
| :-- | :-: | :-: | :-: | :-: |
| Mg II 2796 | 26.1% | 2.3% | **4.2 R_p** | 0.182 |
| Ca II K | 5.5% | 0.33% | 1.9 R_p | 0.199 |
| Hα | 3.9% | — | 1.6 R_p | 0.201 |
| Na D2 | 0.8% | 0.06% | 0.7 R_p | 0.152 |

Two real gaps surface, both informative: (1) **Mg II already reaches Huang's extent**
(line-center 26% → R_eff 4.2 R_p, vs Huang's clustered ~3.6 R_p) — the metal treatment is
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
the Eggleton Roche-lobe radius. `EXHALE_transit.py` gained a `geometry='triaxial'` switch
(`triaxial_depth()`): the transit LOS runs along +x; the state at each 3-D point is the
substellar state at `r_eff`; and the **LOS velocity in each sector** `v_sub(r_eff)·x/r − Ω·y`
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

**Status:** the full Phase-5 pipeline (metal + Balmer lines, spherical & triaxial geometry,
wind + rotation broadening, triaxial-radii prediction) is implemented and runs end-to-end;
the hydro/Ṁ reproduces Huang (Case A 0.056, Case B 0.337); and the metal transit radii match
Huang once compared in his own bins — Mg II (4 Å) 0.22 vs 0.182, Na 0.175 vs 0.152. Remaining
items are minor: **Ca II ~1.5× high** (line-center; plausibly Ca II→Ca III or the Ca
abundance), the slightly-high Balmer n=2 (TPM-internal J̄_Lyα estimate), the `R_py·R_pz = R_p²`
normalization (methodological), a terminating Case D run, and the `phase5_transmission.ipynb`
packaging. Note: comparisons must use Huang's definition for each line (4 Å bin for the NUV metals,
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
   re-solves the energy/ionization balance in each cell assuming an outflow (upwinds
   from the next-inner cell). In the dense base the non-monotonic metal line-cooling
   curve gives the energy equation a second, spurious *hot* root that MINPACK can
   land on; the upwind coupling cascades one bad cell outward → sawtooth + spike.
   The pre-existing 2× band guard did not catch the ~1.4× spike.

**Fix (option c — eq fallback where the flow is not a clean outflow).** In
`post_process_adv.f90`, wherever `v ≤ 0` (the breathing inflow base) the
post-processor now *skips* the advection correction and keeps the converged `eq`
ionization and temperature. Three guards, all gated on `pp_metal_on` (metals-off is
byte-identical): the no-He and He H/He advection loops, and the cell-by-cell T solve.
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
metal unknowns only `if (thereis_metals .and. .not. thereis_HeITR)`. Modeling HD189733b
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

**Convergence threshold and `Reconstruction scheme:`-controlled two-stage reconstruction (factual).**
- `du_th` is now a runtime variable (was a compile-time `parameter`) with default
  `1.0e-3` — the original ATES-Code-main value. EXHALE had carried `2.0e-2`, at which
  a run flagged "converged" can still show ~2% spatial spread in the supersonic mass flux
  ρvr².
- The single-stage vs. two-stage choice is set by the `Reconstruction scheme:` line.
  `Reconstruction scheme: PLM+WENO3` selects a two-stage run: the input line
  `du_th [PLM,WENO3]: <du_plm> <du_final>` supplies both thresholds (`du_th_plm` and
  `du_th`), and `EXHALE_main.f90` starts in PLM and switches `rec_method` to WENO3 once
  `du < du_plm` (or PLM stalls), then converges at `du < du_final`. This automates the
  manual PLM → `Load IC` + WENO3 workflow recommended in the ATES README.
  `Reconstruction scheme: PLM` or `WENO3` is single-stage and uses only the first
  `du_th` value (the second is ignored).
- `CFL` is now runtime and settable via the input line `CFL: <value>`.

**Other source-term changes (kept, with caveats).** `energy_semi_implicit.f90` now damps
the implicit energy update with `dF/dT = 1 + c·|dC/dT|` (was `max(0, dC/dT)`), adding
damping on the falling cooling branch; it removed one oscillation in testing but did
*not* by itself fix the convergence behavior below. `Cool_coeff.f90` now interpolates
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
oscillation resembles the known base-breathing behavior, and the lower boundary
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
oscillation) is local time-stepping in each cell, since the global `dt` is currently set by
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
AIOLOS-hard-coded (default off), and the scalar `T_equation` for each cell uses
`pp_nm_cell` only in the post-process path (the time-marching cell-by-cell solve
receives n_e via the now-metal-aware `energy_semi_implicit`).

---

## 18. Automatic initial-condition selection (`IC mode: auto`)

*Added 2026-06-13.*

**Motivation.** The code now carries three IC families — cold hydrostatic
(default), transonic isothermal wind (added for the Phase-4 deep-RLOF Case D,
where the cold IC false-converges), and the hot-Parker warm seed (benchmarked:
no speedup) — and the choice was manual. Kubyshkina+2018 build an IC for each planet
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
report (`EXHALE_setup.out`) records the family actually in effect, with a
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
sum of Cloudy's `he_iso_recomb.dat` radiative recombination for each level gave a total
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

**Implementation.** Q31 became an array over the cells (was a scalar). New `penning_HeI_23S(T,coeff)`
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
`src/modules/post_process/post_process_adv.f90` (Q31 array, `params` assignment in each cell).

## 25. He/H diffusive separation — Phase 1 (2026-07-01, experimental, default OFF)

**Motivation.** EXHALE is single-fluid and re-solves composition by *local* ionization
equilibrium (which conserves the element ratio), so He/H was frozen at the input `HeH` at all
radii — EXHALE could not represent the diffusive He/H separation seen in Taylor et al. (2025)
(8%→2.5%) and Xing et al. (2023). This is the top recommendation of
`docs/methodology_aiolos_taylor_xing.md`.

**Method (Phase 1).** New module `species_diffusion.f90` (deleted 2026-08-25, section 69;
replaced by `src/modules/functions/binary_element_diffusion.f90`) transports the He element ratio
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
- **P2d — metal diffusion for each element (default OFF, `He_metal_diffusion`).** Each trace metal
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

## 28. Lower-atmosphere connection: Tier 1-3 (2026-07-02, all opt-in)

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
- **`backup/regression/run_fcheck.sh`**: the periodic runtime-checked regression recommended by
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

## 32. VULCAN as a subroutine-style pre-step (2026-07-06)

VULCAN (with FastChem inside, as shipped) is now included under `EXHALE/VULCAN/` and can be
invoked by EXHALE itself: `Lower atmosphere: vulcan <R_1bar>` in input.inp runs the new
`src/utils/vulcan_driver.py` at startup when no base.inp exists (builds the planet's
Guillot T(p)/Kzz atmosphere, picks a stellar UV spectrum by host Teff, compiles FastChem
once, runs VULCAN to steady state -- hours on first run, cached .vul afterwards -- and
converts to base.inp), then proceeds with the wind solve. `Lower atmosphere: analytic
<R_1bar>` invokes the fast equilibrium column instead; no key = classic base (fully
optional). EXHALE_ROOT env overrides the code root for relocated installs. The download URLs, required citations (Tsai+2017,2021; Stock+2018,2022 for
FastChem; github.com/shami-EEG/VULCAN, github.com/NewStrangeWorlds/FastChem) and the exact
modifications setup_vulcan.sh applies are in README "Obtaining VULCAN and FastChem". Tested:
analytic branch end-to-end, VULCAN branch with cached .vul (photochemical HD189 base
T=863 K applied), and opt-out regression.

## 33. TPM.py renamed to EXHALE_transit.py (2026-07-06)

The transmission-spectrum post-processor `TPM.py` (the acronym clashed with "Trusted Platform
Module" and did not describe the tool) is renamed `EXHALE_transit.py`, matching the
`EXHALE_*` convention. Run-time overrides are now `EXHALE_TRANSIT_<NAME>` (PATH, SAVE_PREFIX,
HE_LMIN/LMAX/N, RSTAR_RSUN, ROTP, TSTAR); the old `TPM_<NAME>` names remain as
backward-compatible fallbacks (a `_tenv()` helper checks the new name first, then the old).
Output filenames are unchanged (`<prefix>tpm_He10830.txt` etc.) so downstream notebooks and
figure scripts keep working. All in-repo references (scripts, run_tpm_all.sh, examples,
docs, README/HOWTO/manual) were updated; VULCAN/ untouched.

## 34. Sub-level-resolved Balmer opacity + LaRT 2s inclusion (2026-07-06)

The H-alpha (and H-beta) transmission was computed as `f_multiplet * n_2` with the
statistical multiplet oscillator strength (Ha `f_23 = 0.6407`). This implicitly assumes the
`2s:2p` populations sit in the `1:3` statistical ratio. They do not: the metastable `2s`
(`A_2s1s = 8.26 s^-1`, two-photon) and the Ly-alpha-pumped `2p` depart strongly from it. The
`2s` and `2p` sub-levels also have **different Balmer cross sections** (only `2s->3p` is
allowed from `2s`; `2p` goes to `3s`+`3d`). Fixed to sum per sub-level,
`tau ~ f_2s n_2s + f_2p n_2p`, with NIST/Wiese absorption oscillator strengths
`Ha: f_2s=0.4349, f_2p=0.70941` (`=2p->3s 0.01361 + 2p->3d 0.69580`);
`Hb: f_2s=0.1028, f_2p=0.125886`. The fine-structure components share the Doppler-dominated
profile, so this is an exact re-weighting of the same line.

- **`EXHALE_transit.py`** (escape-prob path): now carries `data_n2s`/`data_n2p` separately and
  weights Ha/Hb by `f_2s n_2s + f_2p n_2p`. Effect is small (< a few percent) because with the
  in-code Jlya the `2s:2p` ratio is near-statistical in the Ha-forming region.
- **`examples/tpm_halpha_lart2d.py`** (LaRT path): previously used the `2p`-only population
  (`_, n2p, _ = n2_populations(...)`) with `f_23` — it **discarded the `2s` absorption
  entirely**. Now `build_n2p_2d` returns both `n2s` and `n2p` (masked to the physical
  atmosphere, since `2s` is recombination/collisionally fed everywhere, not only where the
  LaRT scattering rate `Pa>0`), and `halpha_transmission` uses the sub-level sum. The LaRT
  Ly-alpha field correctly pumps only `1s->2p` (the pump enters the `S_2p` source term in
  `n2_populations`, not `S_2s`). Effect is large: +21% (HD189733b) to +67% (WASP-121b).

Verified that the **Ly-alpha absorption line itself is unchanged and already correct** — it
is `1s->2p` only (`f_la = 0.41641`, applied to the ground-state `n_HI`); `1s->2s` is a
forbidden two-photon transition and is not part of the line.

**hd209 benchmark input fix:** `benchmarks/hd209/input.inp` was missing `Stellar Lya flux`,
`Stellar Teff`, and `Stellar radius` (a stray `t` line was also removed). Without the Lya
flux the LaRT converter read `F_lya = 0 -> L_lya = 0 -> Palpha = 0 -> n2p ~ 0`, giving a
spurious `Ha = 0.000%`. Added `Stellar Teff: 6065`, `Stellar radius: 1.155`, and
`Stellar Lya flux: 3.5e3 erg/cm2/s` (a quiet-solar-analog value scaled by R*^2 for this
inactive G0; the absolute LaRT depth scales with this input, so it is provisional).

**Impact on the paper (Table `tab:lart`, Sec. `sec:lartcomp`).** With both paths on the same
VULCAN-BC winds and the sub-level opacity, the escape-prob/LaRT Ha ratios are:
HD209458b 0.84/0.43 (2.0), HD189733b 1.23/0.91 (1.4), WASP-121b 4.41/2.33 (1.9),
WASP-52b 5.69/2.20 (2.6). The inline closure thus overestimates Ha only by ~1.4-2.6x, far
less than the earlier ~4-30x (which combined pre-VULCAN winds with the 2p-only LaRT). Paper
abstract, Sec. 3.2 (new `sec:n2opacity`), Sec. `sec:lartcomp`, summary, and
`figs/make_lart_fig.py` (both curves now from `benchmarks/<p>`) updated accordingly.

## 35. Molecular equilibrium seeding fix, molecular-base coupling, Tier-2 gate refresh (2026-07-16)

- **Frozen-atomic-base bug (diagnosed and fixed).** The 8-unknown molecular
  equilibrium (`ion_system_HeH_mol`, hybrd1) is bistable in its initial guess: from a
  zero-molecular warm start it fails to converge at the dense, optically-thick base,
  and the `thereis_mol` branch never checked `info` -- the failed state was silently
  adopted and then frozen indefinitely (x_H2 = 0.0122 at the base, flat across 130k
  steps and independent of T; numerically NOT a root of the coded network -- it would
  require an H2 photoionization rate 15x the unattenuated stellar bound). The true
  dark-base equilibrium is strongly molecular (3-body formation R15 vs thermal
  dissociation R12), which the earlier gate runs had reached only because the
  then-uninitialized molecular `f_sp` columns acted as an accidental nonzero seed;
  the composition-audit zero-init exposed the latent fragility. Fixes: (i) `info` is
  checked and a failed solve retries ONCE from a physically-informed guess seeded by
  the chemical-equilibrium H2 fit `q_h2_equilibrium(p,T)` (still-failing cells keep
  the previous state and are counted in a one-line summary warning); (ii) `set_IC`
  seeds the molecular H2 column from the same fit (H nuclei conserved) instead of
  zero.
- **Molecular base coupled on.** Solving molecular chemistry now implies the
  molecular base particle count (`Molecular base: True` behavior): with an atomic
  `ntot_bc` the base temperature is inflated by the particle-count ratio (~1.8x for
  a fully molecular base), thermally dissociating the H2 layer the chemistry just
  built. A note is printed when the coupling engages.
- **Tier-2 gates re-established (deterministic; inputs pinned at
  `docs/lower_atmosphere_figs/data_g*/input.inp`).** Gate 0 (HD 209458 b atomic,
  3000 steps): log Mdot = 8.97 unchanged. Gate 1 (HD 209458 b molecular): fully
  molecular base (x_H2 = 0.997, base T unchanged at ~1464 K), sharp H2->H front at
  r = 1.019 Rp, log Mdot 8.87 vs 8.97 atomic, He 2^3S peak +1.3%; stable over 12000
  steps with zero solve failures. Gate 2 (hot-Uranus-like: 0.0457 M_J, R_p = 0.49 R_J,
  T_eq = 1140 K, HD209 orbit/spectrum -- input recovered by matching the archived
  domain fingerprint r_max = 4.814118 to 6 decimals): fully molecular base, front at
  r = 1.166 Rp, H3+ active in the molecular layer. Atomic runs byte-identical
  throughout (make check + OMP identity).

**Gate refresh (2026-07-23, current production defaults).** After the 2026-07-17
onward changes (staged secondary ionization, `He_rec_coupling` ON, He-H charge
exchange ON; §36-39) the three gates were re-recorded at a uniform 12000-step
relaxation-snapshot convention (the earlier gate numbers mixed cap conventions).
Key result: the H2->H fronts and molecular base composition are essentially
unchanged (front 1.019 -> 1.020 for HD209, 1.166 -> 1.156 for the hot Uranus;
base x_H2 0.996 / 0.9998), so the dissociation-front physics is robust to the new
defaults. The wind Mdot now reads log 10.58 for both HD209 runs at 12000 steps --
these are relaxation snapshots, not flux-flat converged: under the new defaults
the HD209 atomic gate does not converge (du reaches a fixed point ~2.6e-6 but the
rho*v*r^2 mass-flux spread plateaus near 4%), so all three gates use the same
step cap. He 2^3S is now strongly suppressed in the molecular base below the front
(He 2^3S + H2 Penning destruction, ~1e4x at the front, rising back to the atomic
value above ~1.3 Rp); the earlier "+1.3% peak" was a whole-domain-peak artifact of
the outer-wind transient, not a base signal. Pinned `data_g*/` refreshed;
`lower_atmosphere_coupling.{md,tex,pdf}` and `fig_mol_structure.pdf` updated.

## 36. Default H/He rates -> Badnell+Mao / Voronov, van Hoof free-free Gaunt factor, free-free charge fix, and secondary ionization (2026-07-17)

Six coupled changes to the H/He microphysics. Items 1-2 are gated by the new
`legacy_hhe_rates` flag (default `False` = new rates); items 3-4 (free-free)
apply always; item 5 (secondary ionization) is gated by `use_sec_ion` (default
`True`); item 6 rides along with the heating rework. Full cross-code context
and the benchmark table: `docs/atomic_data_EXHALE_vs_MoCHII.{md,tex}`.

### What changed and why

1. **H/He recombination default -> Badnell (2023) RR minus Mao & Kaastra (2016)
   alpha_1, giving case B; Badnell dielectronic recombination added for He II.**
   The pre-update H/He case-B coefficients were unattributed fits inherited from
   ATES; a numerical check showed them to be exactly the Hui & Gnedin (1997)
   case-B fits. The new default `alpha_B = alpha_A(Badnell) - alpha_1(Mao)`
   matches the construction MoCHII already uses, so the two codes now share the
   same H/He recombination coefficients. Setting `legacy_hhe_rates = True`
   restores the Hui & Gnedin fits exactly.

2. **H/He collisional ionization default -> Voronov (1997, ADNDT 65, 1).** The
   legacy H I / He I rates were Abel et al. (1997) (reproducing Janev et al.
   1987) and He II a Hui & Gnedin (1997) fit, none source-annotated. Voronov
   (1997) is now the default, matching the metal rates (already Voronov in both
   codes) and MoCHII. Reverted together with item 1 by `legacy_hhe_rates`.

3. **Free-free Gaunt factor -> van Hoof et al. (2014) thermal average.** The
   frequency-averaged non-relativistic `<g_ff>(gamma^2)` (gamma^2 = Z_ion^2
   Ry/kT) is tabulated at 161 points in log10(gamma^2) from -6 to 10 in 0.1 dex
   steps (generator `cooling_data/gauntff_thermal_avg.py`), replacing the earlier
   unattributed two-branch logarithmic fit. Applies always, independent of
   `legacy_hhe_rates`.

4. **Free-free charge-weighting correction (physical bug).** The previous code
   weighted He II (net charge +1) by Z^2 = 4 and evaluated the metal Gaunt factor
   at the nuclear atomic number (Z = 8 for O, 26 for Fe). Free-free emission
   scales with the ion *net* charge, not the nuclear number: Z_ion = 1 for H II,
   He II, and singly-ionized metals; Z_ion = 2 for He III and doubly-ionized
   metals. The most visible consequence is that the free-free cooling of a
   He II-dominated layer drops by a factor of four. Applies always.

5. **Secondary ionization (Shull & van Steenberg 1985) in the main loop.** A
   photoelectron ejected with E0 = e_v - E_th above ~40 eV does not thermalize
   immediately: it deposits only `f_heat(x)` of its excess as heat and drives
   secondary H I and He I ionizations, where x is the ionized fraction of the
   H+He nuclei. The SvS85 asymptotic fits are applied species-resolved:
   `f_heat(x) = 0.9971 (1 - (1 - x^0.2663)^1.3163)`,
   `f_ion,HI(x) = 0.3908 (1 - x^0.4092)^1.7592`,
   `f_ion,HeI(x) = 0.0554 (1 - x^0.4614)^1.6660`.
   Below 40 eV the photoelectron thermalizes fully. Gated by `use_sec_ion`
   (default `True`); `use_sec_ion = False` is bit-identical to the legacy
   full-thermalization path. This differs from `wind_ae`, which combines the
   same SvS85 partition with Dere (2007) tables tied to its fixed spectral grid;
   the main loop runs on a runtime SED and applies the SvS85 species-resolved
   fits directly.

6. **He I 2^3S triplet photoheating channel added.** The metastable He I 2^3S
   (threshold 4.8 eV) previously contributed only opacity and its photoionization
   rate but no photoheating. Its photoelectron energy now enters the heating and
   absorbed-energy budgets and the secondary-ionization source consistently with
   every other absorber.

### Files touched

- `src/modules/radiation/Cool_coeff.f90` -- `rr_badnell`, `dr_HeII_badnell`,
  `rr_mao`, the assembled `alphaB_*_new`, `voronov_ci` / `ci_*_new`, and the
  161-point `gff_avg` van Hoof table.
- `src/modules/radiation/util_ion_eq.f90` -- free-free net-charge weighting in
  `eval_cool`; secondary ionization and He 2^3S heating in `PH_heat_H` /
  `PH_heat_HHe`.
- `src/modules/nonlinear_system_solver/T_equation.f90` -- free-free net-charge
  weighting.
- `src/modules/radiation/ionization_equilibrium.f90`,
  `src/modules/post_process/post_process_adv.f90` -- callers wired for the new
  heating/secondary-ionization signatures.
- `src/modules/files_IO/input_read.f90`, `src/modules/init/parameters.f90`,
  `src/modules/files_IO/write_setup_report.f90` -- `legacy_hhe_rates` and
  `use_sec_ion` declared, defaulted, parsed, and echoed to `parse_dump.txt`.
- `cooling_data/gauntff_thermal_avg.py` -- generator for the van Hoof table.
- `docs/input_schema.md` -- keys K14b (`Legacy_HHe_rates`) and K14c
  (`Secondary_ionization`).

### Quantitative check

New default rates (`alpha_B`, collisional ionization `ci`; cm^3 s^-1), with the
legacy value at 1e4 K for comparison. The legacy `alpha_B(H II)` at 1e4 K is the
Hui & Gnedin (1997) fit, 2.592e-13.

| Coefficient | 1e4 K | 1e5 K | 1e6 K | legacy (1e4 K) |
|---|---|---|---|---|
| alpha_B(H II)   | 2.606e-13 | 3.070e-14 | 2.163e-15 | 2.592e-13 |
| alpha_B(He II)  | 2.807e-13 | 5.006e-13 | 1.015e-12 | 2.616e-13 |
| alpha_B(He III) | 1.536e-12 | 2.370e-13 | 2.244e-14 | --- |
| ci(H I)         | 7.448e-16 | 3.963e-9  | 3.103e-8  | 7.247e-16 |

The Fortran reproduces an independent Python recomputation to all printed
digits, and the legacy branch reproduces the pre-update values exactly. The new
`alpha_B(H II)` at 1e4 K differs from the legacy Hui & Gnedin value by ~1%; He II
and He III differ more, mainly because the new He II value carries the Badnell
dielectronic contribution explicitly. The van Hoof thermal Gaunt factor at
T = 1e4 K and Z = 1 is `<g_ff>` = 1.266. A smoke run confirms the expected
factor-of-four drop in free-free cooling across a He II-dominated layer, and
`use_sec_ion = False` reproduces the legacy heating and photoionization rates
bit-identically.

### Effect on existing results

Because the default rates change, results computed with the new default differ
from the pre-update golden and gate outputs. `legacy_hhe_rates = True` recovers
the old H/He ionization balance, with two always-on exceptions: the free-free
charge-weighting correction (item 4) and the van Hoof Gaunt table (item 3), both
physical corrections intentionally not tied to the switch. Secondary ionization
(item 5) is a further default change; `use_sec_ion = False` restores the
full-thermalization rates bit-identically.

## 37. He recombination radiation -> H I ionization coupling (`use_he_rec_coupling`, 2026-07-17)

### What changed and why

EXHALE carries no diffuse radiation field, so pure case B implicitly assumes
every He II -> He I recombination photon is reabsorbed locally by helium
(`alpha_B = alpha_A - alpha_1`, ground capture dropped). In H-dominated gas that
is wrong: part of the He recombination radiation (the >= 24.6 eV ground-capture
continuum, the 584 A resonance line, the 19.8 eV 2^3S line, and the >13.6 eV
part of the He I two-photon continuum) ionizes H instead. The physics and the
choice of recipe were worked out in `docs/QUESTIONS_2026-07-17.md`.

This adds an opt-in coupling using Draine's (2011) on-the-spot `y`/`z`
parametrization. New input key `He_rec_coupling` (flag `use_he_rec_coupling`),
default `False` -> pure case B, bit-identical to the legacy path. When on, the He
recombination photons that ionize H are added as an extra H I photoionization
rate with its photoelectron heating, and the He II recombination coefficient
becomes `alpha_B + y alpha_1` rather than pure case B.

Two local parameters set the split:

- `y = 1 / (1 + R n_HeI/n_HI)` (Draine Eq. 14.16), the fraction of the
  >= 24.6 eV ground-capture continuum absorbed by H, with the cross-section
  ratio `R = sigma_He/sigma_H` at the 24.6 eV He I ground edge. Evaluated from
  the code cross sections, `R = 6.004`.
- `z = 0.67 + 0.29/(1 + n_e/n_crit)` (Draine Sec. 15.5), the density-dependent
  fraction of case-B (excited-state) cascade photons that ionize H, with the
  2^3S critical density `n_crit = 1100 e^{1.2/T4} T4^{0.5} cm^-3`. This
  interpolates between Draine's two limits (~0.96 low density -> ~0.67 high
  density).

### The two modes

**Atomic (case-B) mode** (no explicit 2^3S): `alpha_1` is the Mao & Kaastra
(2016) ground (1s^2) capture and `alpha_B` the active He II case B.

- He II recombination coefficient `alpha_eff = alpha_B + y alpha_1`.
- Extra H I photoionization rate `n_HeII n_e [z alpha_B + y alpha_1] / n_HI`.
- Photoelectron heating with ground-channel energy `E_gnd = 24.6 - 13.6 =
  11.0 eV` and a cascade-averaged `Ee_casc_HeI ~ 6.14 eV`.

**He I 2^3S (triplet, TR) mode** (computes 10830): the decay channels are
explicit, so `z` is not used -- a 2^3S atom destroyed by photoionization or
Penning emits no 19.8 eV photon, so the cascade must be summed at the actual
channel rates. `alpha_1` is the 1^1S channel (`rec_HeII_11S`).

- 1^1S channel coefficient `y alpha_1 + 0.25 alpha_B`: the net ground capture
  plus the singlet-excited capture channel (0.25 alpha_B) that the current
  network omits.
- H-ionizing photon production summed over channels: ground `y alpha_1
  n_HeII n_e`; singlet-excited `0.8521 x 0.25 alpha_B n_HeII n_e` (0.8521 =
  2/3 from the 584 A resonance line plus 1/3 x 0.5564 from the 2^1S
  two-photon continuum above 13.6 eV); 2^3S radiative decay `A31 n(2^3S)`
  (19.8 eV line, always ionizes H); 2^3S collisionally converted to the
  singlets `n_e n(2^3S) (0.5564 q31a + q31b)`.
- Channel deposit energies 11.0 (ground), 5.53 (singlet-excited), 6.2 (19.8 eV
  line), and 2.512/7.6 eV (collisionally converted two-photon / 584 A).

The two-photon constants are integrals of the Drake, Victor & Dalgarno (1969)
spectral shape, not round numbers: the 2^1S term emits a photon *pair*
summing to 20.62 eV with a distribution peaking at half that, so over the
window above the H I edge it delivers 0.5564 ionizing photons per decay
carrying 8.9646 eV, i.e. a mean photon energy of 16.110 eV and 2.512 eV of
photoelectron energy each. Treating the in-band photons as uniformly
distributed instead would put the mean at 17.109 eV and the photoelectron at
3.511 eV, 40% high. `f_2q_HeI` and `Ee_2q_HeI` in `util_ion_eq.f90` carry
these, and the 0.8521 / 5.53 / 6.14 combinations are formed from them so the
three cannot drift apart.

The corrections are applied at the lagged (pre-solve) densities in the main
ionization-equilibrium loop, and re-evaluated at the advection-corrected
densities in the post-processor (rate/coefficient before the advection solve,
heating after). The He II recombination *cooling* is left at `kT alpha_B`; the
mismatch is at most `y alpha_1 kT`, negligible.

### Files touched

- `src/modules/radiation/util_ion_eq.f90` -- new `he_rec_coupling` subroutine
  returning the corrected He II recombination coefficient, the extra H I
  photoionization rate, and the photoelectron heating; all zero when off.
- `src/modules/radiation/Cool_coeff.f90` -- new `alpha1_HeII_mao` (the Mao &
  Kaastra 2016 ground 1s^2 capture, i.e. the `alpha_1` subtracted in
  `alphaB_HeII_new`).
- `src/modules/radiation/ionization_equilibrium.f90`,
  `src/modules/post_process/post_process_adv.f90` -- callers wired to add the
  coupling to `rcheiiB` / `P_HI` / `heat` (main loop) and to the advection ODE
  and `theat` (post-process).
- `src/modules/files_IO/input_read.f90`, `src/modules/init/parameters.f90`,
  `src/modules/files_IO/write_setup_report.f90` -- `use_he_rec_coupling`
  declared, defaulted `.false.`, parsed, and echoed to `parse_dump.txt`.
- `docs/input_schema.md` -- key K14d (`He_rec_coupling`).

### Quantitative check (HD 209458 b)

Protocol: HD 209458 b `input.inp` + `metals.inp`, the converged profiles loaded
as IC, then four runs `{atomic, He 2^3S} x {off, on}` each relaxed for 15000
steps from the same IC. The IC predates the 2026-07-17 rate update, so both
branches share the same partial re-relaxation and **only the on/off difference
is meaningful** (not the absolute values). Across the wind `y ~ 0.65-0.69` and
`R = sigma_He/sigma_H(24.6 eV) = 6.004`. Script:
`docs/plot_he_rec_coupling.py`; figures
`docs/figures/he_rec_coupling_{atomic,heitr}.pdf`.

| Quantity (on/off) | Atomic mode | He I 2^3S mode |
|---|---|---|
| He I | x1.04 (1.2 R_p) -> x1.25 plateau (1.5-3.5 R_p) -> x1.32 (3 R_p) | x1.01-1.16 |
| He II (minimum) | x0.72-0.74 (1.1-1.2 R_p) -> ~1.0 (2 R_p) | x0.885 |
| n(2^3S) | --- (not computed) | x0.89 (1.05), x0.92 (1.2), x0.98-1.01 (1.5-2 R_p) |
| Temperature | <= +1% (outer +3%) | <= +0.3% |
| log Mdot | 9.58 -> 9.65 (+0.07 dex, ~+17%) | 9.57 -> 9.61 (+0.04 dex, ~+10%) |

The atomic mode carries the larger He ionization-structure change (case B
alone, `y = 0`, under-recombines He), and its extra H photoionization raises the
mass-loss rate by ~17%. In the 2^3S mode the 10830-forming region (1.5-2 R_p)
sees n(2^3S) within ~1-2% of the off case, so the He I 10830 line is essentially
unchanged, while the base 2^3S drops ~11% (the TR-mode 1^1S channel coefficient
change feeds directly into `n_HeII`). The temperature response is small in both
modes.

### Effect on existing results

Default off is bit-identical to the pre-update path; existing goldens and gates
are unaffected. Turning `He_rec_coupling = True` on changes the He ionization
structure and the mass-loss rate as tabulated above.

## 38. Staged activation of secondary ionization (2026-07-23)

### The regression

With secondary ionization on by default (section 36, item 5), the high-gravity
WASP-121 b cases (`backup/regression/wasp_full`, `wasp_he23off`) crashed with a
NaN at ~step 587: negative densities first, then a NaN base pressure and
temperature. The coupling itself is physically correct in magnitude -- at the
shielded base the SvS85 secondary channel raises the H photoionization rate by
~35x over the direct rate, which is the expected size for the ~keV
photoelectrons there. The failure is a transient-path problem, not a steady-state
one: the extra base ionization deepens the base cooling valley and amplifies the
base acoustic mode of the marginal cold initial condition until the hydro
diverges. Two facts pin this down: (a) a run with the coupling off converges
cleanly (du < 1e-3, Mdot 13.30); and (b) restarting with the coupling on from
that converged state is benign -- no NaN, du stays <~ 4e-2 during a ~2000-step
adjustment and re-converges. The steady state with the full physics exists and
is stable; only the cold-start path to it is fatal.

### The fix: converge first without the coupling, then switch it on

A runtime state flag `sec_ion_active` now controls whether the SvS85 partition is
applied, separately from the `use_sec_ion` enable flag. The photoionization /
photoheating routines apply the coupling only when `use_sec_ion .and.
sec_ion_active`; when it is off the code path is bit-identical to the
`use_sec_ion = False` (full-thermalization) path. The main loop keeps
`sec_ion_active = False` while the wind first relaxes; the moment that first
relaxation reaches a convergence/stall criterion, the flag flips to `True`, the
convergence logicals and the mass-flux level-stability history are reset, and the
run continues and re-converges with the full physics. Stops are additionally
held for `N_stall` (default 2000) steps after the flip: du reacts only once the
base adjustment wave driven by the new coupling has formed, so without the hold
the run stopped one step after activation (activation at step 7144, spurious
stop at 7145) with a state that had not yet incorporated the physics. The
Newton ("Solver:
Newton") finish and the direct-solve diagnostic modes (which bypass the marching
loop) switch the coupling on before their solve, and the flag is forced on after
the loop so post-processing and final outputs always carry the physics even if
the loop hit the iteration cap before the flip.

### The `Immediate` override

`Secondary_ionization: Immediate` sets `use_sec_ion = True` and applies the
coupling from step 0, restoring the pre-staging behavior. It is intended for A/B
tests only; on the high-gravity cold-start cases it reproduces the ~step-587 NaN
crash, confirming that the staged path is what removes it.

### Files touched

- `src/modules/init/parameters.f90` -- `sec_ion_active`, `sec_ion_immediate`.
- `src/modules/radiation/util_ion_eq.f90` -- `PH_heat_H` / `PH_heat_HHe` gate on
  the local `sec_on = use_sec_ion .and. sec_ion_active`.
- `src/EXHALE_main.f90` -- init, stage flip on first convergence, Newton-finish
  and diagnostic-mode activation, post-loop force-on.
- `src/modules/files_IO/input_read.f90` -- `Immediate` value parsed.

### Effect on existing results (wasp_full validation)

For `Secondary_ionization: False` the run is bit-identical to before
(verified byte-identical on `wasp_full` and `wasp_he23off`). For the default
(staged on) `wasp_full`: activation at step 7144, re-convergence at step 13489
(du = 9.98e-4), and log10 Mdot moves 13.30 -> 13.22. The ~-17% mass-loss rate
is the real secondary-ionization physics: the SvS85 partition diverts most of
the hard-photon photoelectron energy from heat into ionization and excitation
in the weakly ionized base, so less energy drives the wind. `Immediate`
reproduces the old (crashing) behavior on the high-gravity cases.

## 39. He-triplet physics series from the Falorca & Vidotto (2026) review (2026-07-23)

The line-by-line comparison against Falorca & Vidotto (2026, arXiv:2607.18193)
and Garcia Munoz (2025, A&A 698, A199) — full analysis in
`docs/falorca2026_he3d_review.md` — exposed six physics gaps, all implemented
and validated in this series (its validation work also surfaced the
secondary-ionization NaN regression fixed in §38). Physical correctness is the
acceptance criterion throughout; impact sizes below are informational.

1. **Penning ionization products** (`heh_tr_rows`, `System_HeH_mol`,
   `adv_implicit_HeH_TR`): He(2^3S)+H0 -> He(1^1S)+H+ + e- now produces the
   H+ (+ its 6.2 eV of electron heating, `e_th_HeI - e_th_HeTR - e_th_HI`) it
   always removed the triplet for. In the EUV-shielded base (g_HI
   exponentially small) this is the DOMINANT H-ionization channel: n_HII rises
   2-7x at r < 1.03 on wasp_full, while the wind and Mdot are unchanged. The
   associative branch (-> HeH+, ~10%, GM25) is folded in as H+ (documented at
   the code site).
2. **Collisional ionization in the TR systems** (`heh_tr_rows` b terms;
   `ci_HeI23S` in Cool_coeff.f90, the Black-1981-derived
   6.41e-21 sqrt(T) exp(-55338/T)/E_23S, k(1e4 K) = 3.3e-10): the TR
   (Oklopcic-form) rows silently dropped ALL electron-impact ionization while
   eval_cool kept charging its cooling — both restored/consistent now,
   including CI of the metastable itself (a ~1-2% 2^3S loss at 1e4 K, its
   4.8 eV cooling in the coio channel, and its He I/adv-row counterparts).
3. **Explicit-n(2^3S) cooling** (eval_cool + the semi-implicit energy solver):
   the 10830 A collisional-excitation channel 1.16e-20 sqrt(T) exp(-13179/T)
   n_e n_23S (Black form with the EXPLICIT metastable density — Falorca &
   Vidotto warn Black's implicit steady-state triplet is wrong here) plus the
   0.80/1.40 eV q31a/q31b conversion ledger. The q13 (19.8 eV) term is
   deliberately excluded: the Cen-1992 He I excitation-cooling term already
   carries that channel. The terms are passed into
   `solve_energy_semi_implicit` (all three eval_cool calls) so they act on T;
   without that wiring they were diagnostic-only. Measured: 1-5% of the total
   cooling on WASP-121b (8-11% of local Lya at the 2^3S peak), locally up to
   ~28% on HD 209458 b.
4. **`He_rec_coupling` default ON** (§37 has the model): with the coupling off
   the TR singlet recombination used alpha_1 alone — neither case A nor case
   B — and the He recombination photons ionized/heated nothing. Validated
   stable on wasp_full: base n_2^3S -10%, 10830-forming region within ~2%,
   Mdot unchanged; on HD 209458 b (§37) Mdot rises +0.04-0.07 dex.
   `He_rec_coupling: False` restores the legacy path (byte-identical).
5. **He<->H charge exchange in the default set** (`he_h_charge_exchange`,
   default on; Koskinen 2013 rates via Huang 2023 Table 4): the group-B pair
   moved out of the cx_full-gated table assembly into dedicated
   `he_h_cx_fvec/jac/fvec_adv` routines called by EVERY ionization system
   with He (previously only the metals systems could reach it), including the
   advection pair and the analytic Jacobians of the Newton systems. Effect is
   confined to the neutral-H base — He II x0.59 at r = 1.05, -1% by
   r >= 1.2, Mdot unchanged — exactly the He+ + H0 -> He0 + H+ sink where
   n_H0 dominates. `cx_full` still gates groups C/D. The new input keys of
   §36-§39 were also added to the parser's known-keys list (they parsed
   correctly but tripped the unknown-line warning).
6. **He(2^3S)+H2 Penning ionization** (`penning_HeI23S_H2`;
   `System_HeH_mol` rows 4/5/8 + 4.4 eV heating): GM25 Table A.5
   (Cohen & Lane 1977) fitted as 5.3791e-12 T^0.676 exp(-695.21/T)
   (<=0.13% over 500-10^4 K). On the pinned molecular gates the metastable
   drops 39x (HD209 base) and 430x (hot-Uranus, out to r ~ 1.1) exactly where
   x_H2 -> 1, with the upper wind and Mdot unchanged — GM25's dominant deep
   metastable sink, now present for Tier-2 + triplet runs. The `_adv`
   post-process is atomic-only (documented) and needs no counterpart.

Regression state after the series: `Secondary_ionization: False` +
`He_rec_coupling: False` + `He_H_charge_exchange: False` reproduce the
respective legacy paths byte-identically; the wasp goldens and the parse-dump
corpus were re-snapshotted once at the end of the series (defaults: staged
secondary ionization, He recombination coupling on, He<->H charge exchange
on). Converged wasp_full reference: activation at step ~7144, final count
~13486, log10 Mdot = 13.22.

---

## 40. `du` triggers require a descending crossing (2026-08-11)

The `du < du_th` convergence stop and the JFNK hand-off (with the
secondary-ionization flip that shares its `du` test) now fire only after `du`
has been seen at or above their own threshold; they start armed for a
`Load IC? True` run and unarmed for a generated IC. A transonic IC reported
`du(1) = 1.85e-05` on WASP-121 b, which tripped both on step 1 and ran away to
NaN by step 1220; with the guard the same input converges (step 8706,
`log10 Mdot = 13.21`). The PLM -> WENO3 switch is not guarded (the regression
cases legitimately switch on step 1 at `du(1) = 0.37`). No new input keys.
Details: `docs/newton_scaling_and_base_wall.md` §9.

---

## 41. JFNK line search decides on the true residual (2026-08-11)

The line search in `solve_steady_jfnk` accepted steps on the residual with the
WENO3 weights frozen at the previous iterate, compared it against a memory that
mixed that measure with the true one and was refreshed only on acceptance, and
reported `||R||` and the `info = 0` verdict from the frozen residual as well.
*Measured* on the HD 189733 b hand-off state: 7 of 14 accepted steps raised the
true residual, by factors 1.09 to 9.05, and the reference stood 2.4x below the
true merit of the iterate the search was starting from, which is what produced
the 12-consecutive-failure abort. Trials are now evaluated with the weights
recomputed at the trial state, the merit memory holds true merits and is written
at the start of every outer iteration, and the diagonal scaling is formed before
the merit it scales. The frozen weights still define the inner Newton model.
WASP-121 b now reaches `info = 0` (`||R|| = 5.72e-05`, `log10 Mdot = 13.17`,
superseding the `du`-stop 13.21 of §40); `photo_deep_secion_cont` converges in
24 outer iterations instead of 59; no case exhausts its backtracks any more.
No new input keys, no golden changed. Details:
`docs/newton_scaling_and_base_wall.md` §10.

---

## 42. Metal-line cooling at a cold base: line trapping and the coronal-fit validity floor (2026-08-11)

Two physically wrong statements in the metal cooling were corrected. Both bite
only where the gas is cold and dense, i.e. at the base of a planet whose base
falls out of the coronal regime; the wind is unaffected.

**(a) The escape probability was hardwired to 1.** `eval_cool` built a
`beta_esc` from the lowest XUV-band *continuum* opacity across *one cell*, then
overwrote it with `beta_esc = 1.0` and applied the full optically thin
metal-line cooling everywhere. (AIOLOS `chemistry.cpp:1006` instead scales the
same gray depth by an arbitrary `1e8`, which drives `beta -> 0` and switches
metal-line cooling off; EXHALE had replaced that with the opposite limit.)
Neither quantity is a line optical depth, and neither is grid-independent.
Measured on the converged HD 189733 b base, the dominant coolant line
`[O I] 63um` carries a line-center optical depth `tau = 2.7` through the cold
layer, i.e. an escape probability of `0.16`, not 1.

`beta` is now the line-center escape probability of the line itself, from the
*column* between the emitting cell and the top of the domain
(`Cool_coeff.f90`: `kappa_OI63`, `kappa_CII158`, `line_escape_probability`,
`fine_structure_escape`). The shape is the plane-parallel Doppler result used
by Hollenbach & McKee (1979) / de Jong, Boland & Dalgarno (1980),
renormalized by a factor 2 so that `beta(0) = 1` *exactly* (the published form
tends to 1/2 because it counts escape through one face of a slab; here the
other direction is absorbed by the lower atmosphere, so it removes energy from
the modeled gas either way). The two branches are switched at
`tau_c = sqrt(pi) exp(a^2/4) = 6.967`, where they cross, so the switch is
continuous in value. In the optically thin wind `beta -> 1`, reproducing the
previous behavior.

`beta` enters as `A_ul -> beta*A_ul` INSIDE the two-level solution, not as a
factor on its result: in the subcritical limit the cooling is set by the
collisional excitation rate and must be independent of `beta`, and multiplying
outside would be wrong by a factor `beta` there. Trapping is applied to
`[O I] 63um` and `[C II] 158um` only, the two lines for which the code carries
an explicit two-level solution; everything else stays optically thin, which is
right in the wind and is the residual approximation at the base (< 0.3% of the
base cooling there). A thick resonance line in a metal-rich wind (Mg II h&k)
is still treated as thin — recorded here, and measured in §46, which finds
that to be the correct effective treatment.

**(b) The CHIANTI coronal fits were evaluated far below their validity
floor.** Every metal cooling coefficient in `Cool_coeff.f90` is a fit or a
table built over `1e3-1e5 K`. At the 236 K HD 189733 b base, 99.7% of the
`[O I]` rate came from the fit's softest exponential, `exp(-930.111/T)`, whose
930 K corresponds to no `[O I]` ground-term splitting at all (the splittings
are 227.7 K and 326.6 K); C I is the same case (`exp(-2351.38/T)` against
splittings of 23.6 K and 62.4 K). Below the floor the only excitable metal
transitions are the ground-term fine-structure lines, which the explicit
two-level terms already carry, saturation included.

`coronal_excitation_cutoff(T)` now removes the coronal part below 1e3 K as
`exp(-((T_floor/T - 1)/w)^2)`, `w = 0.5`, and is exactly 1 above it. Value and
`dT`-slope are continuous at the floor (the Brent energy solve and the
semi-implicit update differentiate the cooling in T), and every coefficient is
**bit-identical at and above 1e3 K**. The factor is applied to all
CHIANTI-derived coefficients — the analytic C/N/O and Mg/Ca/Na/Fe forms and
the 1-D/2-D tables, whose `log10 T` axis starts exactly at 3.0 and which
otherwise hold their edge value indefinitely below it — and to the coronal
remainder of `cool_OI_ne_func` / `cool_CII_ne_func`, but NOT to their
two-level parts. The legacy AIOLOS branch (`cno_cool 0`) is deliberately left
alone: its constant floors are crude fine-structure stand-ins, not
extrapolated coronal fits. `w = 0.5` is a modeling choice, not a measurement.

**Measured.** HD 189733 b, cell 1: total cooling `3.40e-05 -> 2.68e-08` erg
cm^-3 s^-1 against a heating of `7.0e-06`, so the base net rate changes sign
from `-2.65e-05` (cooling) to `+6.97e-06` (heating). Almost all of that comes
from (b); (a) alone is a factor ~5 on the `[O I]` two-level term. Sweeping the
code's own cooling assembly in T at that cell's frozen densities, against the
heating measured there, the local radiative balance temperature moves from
**168 K to 486 K** — 168.2 K from (a) alone, 485.5 K from (b) alone, so (b)
does essentially all of the work at the base — and the two cooling curves are
identical from 1e3 K up. A
side-by-side marching test does NOT resolve this: on a 5e4-step scale the base
is riding an inflow transient whose adiabatic heating dominates the radiative
term in both binaries. WASP-121 b
(base at 2358 K, `tau([O I] 63um) = 0.039`, `beta = 0.956`) is unchanged:
`log10 Mdot = 13.17` before and after, outer mass flux to `3e-7` relative,
profile max relative difference `8e-5`. `make check`: `mol_base_handoff`
(metals off) stays byte-identical (PASS); the two metal cases move by the
intended physics with `log10 Mdot = 13.22` unchanged in both — `wasp_full`
max relative difference `3.7e-7` (density), outer mass flux to `2e-8`;
`wasp_he23off` `1.6e-3` (density, temperature) in the narrow band at
`r = 1.012-1.014`, outer mass flux to `1.2e-5`. Goldens NOT re-snapshotted. Separately noted while
comparing: both goldens already contain **negative species densities** in that
band (H I, O I, O II, S II; 17 cells in `wasp_full`, 3 in `wasp_he23off`) — a
pre-existing defect, unrelated to this change.

**Not fixed by this.** The model still has no stellar optical/IR absorption and
no thermal background, so nothing sets a radiative floor at `T_eq` for a
shielded layer; and below ~700 K the guard leaves ions without an explicit
two-level term (C I, N, Mg, Ca, Na, Fe) with no cooling at all — for C I that
omits the real `[C I] 609/370um` lines, whose LTE rate at the HD 189733 b base
is ~1e-10 erg cm^-3 s^-1, i.e. 1e-4 of the local heating. Files:
`src/modules/radiation/Cool_coeff.f90`, `src/modules/radiation/util_ion_eq.f90`,
`src/modules/nonlinear_system_solver/T_equation.f90`,
`src/modules/post_process/post_process_adv.f90`. Details:
`docs/hd189_base_checkerboard.md` §10.

---

## 43. Base grid resolution as an input key (`Base grid [dr,cells]`, 2026-08-11)

The resolution of the uniform region of the `Mixed` grid was two hardcoded
locals in `define_grid.f90` (`N_low = 50` cells of `drc = 2.0e-4` R_p). It is
now an input key, because that spacing is what decides whether the base
carries the stationary 2*dr* entropy mode diagnosed in
`docs/hd189_base_checkerboard.md`: the controlling parameter is the number of
cells per base density scale height, `H/dr` with `H = kT/(mu g)`, and a 4x
refinement was measured to remove the mode.

```
Base grid [dr,cells]:  2.0e-4 50     # the default -- the historical grid
Base grid [dr,cells]:  5.0e-5 200    # 4x refinement at the same 0.01 R_p extent
```

The two numbers share one line, as `du_th [PLM,WENO3]` does, because they are
not independent: their product is the radial extent of the uniform region
(0.01 R_p by default). Refining at fixed extent means dividing `dr` and
multiplying `cells` by the same factor; changing only one moves the junction
with the stretched region. The second value may be omitted (the count then
keeps its default). The key applies to `Grid type: Mixed` only; `Uniform` and
`Stretched` build their grids from `r_max` and `N` alone, and
`EXHALE_setup.out` echoes that the key is ignored for them.

`EXHALE_setup.out` also reports the resolution actually achieved,
`Base scale-height resolution: H(T_eq)/dr = 1/(b0*dr_j(1))` (`b0` is the
surface Jeans parameter, `dr_j(1)` the first cell after the Mixed-grid
smoothing), and warns below 10 cells. On the default grid this is 26.9 cells
for HD 189733 b against 102.3 for WASP-121 b — the margin the high-gravity
planet does not have. `Grid type: Stretched` gives 1.7 cells for HD 189733 b,
i.e. it is unusable at the base for this class of problem.

The smoothing pass at the end of `define_grid` started at the literal `50`,
which was `N_low` written out; it now follows `N_low` so that the junction
between the uniform and stretched regions is smoothed wherever it is.

**Two costs, both structural.** The CFL step is set by the smallest cell, so a
4x finer base needs about 4x more steps for the same physical time. And the
total cell count `N = 500` is a compile-time parameter, so the cells given to
the base come out of the stretched region: on the 4.43 R_p HD 189733 b domain
the stretch ratio goes 1.0119 -> 1.0252 and the cell size at 1.5 R_p grows
from 6.1e-3 to 1.3e-2 R_p. Refining the base therefore coarsens the wind, which
is where `Mdot` and the transmission spectrum are formed. Making `N` runtime-
settable would need the static `(1-Ng:N+Ng)` arrays to become allocatable
throughout; it was not done here.

**Default is bit-exact.** `dr_base` is declared with a *default-real* literal
(`2.0e-4`, not `2.0d-4`) because that is what the old local was: the stored
value is the single-precision neighbor of 2e-4, 2.5e-8 relative below it.
Declaring the exact double moves every base cell by that amount; measured on
`wasp_full` it changed the converged density by 5e-10 relative (step count and
`log10 Mdot = 13.22` unchanged), which is not a physics difference but is not
byte-identity either. The literal is kept so that an `input.inp` without the
key reproduces earlier runs exactly. The consequence to know: writing the
default out explicitly (`Base grid [dr,cells]: 2.0e-4 50`) does *not* reproduce
it, because the list-directed read into `real*8` gives the exact double.

**HD 189733 b, what is and is not shown.** Splitting the alternating-amplitude
window separates the ghost-to-cell-1 boundary jump (cells 1-2) from the 2*dr*
mode proper (cells 3-12). Measured on cold-start re-convergences at the same
stage of their transient, against the pre-fix converged reference: the interior
`A(ln rho)` is 9.0e-3 (reference, base at 236 K), 1.6e-3 on the default grid
with the base now at 601 K (the metal-cooling fix of item 42 alone, a factor 6),
and 9e-5 / 5e-5 at `1.0e-4 100` / `5.0e-5 200` (a further factor 17-31, into the
regime where the mode is gone). The
1-12 window stays at 1.3-2.4e-2 in all cases because it is measuring the
boundary jump, which grid refinement does not touch. **These runs had not
converged**, so no `Mdot` is quoted from them and `HD189733b/` was left alone.

Recorded while doing this: `count_max = 1000000` is a compile-time parameter and
the `EXHALE_MAXSTEPS` hook only lowers it, so a 4x-refined base - whose CFL step
is 4x smaller - can reach at most 1/4 of the physical time the default grid
reaches within the cap. The pre-fix HD 189733 b run used its full 1e6 steps, so
4x is not reachable to convergence there without raising `count_max`; 2x is, and
2x is already where the mode goes.

**Checks.** Byte-identity was tested differentially, against the binary built
immediately before this change rather than against `backup/regression/golden/`
(the goldens are stale by the intended metal-cooling physics change of item 42).
All output files of `wasp_full`, `wasp_he23off` and `mol_base_handoff` are
byte-identical including headers, with the same step counts (13488 / 13482 /
12000) and `log10 Mdot` (13.22 / 13.22 / 10.58). WASP-121 b run to completion
with the key absent: `log10 Mdot = 13.17`, JFNK `info = 0`,
`||R|| = 5.7e-5` — unchanged from item 42. Malformed values are rejected in
`define_grid` (a cell count outside `[2, N-10]`, or a uniform region that does
not fit inside `r_max`).

Files: `src/modules/init/parameters.f90`, `src/modules/init/define_grid.f90`,
`src/modules/files_IO/input_read.f90`,
`src/modules/files_IO/write_setup_report.f90`.

---

## 44. Base ghost temperature, runtime step cap, coronal-cutoff width as input keys (2026-08-11)

Three settings that were compile-time constants or an unexamined default became
`input.inp` keys. All three default to the previous behavior, so a file without
them runs exactly as before (`make check`: all three golden cases
byte-identical).

### 44.1 `Base ghost temperature: isothermal | continuous`

The lower ghost cells pin the density to `rho_bc` (the mass reservoir) and, in
the legacy closure, the pressure to `ntot_bc + dp_bc`, i.e. the temperature to
`T0`. **The `T0` pin has no physical backing**: the radiative equilibrium of the
lower atmosphere is outside the model, so nothing in the code determines the
ghost temperature, and `T_eq` is a choice, not a boundary condition. It becomes
actively wrong once the first interior cell settles far from `T0`. With the
metal-line cooling of item 42, cell 1 on HD 189733 b sits at 494 K against
`T0 = 1183 K`: a factor-2.4 contact discontinuity, plus a 2.3x density
inversion, held permanently on the boundary.

`continuous` imposes `dT/dr = 0` at the base face instead. The ghost keeps the
base composition -- `ntot_bc` nuclei and `dp_bc` electrons at `rho_bc`, exactly
the particle count the isothermal pin uses -- and carries the cell-1
temperature:

```
p_ghost = (ntot_bc + dp_bc) * T_1 ,   T_1 = p_1 / (n_tot + n_e)_1 .
```

`(n_tot + n_e)_1` comes from `get_species_densities`, the single point where the
code decides what counts as a particle, and is stored in `n_part_cell1`. The
ghost pressure therefore remains a differentiable function of the interior
pressure (what the JFNK line search needs), while the ionization state it
divides by is lagged exactly like every other composition quantity across a
hydro step. `BC_component_constrho` is the only place the lower ghost is built,
and marching, the JFNK residual and the reconstruction boundary all go through
it, so the closure cannot differ between code paths.

The key sets the same quantity as `Hydrostatic base: True` (the ghost pressure,
which that key extrapolates from the interior gradient instead), so the two are
mutually exclusive: `Hydrostatic base` is tested first and `input_read` warns
when both appear. With `Base BC: pressure` the microbar target is imposed at
`T0` when `n0` is derived, so under a continuous-`T` ghost only the base
*density* stays anchored; `input_read` notes that too. The closure in effect is
echoed in `EXHALE_setup.out`.

**Measured** on HD 189733 b (production configuration restarted from the
converged production state and Newton-finished, in a scratch copy; the planet
folder was not written to). The baseline restart reproduces the stored state
digit for digit, so the comparison is like for like:

| | `isothermal` | `continuous` |
|---|---|---|
| `T_ghost` / `T_1` [K] | 1183.0 / 493.6 | 539.9 / 540.2 |
| `rho_1/rho_ghost` | 2.33 | 0.95 |
| `A(ln rho)` cells 1-12 | 0.0264 | 0.0092 |
| `A(ln rho)` cells 3-12 / 5-14 | 4.41e-3 / 1.62e-3 | 9.28e-3 / 8.95e-4 |
| JFNK | `info = 0`, `\|\|R\|\| = 5.37e-4` | `info = 0`, `\|\|R\|\| = 2.76e-4` |
| `log10 Mdot` | 9.04 | 9.05 |

The boundary jump and the density inversion are gone; the extended alternating
tail is gone as well (beyond cell 7 the amplitude drops 13x to 420x), and what
remains is a stronger disturbance confined to cells 2-4. Three further
restart+Newton cycles reproduce each state digit for digit, so both are exact
fixed points and no further decay is available from cycling. Details, window
sweep and scope: `docs/hd189_base_checkerboard.md` §13. **Not** measured at that
point: a cold start under the new closure (6 h, 1e6 steps), other planets,
transit observables.

**Adopted for production later the same day.** `HD189733b/input.inp` now carries
`Base ghost temperature: continuous`, and the folder was re-converged with it
together with the ionization-root validation of item 45
(`run_20260811_rootfix.log`): `log10 Mdot = 9.05`, JFNK `info = 0`,
`||R|| = 1.724e-04`, 2002 steps. On that converged state
`T_ghost/T_1 = 529.4/529.6 K`, `rho_1/rho_ghost = 0.95`, and `A(ln rho)` over
cells 1-12 is `1.9e-3`, well below the 0.0092 of the isolated A/B above --- the
root validation appears to remove a further part of the base disturbance
(`docs/hd189_base_checkerboard.md` §14). The paper and `python/paper_data.py`
carry 9.05 for this planet, with He 10830 2.485%, H-alpha 0.524% and Ly-alpha
35.9% line-center depths (2.467 / 0.520 / 35.74 before).

### 44.2 `Max steps: <N>`

`count_max` was an `integer, parameter = 1000000`. Item 43 recorded the
consequence: a 4x-refined base, whose CFL step is 4x smaller, cannot reach a
converged state within a cap that the default grid already used in full, and the
`EXHALE_MAXSTEPS` hook only lowers the cap. `count_max` is now a runtime
variable with the same default. `EXHALE_MAXSTEPS` keeps its separate meaning (a
deterministic mid-loop exit used by the regression harness); it does not raise
`count_max`.

### 44.3 `Coronal cutoff width: <w>`

The coronal-excitation guard of item 42 multiplies every CHIANTI-derived coronal
coefficient by `exp(-x^2)` with `x = (T_floor/T - 1)/w` below `T_floor = 1e3 K`.
`w = 0.5` was a hardwired `parameter`. It is a modeling choice, not a measured
quantity, and the temperature the base settles at depends on it at the ~100 K
level, so it is now settable, still defaulting to 0.5 (division by a variable
holding 0.5d0 is the same arithmetic, hence byte-identical). Values `<= 0` are
rejected at parse time. The manual states the sensitivity where the guard is
documented.

### 44.4 Checks

`make check`: `wasp_full`, `wasp_he23off`, `mol_base_handoff` all byte-identical
with the keys absent. WASP-121 b run to completion with the keys absent:
`log10 Mdot = 13.17`, JFNK `info = 0`, `||R|| = 5.686e-05`, 8228 steps, and its
output is byte-identical to the same case run with a binary built from `HEAD`
without these changes. (The stored `WASP-121b/output/` differs slightly from
both -- it predates this comparison, not these keys.) The parse-dump corpus
(`backup/regression/run_parse_corpus.sh`) gains three lines per case and needs a
deliberate `golden` refresh; it was not refreshed here.

Files: `src/modules/init/parameters.f90`, `src/modules/states/Apply_BC.f90`,
`src/modules/functions/composition.f90`,
`src/modules/radiation/Cool_coeff.f90`,
`src/modules/files_IO/input_read.f90`,
`src/modules/files_IO/write_setup_report.f90`.

---

## 45. Ionization-equilibrium roots validated against the physical simplex (2026-08-11)

`Ion_species.txt` carried negative number densities (H II, O II, O III, Fe I)
in a thin band just above the base: 17 cells at r = 1.0106-1.0140 R_p in the
`wasp_full` golden, 15 cells in the WASP-121b production run, 5 in HD189733b.

The equilibrium systems are polynomial in the stage fractions and possess
roots outside the physical simplex. At the step where secondary ionization is
switched on (step 7186 of `wasp_full`) MINPACK `hybrd1` returned `info = 1` on
such a root -- x(H II) = -5.2e-6, the Fe stage fractions summing to 1.2209,
i.e. an Fe I density 22% below zero. Nothing tested the root, and the next
step's warm start re-seeded the solve from the same values, so the state
reproduced itself for the remaining 6302 steps. The heating paid for it twice:
the metal photoheating channel went negative (-6.6e-7 erg cm^-3 s^-1, which is
impossible), and the ionized fraction, clamped into [0,1], sent the Spitzer &
Scott (1985) heating fraction f_heat to zero, switching off every photoheating
channel above the photoelectron threshold. Total heating in the band sat 93%
below its neighbors.

Every atomic cell solve now validates its root -- each stage fraction >= 0,
and each element's ionized stages summing to at most its nucleus total -- and
on failure restarts from physically defined starting points: the uncoupled
ionization balance at the incoming electron density
(`ionization_balance_at_fixed_ne`, admissible by construction), then the
optically thick fully neutral limit. A stored state that is not physical no
longer seeds the next solve, which is what breaks the self-sticking. A cell
where no starting point produces an admissible root is left on the ionization
balance and reported rather than accepted. A first attempt that converges to a
physical root is accepted unchanged, so healthy cells follow exactly the same
solver path as before. `q_abs`, an absorbed energy rate, must now be positive
before it normalizes the heating efficiency.

Gates (single-threaded; a baseline `make check` on the same tree passed all
three cases, so the differences are this change alone). `wasp_full` FAILs by
intent: negative densities 17 cells -> 0, Mdot 13.22 -> 13.22, T/rho/v/p
within 0.89%/0.77%/1.3%/0.13% of the golden away from the band, metal
photoheating positive everywhere (domain minimum -8.9e-7 -> +1.03e-7).
`wasp_he23off` also FAILs: its golden looks clean but the run reports 35 roots
outside the simplex over 13478 steps, previously accepted in silence; Mdot
unchanged, T/rho/p within 0.28%/0.086%/0.039%. `mol_base_handoff` is
byte-identical. WASP-121b (13.17) and HD189733b (9.05) keep their mass-loss
rates and step counts and lose their negative densities entirely. Goldens were
not re-snapshotted.

Files: `src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/radiation/util_ion_eq.f90`,
`src/modules/nonlinear_system_solver/newton_solver.f90`,
`src/EXHALE_main.f90`. Full account: `docs/ionization_root_validation.md`.

---

## 46. Line trapping in the thick metal resonance lines: measured, no change (2026-08-11)

§42 introduced `beta(tau)` for `[O I] 63um` and `[C II] 158um` and recorded the
thick resonance lines of a metal-rich wind (Mg II h&k, Ca II H&K, Na I D, the
Fe II UV multiplets) as an untreated case. They were measured on the converged
WASP-121 b, HD 209458 b and HD 189733 b runs. **No code path changed.**

The lines are thick — Mg II k reaches `tau0 = 7.6e4` at the WASP-121 b base and
stays above `1e2` through the region where Mg II carries a third of the local
cooling, and over 99% of the integrated Mg II / Ca II / Fe II / Mg I cooling
comes from gas with `tau0 > 1`. That is not the test, though. Trapping
lengthens the random walk of a resonance photon; it does not destroy it. In
two-level equilibrium the correction to an optically thin coronal fit is

```
S = beta A_ul / ( beta A_ul + ne q_ul ),   q_ul = 8.629e-6 Ups/(g_u sqrt(T))
```

— not a factor `beta` — so it bites only above `n_crit,eff = beta A_ul/q_ul`.
For these permitted lines `A_ul ~ 1e8 s^-1`, so even at `beta ~ 1e-6` the
escape rate keeps `n_crit,eff` at `1e9-1e15 cm^-3`, against a maximum `ne` of
`4.3e9 cm^-3` in the WASP-121 b run. The forbidden fine-structure lines are the
opposite case (`A_ul ~ 1e-5-1e-3 s^-1`, `n_crit,eff ~ 1e0-1e5 cm^-3`), which is
why they, and only they, need `beta`.

Measured `Delta Lambda` integrated over the wind: **0.0214%** of the total
radiative losses of WASP-121 b, 0.0022% of HD 209458 b, 0.0015% of HD 189733 b.
The largest local effect is 1.65% in the first WASP-121 b cell, whose `T` is
pinned at `T_eq` by the boundary condition; the loss exceeds 0.1% only below
`r = 1.007 R_p`. Per line it is 0.1% of the Mg II cooling, 0.015% of Ca II,
1.4% of Na I (its `A_ul` is the smallest of the set) and 0.02-0.8% of Fe II.
The measurement is a conservative bound: it uses the Doppler-core escape
probability, which understates `beta` by 4-10x for these Voigt profiles, and no
turbulent broadening.

Checked in passing: the Fe II `a6D` fine-structure lines at 25.99/35.35 um, the
one metal transition with an `A_ul` small enough to trap, reach only
`tau0 <= 0.06`, so the `cool_FeII_ne` statistical-equilibrium table needs no
escape probability either. H I Ly-alpha, outside this scope, has
`tau0 = 1.3e8` and a cooling-weighted `S = 0.9994` — two-level trapping is not
what suppresses Ly-alpha cooling, and the destruction channels that could
(photoionization of `H(n=2)`, collisional `2p -> 2s`) were not evaluated.

`beta = 1` for the resonance lines is kept and is now documented at the code
site as a justified effective treatment with its validity condition, in the
`SCOPE` comment of `src/modules/radiation/Cool_coeff.f90`, instead of an
acknowledged omission. Full measurement:
`docs/resonance_line_trapping.md`.

---

## 47. Metal electrons and a validity range for the post-process advection correction (2026-08-12)

Two defects in `post_process_adv`, both confined to the advection-corrected
`*_adv.txt` output (the post-processor is the last call of the run and its
arguments are `intent(in)`, so the converged wind and the equilibrium `*.txt`
files cannot be affected by anything here).

**The electron density inside the advection residuals was wrong.** It counted
only the H and He electrons while the equilibrium residual it corrects counts
the metal electrons as well (`metal_electron_sum`), and the *same* routine
already used the metal-inclusive `n_e` for its heating, cooling and temperature
solve. In a shielded base the metals are the dominant electron donors --
measured `n_e,metal/n_e = 0.48-1.00` in the base cells of all four paper
planets -- so the recombination terms were low by up to six orders of
magnitude. Solved to machine precision with Brent, the residual *as
coded* had its fixed point at `x_HII = 4.74e-8` in the first HD 209458 b cell
against an equilibrium `7.46e-14`. A field `adv_cell%xe_metal` now carries the
metal electrons per H nucleus into all three advection residuals, with the
definition `calc_ne` and `metal_electron_sum` use. It is not conditioned on
`eos_include_metals`: that switch decides whether metals enter the mass and
particle budget, whereas recombination needs the true electron density. Metals
off gives zero identically.

**The correction was applied where it cannot be computed and is not needed.**
The residuals carry the *neutral* fraction and report the ion density as
`(1-x_HI) n_h`, so the ion fraction inherits the solver's absolute resolution
on `x_HI`, `xtol = sqrt(eps) = 1.5e-8`. Where the equilibrium ion fraction is
`1e-13` the returned value is quantized at `1e-8` with an arbitrary sign: the
first HD 209458 b cell got `x_HII = -4.96e-8` where the exact root is
`+4.74e-8`, right magnitude and wrong sign, and the upwind cascade carried the
negative density outward over 14 cells. Those same cells are in local
ionization equilibrium to `Da = (dr/v)(P_HI + alpha_HII n_e) = 350-2900`, so
the equilibrium solution is the answer the ODE would give anyway. The
correction is now skipped, and the converged equilibrium ionization kept,
wherever the cell is (i) inflowing (`v <= 0` on either face, the upwind
discretization has no upstream cell), (ii) equilibrium-dominated (`Da > 100`)
or (iii) unrepresentable (`x_HII,eq < 1e-6`, less than 1% relative accuracy at
that resolution). (i) and (ii) are statements about the flow, (iii) about the
representation of the unknown, so the pre-existing `pp_metal_on` gate on the
inflow condition is removed -- it made a discretization property depend on the
metal-cooling switch. The temperature loop keeps condition (i) alone, likewise
ungated; a thermal Damkohler number is not evaluated. Clipping the extracted
densities at zero was tested and rejected: it breaks the H nucleus budget
`n_HI + n_HII = n_h`.

Gates. `make check`: PASS, all three cases byte-identical (the goldens compare
`Hydro_ioniz.txt` and `Ion_species.txt`, which this does not write). Isolated
A/B, the pristine HEAD binary and the modified one run from the same
configuration: `Hydro_ioniz.txt`, `Ion_species.txt`, `Cooling_breakdown.txt`
and `Excited_H.txt` byte-identical, only the `_adv` files differ.

| planet | Mdot [log g/s] | negative `_adv` entries | cells held at eq (newly) | outermost held | max change, `r > 1.05` |
|---|---|---|---|---|---|
| HD 209458 b | 9.31 -> 9.31 | 48 -> 0 | 103 -> 124 (21) | 1.0316 | 1.4% (He III) |
| HD 189733 b | 9.05 -> 9.05 | 0 -> 0 | 135 -> 212 (77) | 1.0380 | 0.37% (He III) |
| WASP-121 b | 13.17 -> 13.17 | 0 -> 0 | 0 -> 0 (0) | - | 3.8% (He III) |
| WASP-52 b | 11.63 -> 11.63 | 0 -> 0 | 3 -> 57 (54) | 1.0106 | 3.3% (He III) |

Every cell whose equilibrium ion fraction is below `1e-6` is now held at the
equilibrium value, and no newly held cell lies above `r = 1.038 R_p`. The
change in the wind is entirely the metal electrons: a build carrying the guard
alone removes all 48 HD 209458 b negative entries by itself and changes the
profiles above `r = 1.05 R_p` by exactly zero on three of the four planets and
by `4.2e-5` on WASP-52 b. Transmission spectra move by at
most 0.07% in peak depth (He 10830 of HD 209458 b) and 0.08% in equivalent
width (H-alpha of WASP-121 b).

The re-run reproduces the stored `output/` of WASP-121 b and WASP-52 b
byte-identically but not that of HD 209458 b and HD 189733 b, which restart
from `output/*_IC.txt` files that had been refreshed to the previous run's
converged state; re-running continues that convergence (HD 209458 b `du`
1.199e-2 -> 3.13e-3) and moves the profiles by up to 9% in density at unchanged
`Mdot`. That is unrelated to this change --- the pristine binary follows the same
new trajectory --- and is why the table above is an A/B between two binaries.

Files: `src/modules/post_process/post_process_adv.f90`,
`src/modules/nonlinear_system_solver/ion_cell_state.f90`,
`src/modules/nonlinear_system_solver/System_implicit_adv_{H,HeH,HeH_TR}.f90`.
Full account: `docs/postprocess_advection_validity.md`.

Found and **not** fixed, marked at the code site: `System_implicit_adv_HeH`
(helium on, He 2^3S off) writes its electron-impact ionization terms without
the `n_e n_h` factor that `adv_implicit_H`, `adv_implicit_HeH_TR` and the
equilibrium rows all carry, which also makes them dimensionally inconsistent
with the photoionization rates they are added to. None of the four paper
planets uses that path. (Since fixed, 2026-08-12: the terms now carry
`ionhi*xe*n_h`; commit b0d44bc, TO_BE_DONE item (B).)

---

## 48. Ground-term fine-structure statistical equilibrium for C I, C II, N II and O I (2026-08-12)

The CHIANTI C/N/O cooling fits carry the ground-term fine-structure (FS) lines
in the optically thin, *low-density* limit. Their critical densities are of
order `1e0-1e5 cm^-3`, decades below the base density of an irradiated
atmosphere (`n_e ~ 1e9`, `n_HI ~ 1e14 cm^-3`), so at the base the coronal form
overstates what those lines can radiate by 2.9-8.3 decades. Two of them
(`[O I] 63 um`, `[C II] 158 um`) already had a two-level saturation; C I, N II
and the `[O I] 145.5/44 um` channels did not, and C I was in consequence the
dominant coolant of the whole base layer of HD 189733 b at a value comparable
to the local heating rate.

**What replaces it.** Every C/N/O coolant with a split ground term now has that
term solved in exact statistical equilibrium at the local `(n_e, n_HI)`:

```
solve  sum_j f_j R_ji = f_i sum_j R_ij,   sum_i f_i = 1
R_ul = C_ul + beta_ul A_ul,   R_lu = C_ul (g_u/g_l) exp(-E_ul/kT)
C_ul = n_e k_e,ul(T) + n_HI k_H,ul(T)
W_FS = sum_{u>l} f_u beta_ul A_ul k_B E_ul        [erg/s per ion]
```

and the coronal curve is refitted to keep only the channels that *leave* the
ground term (`Lambda_rem`). The coefficient handed to the assembly is
`Lambda_eff = W_FS/n_e + Lambda_rem`, so the `n_e n_ion` prefactor recovers the
H-collision channel exactly. Treated: C I `2p2 3P` (609.1/370.4 um), C II
`2p 2P` (157.7 um), N II `2p2 3P` (205.3/121.8 um), O I `2p4 3P`
(63.2/145.5/44.1 um). That is the complete set -- N I and O II have a
single-level `4S` ground term, Mg I/II, Ca II and Na I a single ground level,
Fe I is built from permitted lines only, Fe II already uses a 2-D
statistical-equilibrium table.

Limits: `n_e, n_HI -> 0` collapses onto the ground level and reproduces the
coronal within-term channel to the accuracy of the collision-strength fits
(1.3-4.4%); high density saturates at the exact multilevel LTE emission
(verified to `1e-14` relative). `beta_ul` enters inside the solution, not as a
factor on the result, which is the only placement that is right in both limits.
All eight lines with a transition probability now get their own line-center
escape probability from the outward column; previously only two did.

**Two defects fixed on the way.** (i) The two-level form multiplied a solution
already normalized to the two-level partition function by the ground-term
Boltzmann fraction of its lower level, so its LTE limit was low by
`1 + (g_u/g_l) exp(-E/kT)` -- 2.7x for `[C II] 158 um` at the 530 K base of
HD 189733 b. (ii) The H-atom de-excitation rates (a constant `4.0e-11` for
C II, `4.2e-11 (T/100)^0.67` for O I, both flagged approximate in the source)
are replaced by fits to published quantum scattering calculations: Yan & Babb
(2023, MNRAS 518, 6004) for C I and N II, Abrahamsson, Krems & Dalgarno (2007,
ApJ 654, 1171) for O I, Barinovs et al. (2005, ApJ 620, 537) for C II. They are
5-20x larger.

**What it does.** The `coronal_excitation_cutoff` guard was standing in for the
missing saturation, so its width `w` behaved as a tuning knob on the base
temperature. It no longer does: the local balance temperature of the base cell
is identical to every digit printed over `w = 0.02-1.2` on all four planet
runs, against a 500-650 K spread before. Converged runs:

| run | base T | Mdot at 2 R_p |
|---|---|---|
| HD 189733 b (warm restart, `w = 0.1`) | 528.5 -> 613.2 K | 2.351e9 -> 2.850e9 g/s (+21.2%) |
| WASP-121 b (cold start) | 2406.9 -> 2438.8 K | 3.106e13 -> 3.370e13 g/s (+8.5%) |
| `wasp_full` golden case | 2340.6 -> 2357.3 K | 3.482e13 -> 3.549e13 g/s (+2.0%) |

Transit depths on HD 189733 b move by +4.4% (Ca II) to +20.7% (H-beta). The
profiles outside `r ~ 1.1 R_p` move by a few percent; the change is a base-layer
change that propagates through the density.

Files: `src/modules/radiation/Cool_coeff.f90` (the SE solvers, the four ion
coefficients, the generalized line opacity and escape routine),
`src/modules/radiation/util_ion_eq.f90`,
`src/modules/nonlinear_system_solver/T_equation.f90`,
`src/modules/post_process/post_process_adv.f90`,
`cooling_data/fit_fs_saturation.py` (generator; prints every coefficient),
`cooling_data/make_cooling_doc_figures.py`, `docs/cooling_formulas.tex`,
`docs/coronal_cutoff_width.md` section 7, `TO_BE_DONE.md` item (C).

Also fixed here, unrelated but found while reading the cooling diagnostics:
`examples/exhale_io.py` named only four fixed channels in
`Cooling_breakdown.txt` while the writer emits six (the collisional-excitation
channel is written split by absorber), so every metal column was labeled two
ions too early. `python/paper_data.py` parses the header and was already
correct.

Still open: with the artifact gone it is plain that the base is not in local
radiative balance at all -- at the HD 189733 b base the total radiative cooling
is `3.5e-8` against a heating of `3.4e-6`. Its temperature is set by the inner
boundary and by the flow.

---

## 49. Base ghost built from the composition of the state it bounds

*Added 2026-08-12.*

`Apply_BC` closes the `Base ghost temperature: continuous` ghost with
`p_ghost = (ntot_bc + dp_bc) * p(1) / n_part_cell1`, where `n_part_cell1` is the
cell-1 particle count written by `get_species_densities`. Two call sites filled
the ghosts before any composition solve had run on the state in question, so
they used whatever `n_part_cell1` was left over -- the `input_read` placeholder
`ntot_bc + dp_bc` in the first case:

- `init.f90`: the initial `Apply_BC`. On HD 189733 b the placeholder differs
  from the true cell-1 particle count by 4.4%, and the finite-volume steady
  residual of the loaded state (which stops right after `init`) reported the
  resulting contact discontinuity at the base face as a cells 1-2 spike of 60
  (mass) and 120 (energy) per sound crossing time.
- `steady_newton.f90`, `eval_residual`: `Apply_BC` after `unpack_U`. There
  `F(Y)` depended on the previous `Y` as well as on `Y`, so a finite-difference
  Jacobian column mixed two states.

Both now evaluate the composition of the interior first. `init.f90` also seeds
`base_flux_const` from the initial state, for the same reason: with
`Base velocity: massflux` that constant starts at `-1`, which `Apply_BC` reads as
"not available yet" before falling back to the valve, so the first step and any
diagnostic stopping right after `init` used a different velocity boundary
condition than the run. Files: `src/modules/init/init.f90`,
`src/modules/time_step/steady_newton.f90`.

Checks. `n_part_cell1` is read only by the continuous branch and
`base_flux_const` only by the mass-flux branch, so a run setting neither key is
unaffected and `make check` is byte-identical; an `isothermal` HD 189733 b
residual reproduces bit for bit. With the key on, the
cells 1-2 mass residual falls from `5.97e+01, 5.92e+01` to `6.41e-02, 7.70e-02`
and the energy residual from `1.12e+02, 1.20e+02` to `1.31e-01, 1.54e-01`, with
cell 3 outward and the wind window unchanged to every digit. Re-converging the
planet from its stored state: JFNK scaled merit `||Fs||_2` `1.43 -> 0.315`,
`info = 0`, `log10 Mdot` `9.15 -> 9.14`. Measurements and term-by-term split:
`docs/hd189_base_checkerboard.md` section 15.

## 50. H(n=2) rate coefficients, the trapped-2p balance, and de-excitation heating (2026-08-12)

*Added 2026-08-12.*

Four corrections to the excited-hydrogen and Ly-alpha modules, plus the
structural change that made the first one possible. Measurement and the full
number tables: `docs/lya_destruction_channels.md` section 11.

**One definition of the n=2 rate coefficients.** The Ly-alpha and n=2 atomic
data and every collisional rate coefficient existed in three copies -- inside
`n2_populations`, again as standalone functions in the same file, and a third
time in `lya_rt.f90`. They are now one module,
`src/modules/radiation/hydrogen_n2_rates.f90` (`c1s2s_rate`, `c1s2p_rate`,
`c2s2p_rate`, `c2s1s_rate`, `c2p1s_rate`, `c2p2s_rate`, `alpha_B_hydrogen`,
`alpha_2s_hydrogen`, `alpha_2p_hydrogen`, `n2p_destruction_rate`, and the line
constants). `excited_hydrogen.f90` and `lya_rt.f90` use it.

1. **Statistical weights shadowed by dummy arguments.** `n2_populations` was
   declared `(T, n1s, ne_l, Jlya, G2s, G2p, n2s, n2p)`. Fortran does not
   distinguish `G2s` from the module parameter `g2s`, so the three
   detailed-balance rate coefficients written out in the routine body evaluated
   their weights as the Balmer photoionization rate the caller passed:
   `2s->1s` and `2p->1s` de-excitation came out 5 to 306 times too small
   (planet-dependent, scaling with the stellar Balmer continuum) and the
   `2p->2s` l-mixing rate 3 times too large. The dummies are now `gam_ion_2s` /
   `gam_ion_2p` and the body calls the shared rate functions, so the weights are
   resolved at module scope and out of a dummy's reach. Two further shadows of
   the same kind, both harmless, were renamed: the frequency-sample count
   `ng = 400` in `gamma_n2_balmer` and `heat_n2_balmer` shadowed the global
   ghost-cell count `Ng`; it is `n_nu`.

2. **Recombination cascade source.** `a2s*ne**2` / `a2p*ne**2` -> `a2s*ne*nHII` /
   `a2p*ne*nHII`. The electron recombines onto a proton; at the base the
   electrons come largely from helium and metals while hydrogen is still
   neutral, and `ne/nHII` reaches 2 on WASP-121 b and 5e6 on HD 209458 b. The
   Python transmission tool (`exhale_transit_lib.py`, `EXHALE_transit.py`) and
   `examples/tpm_halpha_lart2d.py` carried the same expression and were
   corrected the same way; their `n2_populations` signatures gained an `nHII`
   argument.

3. **Destruction channels in the `J_int` closure.** `lya_rt.f90` built the
   trapped field from `n2p = P/(A_2p1s beta)`. The 2p budget also loses atoms to
   collisional de-excitation, n=2 photoionization, and l-mixing followed by
   two-photon decay; the denominator is now `A_2p1s beta + D` with `D` from
   `n2p_destruction_rate`. Where `beta` is smallest the omission over-estimated
   `J_int` by up to 1.5x. In the same module the top-down Ly-alpha line-center
   optical depth, previously two identical inline loops, is the single routine
   `lya_line_center_optical_depth`; all three `jlya_mode` values call it, so the
   `tau_Lya` column of `Excited_H.txt` is filled in the LaRT import mode instead
   of being written as zero.

4. **`Deexc heat` now defaults on.** The comment justifying the old default --
   that collisional de-excitation heating overlaps the H I collisional-excitation
   cooling -- was wrong: `coex_rate_HI` is the one-way Cen (1992) coronal rate
   with no de-excitation term in it, so `Hdx` is the correction to that limit,
   not a double count. Most of the n=2 population it de-excites was pumped by
   Ly-alpha rather than excited by a collision, which makes the term absorbed
   Ly-alpha thermalized by a collision, the same size as the photoelectric
   heating the code already applied unconditionally. The key now also honours
   `Deexc heat: False`, which it previously ignored. **This changes results for
   any run that supplies a stellar `T_eff`/radius and does not set the key.**

Audited and left alone: the energy equation applies no escape-probability factor
to the Ly-alpha cooling. `cool` is built from the bare Cen coefficient
(`Cool_coeff.f90:843` -> `util_ion_eq.f90:640,787`); `beta` appears nowhere in
that path, in any `jlya_mode`. That is the correct form -- 98.2-99.99% of trapped
Ly-alpha photons escape and the collisional destruction that would suppress the
cooling is at most 0.06% -- so the historical expectation of a ~10x transfer
suppression is not in the code and should not be.

**Checks.** `mol_base_handoff` (no stellar `T_eff`, so the excited-H model is
off) is byte-identical at every stage, as is its comparison against the stored
golden. `wasp_full` and `wasp_he23off` move: `log10 Mdot` 13.23039 -> 13.22966
and 13.23106 -> 13.23026, with the largest local temperature change 0.57%
(item 1), 0.081% (item 2) and 0.147% (item 3); both cases pin `Deexc heat: True`
in their `input.inp`, so item 4 moves nothing in them. Their goldens were
already stale against the pre-fix tree and were left as they are. Re-converged
from their stored initial conditions in a scratch copy, WASP-121 b gives
13.2075 -> 13.2030, HD 189733 b 9.1286 -> 9.1392 and HD 209458 b
9.3790 -> 9.4581; the H-alpha and H-beta transit depths move by 0.4-1.7%
(`n(2p)`, which carries the opacity, barely moves, while `n(2s)` halves).
HD 209458 b is the largest move (0.075 dex, from items 1-3) because its
Ly-alpha cooling is a negligible share of its radiative losses while the n=2
heating terms are not; it is also the case `docs/lya_deexcitation_heating.md`
recorded as going NaN at the base with the de-excitation heating on, which no
longer reproduces with the corrected populations (re-converges from its stored
initial condition, and a cold start survives a 20000-step cap). A build at
`-O1 -fopenmp -g -fcheck=bounds,do,mem` running the `wasp_full` configuration is
clean. Files: `src/modules/radiation/hydrogen_n2_rates.f90` (new),
`src/modules/radiation/excited_hydrogen.f90`,
`src/modules/radiation/lya_rt.f90`, `src/modules/init/parameters.f90`,
`src/modules/files_IO/input_read.f90`, `Makefile`, `exhale_transit_lib.py`,
`EXHALE_transit.py`, `examples/tpm_halpha_lart2d.py`.

---

## 51. Molecular chemistry and trace metals solved in one system (2026-08-13)

The parser refused `Molecular chemistry: True` together with a `metals.inp`.
That refusal blocked two things: the deep (1e-3 bar) lower-atmosphere handoff
of `base_composition_handoff_plan.md` §11.7, whose A/B pair carries metals and
whose consistent configuration is Tier-2 molecular chemistry; and trace-metal
work on sub-Neptunes, which have a molecular base by construction.

The two networks cannot be solved apart, and the reason is quantitative rather
than architectural: **in the shielded molecular base the metals supply
essentially all the free electrons.** Measured on the two smoke tests below
(base cell, molecular-base planets, solar metals):

| case | n_e with metals [cm^-3] | metal share | n_e metals-off [cm^-3] |
|---|---|---|---|
| hot Uranus, solar 7 metals | 1.20e8 | 1.000 | 6.29e4 |
| HD 209458 b, 1e-3 bar handoff, C/N/O | 9.42e5 | 0.990 | 9.27e4 |

The reason is that the flux which ionizes the low-IP metals -- Mg I at
7.65 eV, Fe I at 7.90 eV, C I at 11.26 eV -- is below both the H I edge
(13.60 eV) and the H2 threshold (15.4 eV), so it is absorbed by neither and
reaches the base almost unattenuated. On the hot-Uranus case it holds Mg and Fe
about 1.6% ionized at the base, which is already enough to swamp the H/He and
molecular electrons. Every recombination and electron-impact term of the
molecular network scales with n_e, so solving the network at the metal-free n_e
misses a factor 1900 (hot Uranus) to 10 (HD 209458 b at 1e-3 bar) in that
layer.

**The merged system.** `System_HeH_mol_metals` follows the pattern
`System_HeH_TR_metals` established: the molecular unknowns keep their rows, the
metals are appended above them, and both blocks are shared code rather than
re-derived. The eight balance rows of the molecular network moved into
`mol_heh_rows` (still in `System_HeH_mol`, with the Koskinen 2022 rate
coefficients, which stay defined in exactly one module) and take n_e as an
input, exactly as `heh_tr_rows` does; the metal rows, the metal electron sum
and the stage fractions come from `ion_residual_core` unchanged. Layout:

```
x(1..3)   H II, He II, He III        x(4..7)  H2, H2+, H3+, HeH+
x(8)      He 2^3S      (when tracked)
x(mbase..)  X+, X++ for each metal element,  mbase = metal_row_base() = 8 or 9
N_eq = 7 (+1 with the triplet) + 2*n_melem
```

`metal_row_base()` is the single definition of that offset, read by the
residual, by the driver's seeding and root validation, and by `cx_metal_base`
so the Huang Table-4 charge exchange writes to the shifted metal rows. With
`met_nelem = 0` the system reduces exactly to `System_HeH_mol`.

**Root validation and seeding (§45 extended to the molecular branch).** The
molecular cell solve is bistable because the cell is: a molecular basin (the
dense, optically thick base) and an atomic basin (the wind above the H2 -> H
front). It now tries one starting point in each in turn -- the warm start, the
chemical-equilibrium H2 fraction at the local (p, T) with the metal stages from
their own ionization balance at the incoming n_e, then the molecule-free
ionization balance of every element -- and keeps a converged root only if it is
also physical, on the same rule §45 applied to the atomic systems. The metal
part of `ionization_balance_at_fixed_ne` was split into
`metal_ionization_balance_at_fixed_ne` so the molecular retry seeds its metals
from the same balance instead of a second copy of it. The third starting point
is what the front cells needed: on the hot-Uranus case the coupled solve leaves
3 cells at r = 1.21-1.22 (the ionization front) with no admissible root at every
step with two starting points, and 0 cells with three.

**Two defects found in the surrounding code and fixed.**

- The sub-13.6 eV energy grid of a power-law spectrum was an `if
  (thereis_HeITR) ... else if (thereis_lowIP_metal)`, with the comment "HeITR +
  metals is unsupported, so at most one applies". That has not been true since
  the triplet and the metals were merged. With the triplet on, the grid floor
  stopped at 4.80 eV, which is above the K I threshold (4.341 eV), so K I would
  be held spuriously neutral in a triplet run. The floor is now the lowest
  threshold over whichever sources are active, matching what the loaded-SED
  path (`sed_read`) already did. No golden carries K, so this is byte-identical
  there.
- The `pp_metals 2` metal re-solve passed the global `N_eq` to `hybrd1` while
  its residual (`ion_system_metals_pp`) writes only rows 1..3+2*n_melem. With
  the He 2^3S triplet that already left one row of `fvec` unwritten; with the
  molecular layout it would leave four or five. It is now sized `3 + 2*n_melem`
  with a matching workspace, the same reasoning `lwa_adv` already carried a few
  lines above. Both goldens run `pp_metals 1`, so this is byte-identical there
  and changes results only for `pp_metals 2` runs.

**Left in place, with the validity range now written at the code.** `eval_cool`
builds its own electron density without the molecular ions. That is right in
the hot atomic gas it was written for and wrong inside a deep molecular base,
where H3+ can be the dominant ion; there it disagrees with the n_e the
equilibrium solver itself uses (`ionization_equilibrium` does pass `nmol_eq` to
`calc_ne`), so the ne-scaling cooling channels are under-counted in that layer.
Two definitions of the same quantity is one too many, but closing the gap
changes the molecular results and belongs in its own change with its own
measurement; the validity range is written at the call site and the item is
left open below. *Closed in §52.*

**Gates.** `make check` PASSes byte-identical on all three cases
(`wasp_full`, `wasp_he23off`, `mol_base_handoff`); a baseline `make check` on
the same tree before the change also passed, so the comparison is this change
alone. Goldens were not re-snapshotted. A build at
`-O1 -fopenmp -g -fcheck=bounds,do,mem` is clean on all three coupled layouts
for 200-300 steps: hot Uranus with the He 2^3S triplet (N_eq = 28, metals at
row 9), the same case without it (N_eq = 27, metals at row 8), and the
HD 209458 b deep handoff.

Hot Uranus (the Tier-2 gate configuration of
`docs/lower_atmosphere_figs/data_g2`) plus solar C/N/O/Mg/Ca/Na/Fe, 12000-step
relaxation snapshot against the same case metals-off: the H2 -> H front sits at
r = 1.1577 in both (the gate's 1.156), the base stays fully molecular in both,
log10 Mdot = 10.58 in both, no negative densities, no cell left without an
admissible root (the same case with the triplet off, N_eq = 27 with the metals
at row 8, is clean over a shorter 400-step smoke). What the metals change is
the base chemistry: n_e rises 1900x
and H3+ falls by a factor 1350 (7.08e4 -> 5.25e1 cm^-3), its peak moving out
from r = 1.040 to r = 1.130 -- H3+ is destroyed by dissociative recombination
(Koskinen R6/R7, both proportional to n_e), so it survives only where the metal
electron fraction starts to fall. Metal ionization runs the expected way: C II/C
= 8e-10 at the base, 0.22 at the front, 0.43 at 1.9 R_p; Mg III stays below
1e-8 up to the front and reaches 2.4e-3 only at 1.9 R_p.

HD 209458 b at the 1e-3 bar handoff (the §11.7 blocker), C/N/O, He 2^3S on,
4000-step bounded run: it integrates. §11.7 records that the atomic pair built
at this level terminated with a NaN at r = 1.067 after 1846 steps with the time
step collapsed to dtu = 1.2e-6 (that run was not repeated here); this one
reaches the cap with dtu = 0.198 and no NaN, a fully
molecular base, and the front at r = 1.0150 (1.0169 metals-off). It is a
relaxation snapshot at du = 747, nowhere near converged, so no mass-loss rate
should be read off it.

**Open item.** The `_adv` post-process is still molecule-free by design
(`post_process_adv.f90` header): `nh = nhi + nhii` counts only the free H
nuclei and the molecular electrons are absent from its n_e. The metal stages
are carried there as usual under `pp_metals`, but on that molecule-free H/He
background, so `_adv` metal profiles inside the molecular layer inherit the
approximation. With the `eval_cool` electron density closed in §52, this is the
one remaining place where a molecular run and a metal run still meet on the
atomic assumption.

Files: `src/modules/nonlinear_system_solver/System_HeH_mol_metals.f90` (new),
`src/modules/nonlinear_system_solver/System_HeH_mol.f90`,
`src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/radiation/util_ion_eq.f90`,
`src/modules/files_IO/input_read.f90`,
`src/modules/files_IO/write_setup_report.f90` (the startup report now names the
species blocks the coupled system carries and prints `N_eq`),
`src/modules/init/set_energy_vectors.f90`,
`src/modules/post_process/post_process_adv.f90`, `Makefile`,
`examples/16_molecular_metals/` (new). Documentation: `README.md`,
`README_HOWTO.md`, `examples/README.md`, `docs/EXHALE_user_manual.tex`
(the example table also became a `longtable`, which it needed already --- as a
float it ran past the bottom margin into the page footer),
`docs/input_schema.md`, `docs/lower_atmosphere_coupling.md`,
`docs/base_composition_handoff_plan.md`.

---

## 52. One electron density for the cooling and for the ionization solver (2026-08-13)

`eval_cool` rebuilt its own free-electron density and left the molecular ions
out of it. The equilibrium solver in the same sweep does not: it balances
ionization against `calc_ne(nhii,nheii,nheiii,ne,nm,nmol_eq)`, i.e. H+, He+,
He++, the metal stages under the `eos_metals` policy, **and** H2+, H3+, HeH+.
A molecular run therefore carried two definitions of n_e, and the cooling used
the smaller one. That is wrong wherever a molecular ion carries a
non-negligible share of the charge, which is exactly the deep molecular base.
Measured on the metals-off hot-Uranus Tier-2 gate, the two definitions differ by
**168.5x at the base cell** (373 vs 6.29e4 cm^-3, H3+ supplying essentially all
of it), 4.5x at r = 1.05, 2.8% at r = 1.10, and nothing above r ~ 1.15; on the
`mol_base_handoff` case the base ratio is 175.7x.

`eval_cool` now takes an optional `nmol` (H2, H2+, H3+, HeH+, cgs) and passes
it straight to `calc_ne` -- the electron density is not reconstructed anywhere
inside `eval_cool`, it has one definition and one implementation. Callers that
track the molecular network supply it: `ionization_equilibrium` (`nmol_eq`),
the semi-implicit energy solver (`nmol_dim`, the cgs copy of the adimensional
array it already hands `calc_ne`/`calc_ntot`), and the `Cooling_breakdown`
dump. The `_adv` advection post-process deliberately does not: that
reconstruction is molecule-free by design (its module header), and omitting the
argument gives it the atomic charge sum it is built on.
`ionization_equilibrium`'s own `calc_ne` call lost its `thereis_mol` branch at
the same time -- `nmol_eq` is zero for an atomic run, so one call covers both.

Every n_e-scaling channel picks the correction up: recombination, collisional
ionization, collisional excitation and bremsstrahlung are linear in n_e, and
the coronal metal line coefficients enter as `ne*nm*c_metal`. The saturated
ground-term coefficients (Fe II, and C I/C II/N II/O I under `cno_chianti`)
are `Lambda_eff = W/ne`, so the assembly's n_e cancels and they respond only
through the level populations (`C21 ~ ne`) -- the correct behavior, and there
is no double counting of n_e anywhere in `eval_cool`.

**Atomic runs are exactly unchanged.** `nmol_eq` and the molecular `f_sp`
columns are zero, and `calc_ne` adds them as `ne + 1.0*0.0`, which is exact.
`make check` is byte-identical on `wasp_full` (metals + He 2^3S) and
`wasp_he23off`; a full WASP-121b run (metals, He 2^3S, Ly-alpha escape
probability, excited H, transonic IC) reproduces all eight output files
byte-identically -- `Hydro_ioniz{,_adv}`, `Ion_species{,_adv}`,
`Cooling_breakdown`, `Heating_breakdown`, `Excited_H`, `IC_dump` -- at
log10 Mdot = 13.20.

**`mol_base_handoff` changes, and that is the point of the change.** At the
pinned 12000-step snapshot the two runs stop at the same `du = 2.6689` and the
same `dtu`, every hydrodynamic column agrees to <= 4e-8 relative and every
species column to <= 4.2e-6 (that worst case is H3+ at r = 2.96, where H3+ is
3e-19 cm^-3); the one column that moves is the reported cooling rate at the
near-neutral base, 2.799e-20 -> 4.942e-18 erg cm^-3 s^-1 (x176.5, the n_e
ratio, as expected for channels linear in n_e). log10 Mdot stays 10.58. The
golden was NOT re-snapshotted here.

**Why the dynamics barely notice.** A molecular layer cools through the H3+
infrared lines, which do not scale with n_e, so multiplying the electron
channels by 176 at the base moves the cooling from 2.8e-20 to 4.9e-18 against
1.3e-7 erg cm^-3 s^-1 of photoheating there. The correction is real and it is
now the same quantity the solver balances ionization against, but it is
energetically negligible in the layer where it is large.

**Molecular smoke tests** (hot Uranus, the Tier-2 gate configuration,
12000-step relaxation snapshots -- both runs stop at `du` of order unity, so
these are snapshots and not converged solutions):

| | metals off | metals on (solar C/N/O/Mg/Ca/Na/Fe) |
|---|---|---|
| n_e(base) old -> new | 373 -> 6.29e4 cm^-3 (x168.5) | 1.2016e8 -> 1.2016e8 (x1.0000) |
| cooling at the base | 3.18e-20 -> 5.36e-18 erg/cm3/s | +7e-7 relative |
| H2 -> H front | 1.15631 -> 1.15631 | 1.15677 -> 1.15677 |
| H3+ peak | 7.0808e4 -> 7.0807e4 cm^-3 at r = 1.0396 | 52.465 -> 52.465 at r = 1.1304 |
| `du` at step 12000 | 4.3354 -> 4.2914 | 5.8177 -> 5.8177 |
| log10 Mdot | 10.58 -> 10.58 | 10.58 -> 10.58 |

The metals-on column is the expected null: the low-IP metals hold n_e at 1.2e8
at that base (section 51), against 52 cm^-3 of H3+, so adding the molecular
ions moves nothing (max relative change 7e-10 in rho/T/p, 3e-9 in the metal
stages). Metals-off, the snapshot differs by 1% in `du` and by sub-percent in
the bulk state just above the layer (rho -0.7%, T +0.8% at r = 1.19), with the
trace species that depend exponentially on T amplifying that into H3+ -16%,
H2 -13%, He II +10% at the same radius. Those are differences between two
snapshots of an unconverged relaxation, not a converged shift; the gate
observables -- front position, H3+ peak, mass-loss rate -- are unchanged.

On the HD 209458 b 1e-3 bar handoff with C/N/O (4000-step bounded snapshot,
also unconverged) the base picture is the same: the metals dominate n_e, the
base cooling rises 0.45%, the molecular-layer cooling by at most 1.9%, the
front moves 1.01503 -> 1.01518 and the H3+ peak by 0.05%.

An `-O1 -fopenmp -g -fcheck=bounds,do,mem` build runs all three molecular
layouts (hot Uranus with and without metals, HD 209458 b deep handoff) for 300
steps with no runtime check firing.

**Two related items found and recorded, not changed here.**

- The H3+ infrared cooling is added to `cool` in `ionization_equilibrium`,
  *after* `eval_cool` returns. The JFNK steady residual uses that array, but
  the semi-implicit energy update (the default) rebuilds `cool` from
  `eval_cool` alone and overwrites it, so it never sees the H3+ term. For a
  molecular run the marching stage and the Newton stage therefore integrate
  slightly different energy equations, and the marching fixed point is not the
  zero of the residual the Newton solver drives down. The fix is to move the
  term inside `eval_cool` (which already receives `nmol`) with its own
  `cool_chan` column; that changes the `Cooling_breakdown.txt` schema and the
  molecular results, so it belongs in its own change. The `Cooling_breakdown`
  header now states that the channel sum falls short of the `cool` column by
  that term in the molecular layer.
- H2+, H3+ and HeH+ donate their electron to n_e but are left out of the
  bremsstrahlung charge sum, so that sum is smaller than n_e inside the
  molecular layer. They are Z_ion = 1 and belong there formally, but they only
  exist at T ~ 1e3 K, where the hot-plasma free-free expression is an
  extrapolation emitting at frequencies the atmosphere is not thin to, and
  where free-free is negligible against the H3+ and metal line cooling. The
  scope is written at the accumulator.

Files: `src/modules/radiation/util_ion_eq.f90`,
`src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/time_step/energy_semi_implicit.f90`,
`src/modules/post_process/post_process_adv.f90` (comments only).
Documentation: `docs/lower_atmosphere_coupling.md`.

---

## 53. A fourth regression case for the molecular + metals system (2026-08-13)

The coupled molecular + metals system of §51 was not guarded by anything in the
matrix: `mol_base_handoff` runs the molecular network metals-off, and the two
`wasp_*` cases run metals without molecules, so the merged residual, the
`metal_row_base()` offset and the three-starting-point seeding had no golden.

`backup/regression/mol_metals` is that case: the same hot-Uranus Tier-2 gate
configuration as `mol_base_handoff` (`Molecular chemistry: True`, its `base.inp`
handoff, 12000 steps pinned in `maxsteps`) plus a solar-abundance
`metals.inp` with the seven elements C/N/O/Mg/Ca/Na/Fe. It stops at
`du = 3.1770`, `log10 Mdot = 10.58`. `DEFAULT_CASES` in
`backup/regression/run_check.sh` is now
`wasp_full wasp_he23off mol_base_handoff mol_metals`, and the case list in the
script header describes what each one guards.

Goldens for all four cases were re-snapshotted at the end of this change
series, so they carry §51 and §52. The two atomic cases are byte-identical to
their previous goldens (both changes are exact no-ops without molecules);
`mol_base_handoff` moves only in its cooling column, by the factor 176.5 of
§52.

One harness trap, hit while doing this and now written into the script: `golden`
does **not** re-run the cases, it copies whatever `output/` currently sits in
each case directory. Snapshotting straight after a code change baselines the
*previous* binary. The order is `check` (which rebuilds and re-runs), then
`golden`, then `check` again.

Files: `backup/regression/run_check.sh`, `backup/regression/mol_metals/` (new),
`backup/regression/golden/` (re-snapshotted). Documentation: the workspace
`CLAUDE.md`, `README_HOWTO.md`.

---

## 54. H3+ infrared cooling inside `eval_cool`: one energy equation for a molecular run (2026-08-13)

Closes item (F) of `TO_BE_DONE.md`, opened by the audit in section 52.

The H3+ infrared cooling was added to `cool` in `ionization_equilibrium`, *after*
`eval_cool` returned. Only some of the consumers of `cool` saw it:

| term in `cool` | rebuilt by `energy_semi_implicit` (the default marching temperature update) | returned by `ioniz_eq` (JFNK steady residual, the explicit-energy option, the in-loop `||R||` monitor, the `Hydro_ioniz.txt` cool column) |
|---|---|---|
| H/He recombination | yes | yes |
| collisional ionization (H I, He I, He II, He 2^3S) | yes | yes |
| collisional excitation (H I Ly-alpha, He I, He II, He 2^3S 10830 + 2^1S/2^1P conversion) | yes | yes |
| bremsstrahlung (H+, He+, He++, metal ions) | yes | yes |
| metal line cooling (27 ions, saturated ground-term fine structure included) | yes | yes |
| **H3+ infrared (Miller+2013)** | **NO** | **yes** |

`heat` has no such split: it is assembled once in `ioniz_eq` (photoheating per
absorber including H2 and the metals, Shull & van Steenberg secondary-ionization
partition, Balmer photoelectric and Ly-alpha de-excitation, He-recombination
coupling, He 2^3S + H and He 2^3S + H2 Penning) and passed into
`energy_semi_implicit` as `intent(in)`, never rebuilt. H3+ was the only
asymmetric term, and it was the whole asymmetry.

Three consequences followed from that one row.

- **The marching fixed point was not the zero of the residual the Newton solver
  drives down.** The marching relaxed `heat = cool_atomic+metal`; the residual
  demanded `heat = cool_atomic+metal + Lambda_H3+`.
- **The two energy-solver options integrated different equations.**
  `Energy solver: Explicit` used the `ioniz_eq` array directly and did carry the
  term; the semi-implicit default (the production setting) did not.
- **The `cool` column of `Hydro_ioniz.txt` meant different things depending on
  how the run stopped.** A `du`-stop left the semi-implicit array (no H3+) in
  place; the JFNK finish calls `ioniz_eq` after `solve_steady_jfnk` and left the
  array that has it.

**Size of the missing term.** On the metals-off hot-Uranus Tier-2 gate
(`backup/regression/mol_base_handoff`) the H3+ lines carry **more than 99% of
the total radiative cooling from the base out to the H2 -> H front**
(`Cooling_breakdown` H3p_IR / cool_total > 0.5 over r = 0.9998-1.1856, > 0.99
over 0.9998-1.1557), peaking at 9.35e-7 erg cm^-3 s^-1 at r = 1.0413, which is
**7.3x the local photoheating**. Below the front every atomic channel is
exponentially suppressed at T ~ 1.3e3 K, so the marching stage was integrating
that layer with essentially no radiative loss at all: the cooling it used was
smaller than the true total by a factor 1.1e2 to 1.2e11 over r = 1.00-1.156.
With solar trace metals (`mol_metals`) the metal lines dominate instead and H3+
reaches only 22% of the local cooling, in a thin band r = 1.066-1.193.

**The fix.** `eval_cool` computes the H3+ cooling itself, from the `nmol` array
it has taken since section 52, using the same `h3p_cooling` module (one
definition of the Miller+2013 fits and the Table-6 non-LTE factor; `eval_cool`
`use`s it, nothing is re-derived). The post-return addition in
`ionization_equilibrium` is gone, and a comment at the `eval_cool` call states
why nothing may be added there. `cool_chan` gains column 7 for the channel, so
the decomposition still reproduces `cool` exactly (measured
`max |sum(channels)/cool - 1| = 3.3e-16` on the molecular gate, 5.5e-16 on
`wasp_full`). Callers that model a molecule-free gas (`post_process_adv`, and
`T_equation` which solves for the `_adv` temperature) omit `nmol` and therefore
still get the atomic charge sum *and* no H3+ term -- the same approximation,
stated at both call sites.

A side effect worth naming: `energy_semi_implicit` builds its damping
derivative `dC/dT` by finite-differencing `eval_cool` in temperature, so the
steep temperature dependence of the H3+ emission now enters the implicit
temperature update instead of being invisible to it.

**Output schema.** `Cooling_breakdown.txt` gains `col11 H3p_IR`; the 27
metal-ion columns move from 11-37 to 12-38. It is written for every run and is
zero without the molecular network. `python/paper_data.py` reads the column
names from the header and follows automatically; `examples/exhale_io.py` has a
fixed list, which now carries `H3p` (and is named `COOL_GAS_CHANNELS`, since it
is no longer only the H/He channels).

**Atomic runs are exactly unchanged.** `nmol(:,3)` is zero, the term is skipped
cell by cell, and the total is `... + 0.0`. `make check` is byte-identical on
`wasp_full` (14060 steps, metals + He 2^3S) and `wasp_he23off` (14037 steps),
with the same step counts as their goldens.

**The two molecular cases change, which is the point.** Both were re-run at
their pinned 12000-step convention; the goldens were **not** re-snapshotted.

| | metals off (`mol_base_handoff`) | metals on (`mol_metals`) |
|---|---|---|
| `du` at step 12000 | 2.6689 -> 2.6274 | 3.1770 -> 3.1769 |
| log10 Mdot | 10.9043 -> 10.9040 | 10.9105 -> 10.9105 |
| H2 -> H front (2 n_H2 / n_H = 0.5) | 1.16087 -> 1.16039 | 1.16246 -> 1.16246 |
| H3+ peak | 7.0237e4 -> 7.0285e4 cm^-3 at r = 1.0419 | 52.139 -> 52.140 at r = 1.1355 |
| base T | 1213.42 -> 1213.42 K | 1213.46 -> 1213.46 K |
| max T change, whole grid | +1.4% at r = 1.93 (wind transient) | +5.2e-5 at r = 1.20 |
| T change where H3+ dominates | -0.53% at r = 1.078, mean -0.20% over r <= 1.156 | -3.0e-5 at r = 1.183 |
| reported cooling at the base | 4.94e-18 -> 5.97e-7 erg cm^-3 s^-1 (x1.2e11) | unchanged |
| max cooling change | x1.23e11 at r = 1.0021 | +28% at r = 1.1597 |

**Why the gate observables barely move, and why that is not evidence the term is
small.** The molecular layer is thermally massive: its H3+ cooling time
`(3/2) n k T / Lambda_H3+` is 3.7e6 s at r = 1.04 and 2.2e7 s at the base, i.e.
**200-2000 sound-crossing times `t_s = R0/v0`**. A 12000-step snapshot advances
a small fraction of that, so the layer has only begun to respond: T has fallen
by 0.1-0.5% through r = 1.02-1.12 and is still falling. The front, the H3+ peak
and Mdot are unchanged to < 0.1% *at this step budget*, so the recorded Tier-2
gate numbers stand -- but the converged molecular solution now has a different
thermal structure than the pre-change code would have produced, and reaching it
needs a step budget the gate convention does not provide. When this was written no
molecular configuration in the tree had converged (`examples/15_molecular` had
run 135694 steps to `du = 2.9e-2`), so it was stated as a limit, not a result;
the paragraph below closes it.

**What was verified and what was not.** The residual/Newton path is unchanged by
the move: `EXHALE_RESIDUAL=1` on one frozen molecular state gives a
**bit-identical** `output/residual_profile.txt` before and after, so only the
marching side moved. What could *not* be shown is the H3+ term disappearing from
the steady energy residual: on the only molecular states available -- 12000-step
relaxation snapshots with a base-breathing transient, volume-weighted
`||R|| ~ 3` -- the cell-by-cell energy residual in the molecular layer is 1.5-7x
larger than `Lambda_H3+/q0`, so the systematic offset is buried in the transient.
That check needs a converged molecular run.

**Closed 2026-08-13: the check on a converged molecular run.** Restarting the
metals-off hot-Uranus snapshot with `Load IC? True`, `Solver: Newton 5.0e-2` and
`Max steps: 150000` reaches `info = 0` after a JFNK finish at step 22556 -- the
first Newton-grade molecular solution the code has produced. In the 180 cells
where H3+ carries more than half the cooling (the cells at `r <= 1.001`, whose
residual sets `||R||`, excluded) the ratio of the cell-by-cell energy residual to
`Lambda_H3+/q0` falls from **4.42 at the 12000-step snapshot (sign 89+/166-) to
0.70 at the default `||R|| < 1e-3` solution (173+/7-) to 0.0017 once the same
state is re-solved at `Resid tol: 1.0e-5` (25+/79-)**: at the converged solution
the residual carries no offset of the size of the H3+ term, and no preferred sign.
Two things follow. First, the default `||R|| < 1e-3` is not a converged molecular
run -- the norm is set by two or three cells above the base that carry ~10^3 times
the residual of anything else, so the molecular layer is still cooling when the
solver stops; `Resid tol: 1.0e-5` is required, and a third solve from that state
reproduces it (T to 0.14%, front to five digits, `||R||` 6.7e-7). Second, an A/B
against the pre-change binary isolates the term directly: marching both binaries
12000 steps from the same converged state, the temperature difference
`T(before) - T(after)` is proportional to the local `Lambda_H3+/(3/2 n k)` with a
single constant: over the 155 cells where `Lambda_H3+` is within two decades of
its peak the ratio has median 1.27 and 10-90% range 1.13-1.32, and over all 187
cells where H3+ carries more than half the cooling the sign is right in 181. The
whole difference between the two marching energy equations is the H3+ term and
nothing else. The physics of the converged solution, including the 190-300 K collapse of
the molecular layer that the optically thin line cooling produces once the
marching feels it, is in `docs/lower_atmosphere_coupling.md`.

An `-O1 -fopenmp -g -fcheck=bounds,do,mem` build runs both molecular
configurations for 300 steps with no runtime check firing.

Files: `src/modules/radiation/util_ion_eq.f90`,
`src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/post_process/post_process_adv.f90` (comments only),
`examples/exhale_io.py`, `python/paper_data.py`,
`python/make_struct_figures.py`. Documentation: `README.md`,
`docs/EXHALE_user_manual.tex` (section 4.3 rewritten -- it still described four
H/He channels and columns 9-35 for the metals, which had been wrong since the
collisional-excitation split; a radiative-cooling row for H3+ added to the
atomic-data table), `docs/lower_atmosphere_coupling.md`, `TO_BE_DONE.md`.

## 55. The infrared field of the lower atmosphere (`Base IR field`) (2026-08-13)

Narrows item (G) of `TO_BE_DONE.md`, opened by section 54.

The infrared coolants of a molecular base were radiating into a vacuum they are
not in. The ground-term fine-structure lines carried an escape probability built
from the OUTWARD column only, with the downward direction dismissed as "absorbed
by the lower atmosphere, which the model treats as a fixed reservoir"; the H3+
bands carried no transfer at all. Absorption of what that reservoir radiates back
was in neither.

**Measured first.** Line-center depths on the converged metals-on hot-Uranus
solution (level populations and Doppler widths from the cooling module itself):

| direction | C I 609 um | C I 370 um | O I 63 um | O I 145 um |
|---|---|---|---|---|
| outward, base cell | 0.088 | 0.172 | 2.01 | 0.54 |
| outward, coldest cell (r = 1.083) | 1.4e-3 | 2.6e-3 | 0.028 | 4.6e-3 |
| one scale height below the base | 0.046 | 0.092 | **1.08** | **0.32** |
| pressure where the downward tau = 1 | 205 ubar | 106 ubar | **17 ubar** | **37 ubar** |

The base sits at 9.0 ubar, so every line that carries the cooling is black
downward within one to three scale heights of unresolved atmosphere, while
outward the collapsed layer is thin (beta >= 0.90 for C I, >= 0.63 for O I above
r = 1.05). Trapping was never the missing piece; the incident field was. The
ratio of the incident mean intensity to the line source function,
(1/2)B_nu(T_base)/B_nu(T), passes 1 below T ~ 610 K and reaches 3.4-5.7 at the
coldest cell: those lines were cooling a gas they should have been heating.

**And it refutes the diagnosis in item (G) for the metals-on case.** The specific
entropy p/rho^gamma RISES monotonically outward through the whole "collapsed"
layer, 1.00 at the base to 1.67 at r = 1.05 and 5.69 at r = 1.083; the gas is net
heated everywhere and is cold only because rho has fallen by 1e2. Its cooling is 36% of the local
photoheating at r = 1.01 and 2-16% above r = 1.02, and heat - cool is balanced by
rho v [de/dr + p d(1/rho)/dr] to a few percent. And the metals-off solution, with
no fine-structure cooling at all, collapses further (115 K at r = 1.033, total
cooling 3.8e-12 erg cm^-3 s^-1). The metals-on layer is an adiabatic expansion,
not a radiative collapse. Radiation does dominate in one place: the metals-off
base, where H3+ carries 100% of the cooling and the entropy FALLS by a factor 3
between the base and r = 1.02.

**What was implemented.** `Base IR field: True` (default `False`) has the lower
atmosphere radiate `B_nu(T0)` over the sky fraction `f = 1 - sqrt(1-(R_p/r)^2)`
and both coolant families absorb it.

- Fine structure: the field enters INSIDE the ground-term statistical
  equilibrium as the photon occupation number
  `nbar = beta_1(tau_dn) f/(exp(E_ul/T0)-1)`, with radiative rates
  `beta A (1+nbar)` down and `beta A (g_u/g_l) nbar` up, so the returned power is
  the net one. The escape probability becomes two-sided,
  `beta = beta_1(tau_up) + beta_1(tau_dn)`, `beta_1` the single-face
  Hollenbach & McKee (1979) / de Jong et al. (1980) form; with the field off it
  is `2 beta_1(tau_up)`, the previous value, bit for bit. Each line stops cooling
  at its own radiative equilibrium temperature -- 575.8 K ([C I] 609 um) to
  641.7 K ([O I] 44 um) for `T0 = 1140` K and half-sky coverage.
- H3+: Miller et al. (2013) fit the TOTAL emission, so the exchange is closed on
  ONE effective band, the nu2 fundamental at 2521.3 cm^-1 (E/k = 3627.5 K):
  `Lambda_net = Lambda_emit(T)[1 - nbar exp(E/T)]`, equilibrium at 936.1 K at the
  base. The approximations and their range are written at
  `h3p_net_cooling_rate`.
- Not applied, with the reason at the code: Fe II (a precomputed
  statistical-equilibrium table, which cannot take a field as an argument), the
  coronal remainders and every other metal ion (fits to CHIANTI sums with no line
  list; their lowest terms are at E/k >~ 1.6e4 K, and the heating they would add
  was estimated at ~0.4% of the losses of this layer), and the legacy
  `use_2lev_cool` branch.

**Effect.** Restarting the two converged `Resid tol: 1.0e-5` molecular solutions
with the switch on, both returning `info = 0`:

| | metals off, off -> **on** | metals on, off -> **on** |
|---|---|---|
| T at r = 1.005 | 626.8 -> **910.0 K** | 1099.0 -> 1098.9 K |
| T at r = 1.02 | 299.6 -> **859.4 K** | 928.9 -> 928.9 K |
| T at r = 1.03 | 140.2 -> **840.2 K** | 820.0 -> 820.1 K |
| coldest cell | 115.1 K @ 1.0334 -> **229.7 K @ 1.0947** | 187.9 -> 189.1 K @ 1.0831 |
| H2 -> H front | 1.03582 -> **1.08833** | 1.07976 -> 1.07979 |
| log10 Mdot [g/s], spherical | 10.539 -> **10.617** | 10.596 -> 10.596 |

The metals-off layer settles ON the H3+ radiative equilibrium curve: predicted
911/886/874 K at r = 1.005/1.02/1.03 with the local sky fraction, solution
910/859/840 K, just below by the margin the expansion takes out. Its H3+ channel
is now a heating term (-1.1e-7 erg cm^-3 s^-1 at r = 1.02). The solution is a
fixed point: re-solving returns `info = 0` at `||R|| = 8.8e-6` and reproduces T
below r = 1.3 to 1e-4. The metals-on layer does not move, which is the measurement
above, not a failure of the closure.

Atomic runs are unaffected: the switch is off by default, and even on it touches
only channels that carry no cooling in an atomic wind. `make check` is
byte-identical on all four cases.

Files: `src/modules/radiation/Cool_coeff.f90`,
`src/modules/lower_atmosphere/h3p_cooling.f90`,
`src/modules/radiation/util_ion_eq.f90`,
`src/modules/nonlinear_system_solver/T_equation.f90`,
`src/modules/post_process/post_process_adv.f90`,
`src/modules/init/parameters.f90`, `src/modules/files_IO/input_read.f90`,
`src/modules/files_IO/write_setup_report.f90`. Documentation:
`docs/lower_atmosphere_coupling.md`, `docs/EXHALE_user_manual.tex`,
`TO_BE_DONE.md`.

## 56. H2 photodissociation in the Lyman-Werner bands (`Stellar LW flux`) (2026-08-13)

The molecular network had no photodissociation of neutral H2. The Koskinen et
al. (2022) Table-1 photo-rates start at the 15.4 eV photoionization edge, so
below that the only H2 losses were thermal (R12), electron impact (R14) and ion
chemistry. `docs/base_composition_handoff_plan.md` section 5 named the omission
as the reason EXHALE's base comes out more molecular than any photochemical code
gives, and recorded that pinning `q_H2` at the base (Route 1) patches over it.
This adds the missing term.

**The band is not in the code's own radiation field, and had to be.** Checked in
the source: `set_energy_vectors` builds the photon grid from 13.6 eV up and
reaches below the H I edge only when the He 2^3S metastable (4.8 eV) or a low-IP
metal is active, and then only to that species' threshold; `read_sed` selects
SED rows by the same `e_low`. For a numerical SED the 11.2-13.6 eV band is
therefore not read at all, and for a power-law SED (every Tier-2 gate case) it
is the XUV power law extrapolated downward, which has nothing to do with a
star's FUV. Even where the grid does reach into the band, its opacity is
continuum photoionization, whereas Lyman-Werner absorption is a forest of
saturated lines. The band flux is a separate input:

```
Stellar LW flux [erg/cm2/s]: 343.0    # 912-1110 A, integrated, at the planet
```

Default 0 = off; inert (with a warning) when `Molecular chemistry` is off.

**Rate.** Draine & Bertoldi (1996), ApJ 468, 269 (published version). For their
flat-F_lambda spectrum (u_nu ~ nu^-2, their eq. 24) at chi = 1 the band photon
flux is F = 1.208e7 cm^-2 s^-1 (Table 1) and the unshielded dissociation rate is
zeta_diss(0) = 4.17e-11 s^-1 (Table 2, Fig. 7 caption), so the rate per band
photon is an effective cross section sigma_LW = 3.452e-18 cm^2. With the mean
photon energy of a flat-F_lambda band, 2hc/(912 + 1110 A) = 12.2635 eV,

    k_LW,thin = 1.757e-7 * F_LW    [s^-1, F_LW in erg cm^-2 s^-1] .

Their Table 1 was reproduced numerically from their eqs. (22)-(24) to confirm
the flux convention (Habing 1.2220e7 vs 1.222e7; nu^-2 1.2084e7 vs 1.208e7;
Draine 1978 1.2313e7 vs 1.232e7). The band-shape sensitivity is +-25%: the same
arithmetic for the soft Draine 1978 field gives 2.685e-18 cm^2.

**Self-shielding** is their eq. (37), the fit to the full multiline calculation
including line overlap, valid over 1e14 < N_H2 < 3e21 cm^-2, with b the thermal
H2 Doppler parameter (2kT/m_H2)^1/2. N_H2 is the star-ward column, built by the
same `calc_column_dens_one` radial integration and `opa_pf` weighting as every
other absorber column, from the incoming (pre-solve) H2 density. No dust
(EXHALE's metals are atomic and trace), and the trace-metal continuum in the
band is tau ~ 4e-23 N_H, negligible where the shielding matters; the H
Lyman-series lines inside the band are inside DB96's own fit and not treated
separately.

**Chemistry and heating.** H2 + hv -> H + H enters the H2 balance row of the
coupled system next to the H2 photoionization (`mol_heh_rows`, shared by
`System_HeH_mol` and `System_HeH_mol_metals`), so it is solved with everything
else. Both products are neutral H, which the H-nucleus closure supplies; no
other row changes. The fragments carry 0.4 eV of kinetic energy (Black &
Dalgarno 1977, ApJS 34, 405, p. 418, "the yield is about 0.4 eV per atom pair"),
added to `heat`. The 4.48 eV bond energy is paid by the photon, not by the gas,
and is not a sink of this channel.

`output/Lyman_Werner.txt` (r, T, x_H2, n_H2, N_H2, f_shield, k_LW, heating) is
written whenever a molecular run carries a band flux.

**Band flux for a planet.** Integrated externally from the VULCAN stellar
spectra in the tree (`EXHALE/VULCAN/atm/stellar_flux/`, stellar surface) over
91.2-111.0 nm and diluted by (R_star/a)^2. The `VULCAN_run_*` spectra directories
are empty since the 2026-08-09 loss, so no HD 209458 spectrum exists here and the
solar spectrum stands in: 343 erg cm^-2 s^-1 at 0.048 AU behind a 1.155 R_sun
star, k_LW,thin = 6.03e-5 s^-1. The same file returns a solar Lyman-alpha
irradiance of 6.12 erg cm^-2 s^-1 at 1 AU against an observed 6-8, which is the
calibration check. HD 189733 b comes out at 600 (Moses11 spectrum) or 1830
(Bourrier 2020) erg cm^-2 s^-1 -- a factor 3 that dwarfs the +-25% of the band
shape.

**Effect.** A/B on the key alone, hot-Uranus gate (metals off, the
`mol_base_handoff` configuration restarted with `Solver: Newton 5.0e-2`,
`Resid tol: 1.0e-5`, `Base IR field: True`), both returning `info = 0`:

| | LW off | **LW on (343)** |
|---|---|---|
| residual norm at exit | 8.15e-6 | 4.86e-6 |
| H2 -> H front (x_H2 = 0.5) | 1.08806 | **1.08622** |
| x_H2 at the base | 0.99921 | 0.99883 |
| x_H2 at r = 1.10 | 2.54e-2 | **1.51e-4** |
| H3+ peak [cm^-3] | 3.19e5 @ 1.0079 | 3.20e5 @ 1.0025 |
| log10 Mdot [g/s], spherical | 10.6196 | 10.6248 |

and on HD 209458 b (`examples/15_molecular` + `Base IR field: True`, both members
at `Resid tol: 2.0e-5` because the band-on run stalls the JFNK at 1.5e-5 on the
base-adjacent momentum row, item (A) of `TO_BE_DONE.md`; both `info = 0`, both
mass fluxes flat to 4.1e-3): front 1.01578 -> 1.01290, x_H2 at r = 1.02
2.37e-4 -> 4.22e-6, log10 Mdot 10.5711 -> 10.5600.

The self-shielding is doing exactly what it should: at the hot-Uranus base
N_H2 = 4.3e21 cm^-2 gives f_shield = 9.8e-7 and k_LW = 5.9e-11 s^-1, while above
r = 1.15 the column has fallen below 1e10 cm^-2, f_shield = 1.000 and k_LW
recovers the unshielded 6.026e-5 s^-1.

**It does not close the base-composition gap, and the measurement says why.**
The base q_H2 moves from 0.8618 to 0.8612 against the 0.75 the case's `base.inp`
carries -- 0.5% of the gap -- even though Lyman-Werner becomes the largest single
H2 loss at the base (5.9e-11 s^-1 against 2.5e-11 for H+ + H2). The two facts are
consistent because the base partition is a formation-destruction balance in which
the three-body R15 sets n_H2 ~ n_H^2, so n_H/n_H2 scales as the square root of
the destruction rate: the measured 2.3x rise in destruction raises the atomic
fraction only 1.5x. Reaching q_H2 = 0.75 needs a destruction rate 4100x larger.
The unshielded band would supply it many times over; the column stops it. The
extra atomic H photochemical codes find at 1 microbar comes from catalytic cycles
on O, OH and H2O, and EXHALE's H2/H2+/H3+/HeH+ network has nowhere to put them.
Route 1 (`q_H2_base`) therefore stays the only way EXHALE can carry a
photochemical base partition, and is now a documented modeling choice rather
than a patch over a missing term.

`make check` is byte-identical on all four cases with the key absent.

Files: new `src/modules/lower_atmosphere/lyman_werner.f90`;
`src/modules/init/parameters.f90`, `src/modules/files_IO/input_read.f90`,
`src/modules/files_IO/write_setup_report.f90`,
`src/modules/files_IO/write_output.f90`,
`src/modules/nonlinear_system_solver/ion_cell_state.f90`,
`src/modules/nonlinear_system_solver/System_HeH_mol.f90`,
`src/modules/nonlinear_system_solver/System_HeH_mol_metals.f90`,
`src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/lower_atmosphere/mol_rates.f90`,
`src/modules/radiation/util_ion_eq.f90`, `Makefile`. Documentation:
`docs/lower_atmosphere_coupling.md`, `docs/base_composition_handoff_plan.md`,
`docs/EXHALE_user_manual.tex`, `docs/input_schema.md`, `README.md`,
`README_HOWTO.md`.

---

## 57. A fifth regression case for Lyman-Werner photodissociation (2026-08-13)

Nothing in the four-case matrix exercised the code the section 56 key turns
on: with `Stellar LW flux` absent every case leaves `F_LW_star = 0`, so the
star-ward H2 column, the Draine & Bertoldi self-shielding factor, the LW loss
row in the H2 budget, the 0.4 eV fragment heating and the `Lyman_Werner.txt`
writer all ran in no case.

`backup/regression/mol_lyman_werner` is that case: a copy of
`mol_base_handoff` (the hot-Uranus Tier-2 gate, `Molecular chemistry: True`,
the same `base.inp` handoff, 12000 steps pinned in `maxsteps`) plus the single
line `Stellar LW flux [erg/cm2/s]: 343.0` — the HD 209458 b solar-proxy band
flux of section 56. It stops at `du = 2.1841`, `dtu = 4.8072e-3`, clearly
apart from the LW-off `mol_base_handoff` (`du = 2.6274`, `dtu = 5.7385e-3`):
the case does pass through the new code path. `DEFAULT_CASES` in
`backup/regression/run_check.sh` is now `wasp_full wasp_he23off
mol_base_handoff mol_metals mol_lyman_werner`, and the script header describes
what the case guards.

The golden was taken in the order the section 53 trap requires — `check`
(which rebuilds and re-runs the case), then `golden`, then the full `check`
again — and the final run is byte-identical on all five cases.

Files: `backup/regression/run_check.sh`,
`backup/regression/mol_lyman_werner/` (new),
`backup/regression/golden/mol_lyman_werner/` (new). Documentation: the
workspace `CLAUDE.md`, `README_HOWTO.md`.

---

## 58. JFNK convergence is judged on the iterate the solver returns; `load_IC` no longer accepts an all-zero metal column (2026-08-13)

Two things came out of trying to converge HD 209458 b with molecular chemistry
AND solar C/N/O (`examples/16_molecular_metals`), which fails with the
prescription that works for its metals-off twin `examples/15_molecular`. The
diagnosis is in `docs/hd209_metal_stagnation.md`; what follows is what changed
in the code, and one thing that was measured and deliberately *not* changed.

**The measurement.** `solve_steady_jfnk` accepts a step on the whole-domain
scaled merit `||F/D||_2` over cells `1..N` and decides `info` on `||R||`, the
volume-weighted relative residual over the escape region `[j_min:N]` alone. On
this configuration the two move in opposite directions. Warm-restarted from its
own stalled state, the solve reached `||R|| = 3.958e-5` at outer iteration 2 and
then, over 115 further iterations, cut the merit from 5.67e-1 to 1.58e-2 — a
factor 36 — while `||R||` rose to 9.5e-5. The reason is a shell just above the
base that stops flowing once the metal cooling is on: at `r = 1.02 R_p` cooling
exceeds heating by a factor 9.7 with C/N/O against heating exceeding cooling by
a factor 33 without them, the mean speed over `r = 1.015-1.030` falls from 207
to 6-24 cm/s, and an undamped `2 dr` contact mode grows there whose amplitude
reaches several times the mean velocity. That shell dominates the whole-domain
merit and is invisible to `||R||`.

**Aligning the merit's region was tried and rejected, for the second time.**
Replacing `||F/D||_2` over `1..N` by the volume-weighted RMS of `F/D` over
`[j_min:N]` breaks solves that converge today: `examples/15` goes from
`info = 0` at `||R|| = 9.542e-6` in 272 outer iterations to `info = 2` with a
best of 1.171e-4 in 30, and the HD 209458 b `base.inp` continuation
`photo_deep_secion_cont` from `info = 0` at 7.284e-4 in 4 iterations to
`info = 2` at 3.249e-3 in 15. This reproduces on two more cases the measurement
that rejected the same change on 2026-08-11
(`docs/newton_scaling_and_base_wall.md` section 4), and for the reason given
there: the Newton step is computed from the full residual over all rows, so a
merit that ignores most of those rows rejects the steps that step takes. The
merit stays whole-domain, and the code now says so at the site with the numbers.

**What did change: `info` is decided on the iterate that is returned.** The
solver already tracked the best iterate by `||R||`, already returned it when the
solve failed, and already tested `||R|| < resid_tol` at the top of every outer
iteration — so an iterate that satisfies the tolerance is normally caught the
moment it is produced. Two gaps remained. The restore of the best iterate was
gated on `info /= 0`, and the tolerance was never tested after the loop, so a
solve whose last allowed iteration met the tolerance fell out with `info = 1`.
Both close with one line: `if (rnorm < resid_tol) info = 0` after the loop, with
the best-iterate restore made unconditional. It matters because `EXHALE_main`
treats `info /= 0` as a failure, disables the Newton finish for the rest of the
run, and falls back to a `du`-stopped march — discarding a state that satisfied
the convergence criterion.

*Scope, measured, not inferred*: this changes no trajectory, and it is not what
unlocks `examples/16`. The stalled `examples/16` warm restart was re-run at
`Resid tol: 5.0e-5` with the binary before and the binary after: both hand off
at marching step 9363, both take 2 outer iterations, both print the same two
iterate lines to every digit, and both end `info = 0` at
`||R|| = 3.958e-5`. The `examples/16` failures at `Resid tol: 1.0e-5` are
therefore not a bookkeeping artifact either — that configuration genuinely does
not reach 1e-5, for the reasons in the memo. What the fix removes is the
remaining way for a returned state that satisfies the tolerance to be labelled a
failure.

**`load_IC` no longer accepts a metal column that is present but zero.** It
decided whether the restart file carried an element from the presence of its
columns alone. The schema-2 writer emits the metal columns unconditionally, so a
metals-OFF output carries them with the value zero; a metals-ON restart from
such a file therefore ran with **zero metal density everywhere** while the base
boundary condition still counted metals in the mass and particle budget under
`eos_metals`. It was silent, and it produced confident nonsense: restarting the
converged `examples/15` (metals off) solution as a metals-on run reported
`||R|| = 9.54e-6` — it was reproducing the metals-off solution under a
metals-on label. An element is now taken from the file only if every stage has a
column AND the loaded stages are not identically zero everywhere; otherwise it
is built from the abundance exactly as `set_IC` does, `n_X = melem_ab * n_H` all
in the neutral stage, and the substitution is reported on stdout and in
`EXHALE_setup.out`. The rebuilt density is inserted before `calc_rho`, so the
restart keeps the mass closure `sum_i f_i A_i = 1` that a cold start has by
construction; for the same reason the legacy headerless-file path now carries
its metal mass in `rho` as well. *Measured* on that restart: the three elements
are reported, the base density at `j = 1` goes from 2.0618e14 to 2.0804e14
mH/cm3 (0.9 %, the metal mass), and the reported `||R||` goes from the false
9.54e-6 to 7.36e-2 — the honest distance of a metals-off wind from the metals-on
steady state. A metals-off restart, and a restart whose file does carry metal
densities, both give a byte-identical residual profile to before.

**Gates.** `make check` is 5/5 byte-identical (no regression case reaches the
JFNK hand-off and none uses `Load IC`, so neither change can touch them).
`examples/15` re-run from cold with the changed binary reproduces its reference
output **byte-identically** — `info = 0`, `||R|| = 9.542e-6`, 272 outer
iterations, 25837 marching steps, `log10 Mdot = 10.23`, H2 = H I front at
1.0084 R_p. The hot-Uranus molecular + solar-metals configuration
(`backup/regression/mol_metals` run to the Newton finish instead of its
12000-step regression snapshot) converges: `info = 0`, `||R|| = 3.342e-6` in 23
outer iterations, 34537 marching steps, `log10 Mdot = 10.30`.
`photo_deep_secion_cont` converges: `info = 0`, `||R|| = 7.284e-4` in 4
iterations, `log10 Mdot = 9.69`.

And `examples/16` itself now has a Newton-grade solution on the warm path.
Restarted from the converged `examples/15` state — which is exactly the restart
the `load_IC` fix repairs — with `Resid tol: 5.0e-5` it converges: `info = 0`,
`||R|| = 3.337e-5` in 30 outer iterations, `log10 Mdot = 9.65`. At 1e-5 it does
not, and the cold path converges at neither tolerance — 112920 marching steps to
the hand-off against 25837 for the metals-off twin, then 281 outer iterations to
a best `||R||` of 3.623e-4 with the worst scaled residual at `r = 1.021`, inside
the stagnant shell. Both remain open, in `docs/hd209_metal_stagnation.md`.

Files: `src/modules/time_step/steady_newton.f90`,
`src/modules/files_IO/load_IC.f90`,
`src/modules/files_IO/write_setup_report.f90`. Documentation:
`docs/hd209_metal_stagnation.md` (new), `TO_BE_DONE.md` item (A).

---

## 59. The molecular equilibrium no longer inherits a cell's previous composition (2026-08-13)

`ioniz_eq` solves the ionization equilibrium — and, under
`Molecular chemistry: True`, the dissociation equilibrium with it — cell by
cell, and the composition it returns is what the cooling, the photoionization
columns and the steady residual are all built on. When the coupled molecular
network produced no root inside the physical simplex in a cell, it used to hand
that cell its own INCOMING fractions back and warn (`molecular equilibrium
failed at N cells (kept previous state)`).

That is wrong on its own terms. The composition returned is then not the
equilibrium of the cell at the state being evaluated; it is whatever the
previous evaluation happened to leave there, and the error is not small. In the
JFNK finish it also means that the residual changes definition between outer
iterations: `f_sp` is overwritten by each accepted trial, so the map the next
iteration differentiates and line-searches differs from the previous one, by a
finite amount, at exactly those cells (`docs/hd209_metal_stagnation.md`
section 3, where the same section corrects an earlier claim — within one outer
iteration the Jacobian probes and the line-search trials do all start from the
same composition).

**What the cells were failing on.** Measured with a counting build on the
HD 209458 b molecular + solar C/N/O warm restart (40 marching steps,
single-threaded): 224 cell solves ended below the solver tolerance, 223 of them
keeping a physical but unconverged root and 1 inheriting its previous
composition. The rejected roots were not on another branch. At step 1, cell
`j = 117` (`r = 1.0294`, `T = 1652 K`), the first attempt converged
(`info = 1`) to a root violating the simplex only by `x(He II) = -3.8e-6` and
`x(H2) = -4.5e-4`, and the molecular-basin retry returned a strongly molecular
root, `x(H2) = 0.998`, violating it only by `x(He II) = -1.2e-8` — round-off
around a fully neutral helium. What fails is a root resolved onto a FACE of the
simplex, not an unphysical branch.

**Two changes, both in `ionization_equilibrium.f90`.**

1. *The molecular-basin retry seed is now a state.* It used to combine the
   chemical-equilibrium H2 fraction at the local `(p, T)` with the **previous
   state's** H II, He II, He III and 2^3S fractions — a pair that can put more
   H nuclei into H2 and H II together than the cell has, i.e. a starting point
   that is not a composition at all. The new seed,
   `dissociation_ionization_balance_at_fixed_ne`, is the limit of the network
   in which the couplings to the molecular ions are dropped: the H nuclei split
   between H2 and atomic H by the same Koskinen et al. (2022) Eq. 11 mixing
   ratio the molecular base boundary condition uses, the atomic remainder and
   every other element in their own ionization balance at the incoming `n_e`,
   H2+/H3+/HeH+ at zero. It is inside the simplex by construction and depends
   only on the local state.

2. *A cell whose roots all left the simplex keeps the closest one, clamped.*
   `element_budget_violation` measures how far a root lies outside the allowed
   states (the largest negative stage fraction, and the largest amount by which
   one element's tracked stages exceed its nuclei); the least-offending root is
   kept and `clamp_fractions_to_element_budget` moves it onto the boundary —
   negative fractions to zero, and any element still over its nuclei scaled to
   sum exactly to them. Given the measurement above, this is the treatment of a
   root the solver has resolved onto a face and delivered as a small number of
   either sign: the clamped state is that root, to the accuracy the solve
   reached. The cell's previous composition is never used.

The warning is now `every molecular equilibrium root left the physical simplex
at N cell(s), clamped onto the element budget`, and the run-wide count is
reported at the end of the run next to the `ioniz-eq roots` line. A cell that
reaches the solver tolerance on a physical root takes exactly the path it took
before, so the change is inert wherever the solve is healthy.

**One thing was tried and reverted.** `ionization_fractions_physical` rejects a
root for a component below `-1e-10`, and the rejected roots above sit at
`-1e-8`, inside the solver's own `xtol = sqrt(machine epsilon)`. Widening the
test to that band is defensible on paper and worse in practice: it lets a cell
stop on a root just outside a face instead of retrying from another starting
point, and the retried root is the better one. With the band widened,
`examples/15` goes from `info = 0` at `||R|| = 6.105e-6` in 89 outer iterations
to `info = 2` at 1.112e-4 in 133. The band stays at `1e-10` and the measurement
is recorded at the site.

**Gates.** `make check` is 4/5 byte-identical. The fifth, `mol_metals`, moves:
one cell in the 12000-step run takes a retry it did not take before (`ioniz-eq
roots` restarts 9 -> 10), and the difference that survives to the output is at
round-off — at most `1.3e-11` relative in every `Hydro_ioniz.txt` column and in
every H, He and metal density, `1.1e-7` in H2 and H2+, and `2.3e-4` in H3+ in
the far wind where H3+ is `1e-17 cm^-3` against a peak of 52. `log10 Mdot` is
10.58 before and after. The golden was re-snapshotted for that case only, with
`check` -> `golden` -> `check`. No molecular regression case reaches the clamp
(`0 cell(s)` in all three).

`examples/15` (HD 209458 b, molecular, cold) gets *better*, not merely
different: `info = 0` at `||R|| = 6.105e-6` in **89** outer iterations against
`info = 0` at 9.542e-6 in 272, hand-off at marching step 25839 against 25837,
`log10 Mdot = 10.24` against 10.23, H2 = H I front at 1.0098 R_p against
1.0084. The hot-Uranus molecular + solar-metals configuration holds: `info = 0`
at `||R|| = 5.826e-6` in 23 outer iterations against 3.342e-6 in 23, same 34537
marching steps, `log10 Mdot = 10.30` unchanged. Runtime checks
(`-O1 -fcheck=bounds,do,mem`) are clean on `mol_metals`, `mol_base_handoff` and
on the warm restart that does exercise the clamp.

`examples/16` (HD 209458 b, molecular AND solar C/N/O), the case the memo is
about, does **not** improve — it gets worse. The warm restart that reached
`info = 0` at `||R|| = 3.337e-5` in 30 outer iterations now bottoms at 7.931e-5
in 67 and reports `info = 2`, identically at `Resid tol` 5e-5, 2e-5 and 1e-5;
the cold path bottoms at 1.081e-3 against 5.038e-4 before. Each half of the
change on its own converges that restart — retry seed alone at 3.059e-5, clamp
alone at 4.853e-5 — and each half on its own fails `examples/15`, which the two
together converge in a third of the iterations it used to take. Both
HD 209458 b cases are
marginal solves and an O(1) change to a handful of cells' composition moves
them either way; what does not move is that their worst scaled residual sits in
the stagnant shell at `r = 1.021` throughout
(`docs/hd209_metal_stagnation.md` section 6). Removing the inheritance is
justified on its own terms, not by what it does to these two solves.

Files: `src/modules/radiation/ionization_equilibrium.f90`,
`src/EXHALE_main.f90`. Documentation: `docs/hd209_metal_stagnation.md`
sections 3 and 5, `TO_BE_DONE.md` item (A).
## 60. Damping the stagnant layer: a gated fourth difference inside the numerical flux (`Low-Mach damping`) (2026-08-13)

The HLLC flux resolves the middle wave of the Riemann fan exactly. The jump that
wave carries — the entropy jump, at constant pressure and velocity — moves at
the contact speed `S* ~ v`, and the dissipation the solver applies to it is
proportional to `|S*|`. It therefore vanishes as the flow stagnates. The
acoustic families keep their damping at `|v| +- c_s`, so a stationary
cell-to-cell entropy pattern is a discretely undamped mode of the scheme, and
in a near-hydrostatic layer the velocity field is tied to it by continuity and
by the force balance, so it carries the same pattern.

HD 209458 b with molecular chemistry and solar C/N/O is where that stops being
academic. The metal cooling exceeds the photoionization heating by a factor ~10
at `r ~ 1.02 R_p` and stops the flow: the mean `|v|` over `r = 1.015-1.030`
falls to 6-24 cm/s against ~207 cm/s in the metals-off twin (`M ~ 3e-5` against
a sound speed of ~4.5 km/s), and the alternating component of `v` grows to
4.7-10 times the local mean `|v|` where the metals-off solution carries 0.001.
Every failed steady solve of that configuration puts its worst scaled residual
in that shell (`docs/hd209_metal_stagnation.md` sections 1 and 6).

**The term.** `Low-Mach damping: <eps4> [<M_th>]` adds a fourth-difference
(Jameson, Schmidt & Turkel 1981) stress to the numerical **momentum** flux at
each interior face,

```
D_p{j+1/2} = eps4 g(M) rho_f lambda_f (v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1})
lambda_f   = 1/2 [ (|v| + c_s)_j + (|v| + c_s)_{j+1} ]
```

and, because a stress does work, `D_E{j+1/2} = v_f D_p{j+1/2}` to the **total
energy** flux — the same pairing `viscous_conduction` uses for the physical
Navier-Stokes stress. Without it the kinetic energy the stress removes would
come out of the internal energy instead of being converted into it. With it the
pair conserves total energy exactly, and the internal-energy change is
`-D_p dv/dr`, which has either sign cell by cell but integrates to
`+int mu_4 (d^2 v/dr^2)^2 dr >= 0`: kinetic energy into heat, never the reverse.
The undivided third difference is `dr^3 d^3v/dr^3 + O(dr^4)`, so the stress is
`tau = -mu_4 d^3v/dr^3` with `mu_4 = eps4 g rho c_s dr^3` — the acoustic
momentum diffusivity `rho c_s dr` an upwind scheme already carries, reapplied at
fourth order.

**The mass flux is not touched**, and neither is the energy flux other than
through the work term. The plain JST form, a fourth difference on the conserved
variables themselves, is inadmissible here and was measured to be so: the
physical mass flux `rho v` and energy flux `v(E+p)` both vanish with the
velocity while `dr^3 d^3(rho)/dr^3` and `dr^3 d^3(E)/dr^3` do not, because rho
and E are stratified over a barely resolved scale height. Their dissipative
divergences come out 1e4-1e5 times the physical mass-flux divergence and ~5e3
times the radiative source — they would rewrite the base velocity and thermal
structure rather than damp an oscillation. The momentum flux `rho v^2 + p` does
not vanish at stagnation (it tends to `p`), which is why the momentum form stays
small there.

**The gate.**

```
g = [ max(0, 1 - M_f^2/M_th^2) ]^2 ,   M_f^2 = 1/2 (M_j^2 + M_{j+1}^2)
```

is 1 at `M = 0`, exactly 0 for `M_f >= M_th`, and C^1 at the threshold so the
steady residual stays differentiable for the Newton solve. It is written in
`M^2` rather than `|v|/c_s` because `|v|` has a kink at `v = 0` and the stagnant
layer is exactly where `v` changes sign.

**Discrete consistency is the point of the design.** The stress is added inside
`RK_rhs`, the one routine the marching loop and `assemble_residual` — hence the
JFNK steady solver — both call. The equation the Newton residual measures is the
equation the marching loop relaxes, and the fixed point of one is the zero of
the other. This is exactly what the Shapiro filter lacks: it is applied to the
marching state only, so the Newton residual never sees it and the mode it
suppresses during marching is still an undamped mode of the system Newton
solves.

Both fluxes are set to zero at the base face (`j = 0`) and the outer face
(`j = N`), so nothing is injected or removed through the boundaries; those are
also the faces whose four-cell stencil would reach outside the ghost layer. Cell
`j` then depends on cells `j-2 .. j+2`, which is exactly the stencil the WENO3
residual already has (`WL(:,j)` reads cells `j-1..j+1`, `WR(:,j)` cells
`j..j+2`), so the banded Jacobian bandwidth `kl_jac = ku_jac = 8` is unchanged.

**Stability.** On the `2 dr` mode `v_j = (-1)^j a` the fourth difference is
`16 a`, so the mode decays at `16 eps4 lambda/dr`; the explicit step is
`dt = CFL dr/lambda`, so it is damped by `16 eps4 CFL` per step and explicit
stability requires `eps4 < 1/(16 CFL)` — 0.10 at the default `CFL = 0.6`.
`input_read` warns when the bound is violated. The classical JST range
`1/64`-`1/32` sits inside it.

**How big the term actually gets.** A numerical dissipation is admissible only
where it is negligible against the physical fluxes, and neither the fourth
difference nor the gate is a bound on that by itself, so any run that uses the
key now reports it: `contact_mode_dissipation_magnitude` prints the peak
`|D_p|/|rho v^2 + p|` over the interior faces, where it occurs, and the
outermost face at which the gate is still open. Measured on the converged
profiles (recomputing the same quantities in cgs from `Hydro_ioniz.txt`,
independently of the Fortran, agrees to the printed digits):

| run | `eps4` | peak `\|D_p\|/\|rho v^2+p\|` | at | gate open out to | first `M > 0.1` |
|---|---|---|---|---|---|
| `examples/15` (metals off) | 2e-2 | 1.83e-4 | 1.0005 | 1.043 R_p | 1.76 R_p |
| `examples/16` warm | 1e-2 | 8.36e-5 | 1.0005 | 1.108 R_p | 1.99 R_p |
| `examples/16` warm | 5e-3 | 3.94e-5 | 1.0005 | 1.105 R_p | 1.99 R_p |
| hot Uranus (molecular + solar metals) | 2e-2 | 3.70e-5 | 1.0005 | 1.083 R_p | — |

The peak is at the first interior face, where the base boundary condition forces
the steepest velocity gradient in the domain, and it is 0.02% or less of the
physical momentum flux there; through the stagnant shell itself it is `1e-7` or
below. The first cell with `M > 0.1` is at 1.76-2.00 R_p, so the wind, the
`r_esc = 2 R_p` window over which `||R||` is measured and the `N-20` cell at
which `Mdot` is evaluated all sit where the term is identically zero.

**Gates.** `make check` is **5/5 byte-identical** (`wasp_full`,
`wasp_he23off`, `mol_base_handoff`, `mol_metals`, `mol_lyman_werner`); no golden
was touched, and none should have been, since none of those inputs carries the
key. A `-O1 -fcheck=bounds,do,mem` build is clean on two key-ON cases (the
HD 209458 b molecular + C/N/O warm restart, which exercises the open gate, and
the hot-Uranus cold start, which exercises the closed one).

**What it does to the case it was written for.** HD 209458 b with molecular
chemistry AND solar C/N/O, warm-started from the converged `examples/15` state,
`Resid tol: 1.0e-5`, `Solver: Newton 5.0e-2`, `Max steps: 150000`, same binary
and same restart files in every row:

| `eps4` | hand-off step | `info` | best `\|\|R\|\|` | outer it | `log10 Mdot` | H2 = H I front |
|---|---|---|---|---|---|---|
| off | 48060 | 1 | 5.216e-5 | 500 (cap) | — | — |
| 5.0e-3 | 48062 | **0** | **9.934e-6** | 162 | 9.65 | 1.0112 |
| 1.0e-2 | 48007 | **0** | **9.987e-6** | 202 | 9.67 | 1.0117 |
| 2.0e-2 | 65836 | 2 | 2.529e-3 | 237 | — | — |
| 4.0e-2 | 66015 | 2 | 2.870e-3 | 66 | — | — |

This is the first `info = 0` at `Resid tol: 1.0e-5` this configuration has
produced. The two converged rows stop at 9.9e-6 because the solver tests the
tolerance at the top of each outer iteration and exits on the first iterate
below it, so those are first crossings and not a knife edge; the key-off
control, from the same restart, spent its whole 500-iteration budget and
plateaued at 5.2e-5. Nothing diverges in the rows that do not converge: each
falls back to time-marching and runs out the 150000-step cap at a `du` plateau
of 1.0e-3 to 2.6e-3, `log10 Mdot` 9.58-9.66. Those states are marching states,
not Newton-grade ones, which is why the table leaves their `Mdot` blank.

**The coefficient has a working range, and it is the bottom of the classical
one.** At `eps4 >= 2e-2` this solve is *worse*, and the reason shows up before
the Newton phase: the stress changes the base solution enough to slow the `du`
descent (`du` at step 40000 is 0.139 with the key off and 0.360 at
`eps4 = 2e-2`), so the hand-off comes 18000 steps later from a different state.
The key is opt-in and the pair has to be run on any new configuration; the
number 2e-2 is not a recommended default and there is no default.

**Sensitivity, honestly.** Halving `eps4` from 1e-2 to 5e-3 moves `log10 Mdot`
by 0.02 and the H2 = H I front by 0.0005 R_p, and in the wind (`r > 1.2 R_p`)
`rho` by <= 7.4%, `v` by <= 3.7% and `T` by <= 0.9%. Inside the stagnant layer
the two converged states differ a great deal (`rho` by 137% at `r = 1.013`, `T`
by 38% at `r = 1.015`) and differ in what the layer looks like: the alternating
component of `v` over `r = 1.015-1.030` is 0.024 at `eps4 = 1e-2` and 4.29 at
5e-3. Both satisfy `||R|| < 1e-5` on the wind window, which is the point — the
convergence measure does not look at that layer, so it cannot pin it down.

**The configurations that already converge do not move.** Cold starts,
`Resid tol: 1.0e-5`, `eps4 = 2e-2` where on:

| case | key | hand-off | `info` | `\|\|R\|\|` | outer it | `log10 Mdot` | front |
|---|---|---|---|---|---|---|---|
| `examples/15` (HD 209458 b, molecular) | off | 25839 | 0 | 6.105e-6 | 89 | 10.24 | 1.0098 |
| `examples/15` | on | 25839 | 0 | 3.642e-6 | 100 | 10.24 | 1.0095 |
| hot Uranus (molecular + solar metals) | off | 34537 | 0 | 5.826e-6 | 23 | 10.30 | 1.0774 |
| hot Uranus | on | 34537 | 0 | 6.252e-6 | 23 | 10.30 | 1.0774 |

Both hand off at the same marching step and reach the same `Mdot`; on the hot
Uranus the entire profile moves by at most 7.3e-4 relative in `rho` and 7.1e-4
in `T`. On `examples/15` the wind above 1.2 R_p moves by at most 1.2% in `rho`,
0.5% in `v` and 0.16% in `T`; the first cells above the base move more (`rho` up
to 16%, `T` 6.9% below 1.05 R_p), and those are cells that carry an odd-even
velocity oscillation 86 times the local mean with the key *off*. The cold
`examples/16` start (`Resid tol: 5.0e-5`) improves and still fails: `info = 1`
at 9.307e-4 with `eps4 = 2e-2` against 1.081e-3 with the key off, from the
identical hand-off step 112920 — the gate is closed through the cold relaxation,
when the flow is still fast everywhere, and the term only acts at the end.

Files: `src/modules/flux/low_mach_dissipation.f90` (new),
`src/modules/time_step/RK_rhs.f90`, `src/modules/init/parameters.f90`,
`src/modules/files_IO/input_read.f90`,
`src/modules/files_IO/write_setup_report.f90`, `src/EXHALE_main.f90`,
`Makefile`. Documentation: `docs/hd209_metal_stagnation.md` sections 7 and 8,
`docs/EXHALE_user_manual.tex`, `docs/input_schema.md` K31b, `README.md`,
`TO_BE_DONE.md` item (A).

## 61. The residual below the escape radius, and what closes the HD 209458 b molecular + metals case (2026-08-15)

Section 60 left the `examples/16` configuration converging from a warm restart
at two damping coefficients and failing from a cold one, and left open whether
the acceptance window should be widened to include the layer that the failures
live in. Both questions, and the seed the warm prescription depends on, were
measured on 2026-08-15. The full account is
`docs/hd209_metal_stagnation.md` section 9; this section records the code
change and the conclusions.

### The change: the layer's residual is printed, not tested

`resid_relnorm_below_escape` and `write_resid_below_escape` in
`src/modules/time_step/steady_newton.f90` apply the **same** volume-weighted
relative norm that `resid_relnorm` uses for the convergence measure to the
complementary region `[1:j_min-1]`, i.e. below the escape radius, and locate
the cell carrying the largest cell-wise relative residual
`|F_k,j| / |u_k,j|` there, reporting `j`, `r(j)` and which of the three rows it
is. One line is printed at four points: PTC start, PTC end, JFNK start, JFNK
end. The JFNK end line is emitted *before* the best-iterate restore, because
that is the last point at which `F` and `u` are a consistent pair; it therefore
refers to the last iterate visited, which is also the returned state whenever
no restore happens.

The value is never tested and nothing in the solve depends on it. The whole
change is **+95 lines in that one file**, output only, and `make check` is
**5/5 byte-identical** (`wasp_full`, `wasp_he23off`, `mol_base_handoff`,
`mol_metals`, `mol_lyman_werner`) — no output file moves, only stdout grows.

### Widening the acceptance window is rejected by what it prints

With the diagnostic in place the question is a measurement rather than a
design argument. The norm over `[1:401]` (`j_min = 402`) on states that
converge today:

| state | norm below `r_esc` | worst cell |
|---|---|---|
| `examples/16` converged, `eps4 = 5e-3` | 0.195 | `j = 75`, `r = 1.0153` |
| `examples/16` converged, `eps4 = 1e-2` | 2.352 | `j = 5`, `r = 1.0010` |
| `examples/16` converged, `eps4 = 5e-3` + LW | 0.939 | `j = 3` |
| `examples/15` converged (**metals off**) | 9.99 | `j = 2` |
| warm restart at JFNK start | 0.60-1.21 | `j = 6`-15, momentum |
| cold `eps4 = 5e-3`, at failure | 594 | `j = 97`, `r = 1.0218` |
| cold `eps4 = 1e-2`, at failure | 1599 | `j = 112`, `r = 1.0273` |

A converged state carries a relative residual of 0.2 to 10 below the escape
radius, four to six orders above `Resid tol`, and the largest of them belongs
to the metals-off `examples/15` solution that nobody doubts. That is not about
metals: the layer is near-hydrostatic (`|v| ~ 10 cm/s`), so the denominator of
the momentum row's relative residual, `|rho v|`, collapses wherever the flow is
slow, right state or not. Widening the window on this norm would report
`info = 2` for every configuration that converges today without determining the
layer any better. The window stays `[j_min:N]`; the diagnostic makes the
uncontrolled part visible instead. The failure rows sit three orders above the
converged ones, so the number does separate good from bad — it is the absolute
scale that is not usable as a tolerance.

### The cold path is closed; the warm prescription is pinned to its seed

Scanning `Low-Mach damping` on the cold start gives `info = 1` at 1.081e-3
(off), 9.953e-4 (5e-3), 1.025e-3 (1e-2) and 9.307e-4 (2e-2) — one ~1e-3 floor,
a factor 100 above the tolerance, and the two coefficients that work on the
warm restart do nothing here that `off` does not. The cold hand-off state
appears to lie outside the Newton basin the warm restart reaches, and the
coefficient is not the lever. The cold path for this configuration is
abandoned.

The warm prescription needs one thing section 60 did not say: the restart must
come from an `examples/15` state converged with the **current** binary. Two
seeds, same binary and same keys otherwise:

| seed | `eps4` | hand-off | `info` | best `\|\|R\|\|` | `log10 Mdot` |
|---|---|---|---|---|---|
| pre-section-59 `examples/15` output (kept as `output_pre_item3_20260813`) | 5.0e-3 | 50645 | **2** | 1.303e-4 | — |
| the same | 1.0e-2 | 50629 | **2** | 5.770e-5 | — |
| the same, `+ Stellar LW flux: 343.0` | 5.0e-3 | 50189 | 0 | 9.799e-6 | 9.68 |
| `examples/15` re-converged with the current binary | 5.0e-3 | 48062 | 0 | 9.934e-6 | 9.65 |
| the same | 1.0e-2 | 48007 | 0 | 9.987e-6 | 9.67 |
| the same, `+ Stellar LW flux: 343.0` | 5.0e-3 | 47965 | 0 | 9.900e-6 | 9.68 |

The re-converged seed reproduces the section 60 table exactly, hand-off step
included (that re-convergence itself reproduces the section 60 control:
step 25839, `info = 0`, `||R|| = 6.105e-6`, `log10 Mdot = 10.24`, front at
1.0097 R_p). The stored seed does not: the hand-off comes ~2600 steps later,
from a different state, and the line search then collapses (`lambda ~ 1e-6`)
and aborts. The section 59 change moved the `examples/15` solution by more than
this marginal solve tolerates in its starting point, so a restart file written
by an older binary is no longer the same seed. The Lyman-Werner member
converges from both seeds — the LW coupling appears to widen the basin here
rather than narrow it.

### What the weak determination of the layer costs

The two converged solutions (`eps4 = 5e-3` and `1e-2`, both `info = 0` at
`Resid tol: 1.0e-5`) agree in the wind — `Drho/rho` +5.1% at 1.3 R_p, +3.4% at
1.5, +2.6% at 2.0, `Dv/v <= 1.6%`, `DT/T <= 0.4%`, `log10 Mdot` 9.65 against
9.67 — and disagree inside the layer, where the `1e-2`/`5e-3` density ratio
runs 0.78 at `r = 1.005`, 1.31 at 1.013, 1.41 at 1.020, 1.32 at 1.030, `T`
differs by +20% at 1.005, +38% at 1.015 and -7.7% at 1.020, and the alternating
component of `v` over `r = 1.015-1.030` is 0.392 against 0.028 on a mean `|v|`
of 12-16 cm/s.

That does not stay internal. Run through `EXHALE_transit.py`, the theoretical
peak excess absorption is 40.19% against 40.26% in He I 10830 (0.2% relative),
2.098% against 2.199% in H-alpha (+4.8%), 1.046% against 1.134% in H-beta
(+8.4%), with Ly-alpha saturated in both. He I 10830 forms in the wind and does
not care which layer sits underneath it; part of the H(n=2) absorption is
formed in the layer, so the weak determination appears to propagate to the
Balmer depths at the 5-8% relative level. Quote those depths with the pair, not
from a single coefficient.

### Lyman-Werner with metals on

The `Stellar LW flux [erg/cm2/s]: 343.0` A/B on the re-converged seed at
`eps4 = 5e-3`, both members `info = 0`:

| | LW off | LW on (343) |
|---|---|---|
| `log10 Mdot` | 9.65 | 9.68 |
| H2 = H I front | 1.01108 | 1.01026 |
| `T` at `r = 1.02` [K] | 1801 | 1470 |
| `x(H2)` at `r = 1.02` | 4.0e-6 | 2.1e-5 |
| H3+ peak [cm^-3] | 1.93e5 at 1.0006 | 1.99e5 at 1.0006 |
| HeH+ peak [cm^-3] | 0.285 | 0.053 |
| C II fraction at 1.02 / 1.20 | 0.197 / 0.449 | 0.180 / 0.430 |
| H-alpha peak | 2.098% | 2.297% |
| H-beta peak | 1.046% | 1.221% |
| He I 10830 peak | 40.19% | 40.33% |

The self-shielding is the section 56 physics unchanged: at the base
`N(H2) = 2.68e21 cm^-2`, `f_shield = 2.11e-6`, `k_LW = 1.27e-10 s^-1`; above
`r = 1.10 R_p` the shielding is gone (`f_shield = 1.000`,
`k_LW = 6.026e-5 s^-1`). `x(H2)` at `r = 1.02` is *higher* with the band on,
which appears to be the layer temperature rather than the chemistry: that layer
is 330 K colder, and the temperature dependence of the H2 balance outweighs the
added destruction. The front moves down by 0.0008 R_p, consistent with more
destruction where the shielding has lifted. Caveat: the band changes the
H-alpha peak by +9.5% relative and H-beta by +17%, only about twice the
numerical spread above, so those changes are not quotable as observational
predictions without the coefficient pair run alongside.

Files: `src/modules/time_step/steady_newton.f90` (diagnostic only).
Documentation: `docs/hd209_metal_stagnation.md` section 9,
`examples/README.md`, `README_HOWTO.md`, `TO_BE_DONE.md` item (A).

## 62. Runtime grid size, CODATA k_B, the layer residual on a physical scale, and two diagnostic closures (2026-08-15)

One working series, no commit boundaries implied by the ordering below. The
verification was staged so that every piece except the k_B value is proven a
no-op on the regression matrix, and the golden refresh at the end carries the
k_B change alone.

### 62.1 `Grid cells:` — the cell count is a runtime quantity

`N` was `integer, parameter :: N = 500` in `global_parameters`; it is now a
runtime value set by the optional input key

```
Grid cells: 750
```

with default 500, the value that used to be compiled in. The grid-sized
module arrays of `global_parameters`, `lya_rt`, `excited_hydrogen` and
`ionization_equilibrium` are allocated (lower bound `1-Ng` unchanged) once
input parsing has resolved the key; the main program's state vectors likewise
(they were implicitly static and zero-started, and are zeroed on allocation to
reproduce that). The many grid-sized arrays inside procedures were already
automatic and needed nothing. `PH_heat_H`'s zero-filled `parameter` arrays
became locals assigned on entry — the one place a named constant had to give
way; an A/B on an H-only run (the path `make check` does not cover) is
byte-identical. `write_setup_report` echoes the resolved `N`, the parse dump
gains an `N` line, and `load_IC` now counts the IC file's records first, so a
restart against an IC written at a different `N` reports the mismatch instead
of dying on a bare end-of-file.

This closes TO_BE_DONE item (E): the base-refinement pairs
(`Base grid [dr,cells]`, `Grid cells`) = (`1.0e-4 100`, 500),
(`5.0e-5 200`, 658), (`2.5e-5 400`, 916), (`1.25e-5 800`, 1375) hold the upper
stretch at its default-grid value from `input.inp` alone, no rebuild.

Gates: `make` and `make wind_ae_ic` build (the standalone source set does not
include `parameters.f90`); `make check` **5/5 byte-identical** against the
pre-change goldens; a 750-cell tutorial run writes 754 data rows and finishes;
`-fcheck=bounds,do,mem` clean at 500 and 750 cells.

### 62.2 `kb_erg` is the CODATA value; `lyman_werner` local copy synced

`kb_erg = 1.38e-16` (the truncated ATES literal) is now the SI-exact
`1.380649e-16`. The judgment first: the old value was wrong at 4.7e-4
relative, and every thermal quantity in the code passes through it; that it
"only" moves results at that level is not an argument for keeping it. The
`kb_lw` copy in `lyman_werner.f90` (kept dependency-free on purpose) is synced
and its comment now says the mirror is mandatory.

Measured on the regression matrix (single-threaded, against the old goldens):
the two WASP-121b cases re-converge at count 14065/14040 (from 14060/14037),
field-level changes are a few 1e-4 to 3e-3 relative (largest in `cool`), and
`log10 Mdot` moves 13.22 -> 13.23, i.e. the last printed digit. The three
pinned molecular snapshots move at the same order. **The goldens were
re-snapshotted at the end of this series and carry exactly this change**; the
parse corpus was re-snapshotted for the new `N` line as well.

Still deliberately NOT current CODATA, recorded here so the list is explicit:
`kb_eV = 8.6167e-5` (CODATA 8.617333e-5; also internally inconsistent with
`kb_erg*erg2eV` at 5e-4 both before and after), `hp_erg = 6.62620e-27`,
`hp_eV = 4.1357e-15`, `mu = 1.673e-24` (a hydrogen mass, not m_p), `Gc`, and
the astronomical constants. Scope was pinned to `kb_erg`; a full constants
pass would be its own measured rebaseline.

### 62.3 The below-escape momentum residual has a physical scale

`resid_relnorm_below_escape` (section 61) divided every row by `|u|`; for the
momentum row that is `|rho v|`, which collapses in the quasi-hydrostatic
layer and made converged states read 0.2-10. The momentum row is now divided
by the gravitational force density `|rho| |Gphi_i(j)-Gphi_i(j-1)|/dr_j` (the
discrete potential difference of the source term), so the number is the
fractional violation of hydrostatic balance; the print reports the three rows
separately. Print-only; `make check` 5/5 byte-identical.

Re-measured on the section 9.4 states of `docs/hd209_metal_stagnation.md`
(full table and two side observations in its new section 10): the converged
`examples/15` (old scale 9.99), `examples/16` eps4=5e-3 (0.195) and 1e-2
(2.352) states all hold the layer at **1e-5 to 1.4e-4 of rho g**. The layer
is not momentum-unbalanced; what stays uncontrolled between the pair members
is which hydrostatic stratification it settles on.

### 62.4 `Heating_breakdown.txt` is consistent for a molecular run

`write_heat_breakdown_eq` reconstructed the atomic heating only and printed a
NOTE admitting H2 photoheating, He(2^3S)+H2 Penning and Lyman-Werner heating
were missing; its electron density also omitted the molecular ions and its
`xion` the molecular nuclei. It now rebuilds all of them exactly as
`ionization_equilibrium` does (same `isp_*` mapping as the cooling dump, same
`calc_ne`, same rate expressions), passes the H2 density into `PH_heat_HHe` so
the H2 column (col9) and the opacity are the solver's, and appends
`col15 heat_He23S_H2_Penning` and `col16 heat_H2_LW`. Measured on the
`mol_lyman_werner` case: max |sum(cols 5-16)/heat_total - 1| = 2.0e-15, with
all three molecular channels nonzero. Atomic runs are unchanged (the new
terms are exactly zero there); the file is not part of the regression
comparison. The `paper_data.py` fallback column list is synced.

### 62.5 `Lya absorbing bottom:` — a downward loss channel for the escape-probability closure (default off)

The section 30 audit follow-up: Huang et al. (2017) close their Ly-alpha
Monte Carlo with a pure-absorber bottom (H2 accidental-resonance true
absorption at N_H2 ~ 1e14 cm^-2), while our local closure had no lower
boundary, so Jbar near the base was likely overestimated. With

```
Lya absorbing bottom: True
```

a third escape channel joins the wing and Sobolev ones: the same Neufeld form
on the line-center depth from the cell down to the bottom of the domain,
`beta_bot = min(1, pi^-1/4 sqrt(a/tau_bot))`, `tau_bot = tau(bottom)-tau(r)`,
combined as parallel channels. It feeds the existing `beta_tot`, so the loss
acts consistently on the internal field and on the trapped stellar beam.
Default off, and the off path executes the original expression unchanged —
covered by the byte-identical regression pass below.

First look (HD 189733 b benchmark configuration, 300-step snapshot, so
qualitative only): Jbar and n_2p drop together by the same factor — 0 at the
innermost ghost, x0.24 at r = 1.0002, x0.74 at 1.01, back to within 1.5% by
r = 1.023 and exactly 1 above 1.05. The loss is confined to the bottom few
scale heights, as the audit expected; the Balmer-forming region above is
untouched. A converged A/B is the follow-up if a quantitative statement is
ever needed.

### 62.6 windae du-trigger arming, measured

The 2026-08-11 descending-crossing guard (newton_scaling_and_base_wall.md
section 9) had been measured on the transonic and cold-hydrostatic ICs but
never on `IC mode: windae`, whose ~3-step false du stop originally forced the
"always Newton-finish" rule. Measured on `examples/12_windae_ic_hd209`
(2000-step cap): du at step 2 is 1.43e-4 — below every threshold on the
untouched generated IC — and nothing fires; the du stop arms at step 23
(du = 1.027e-3) and the Newton hand-off at step 222 (du = 1.002e-2), both on
the first ascending crossing. The false stop is structurally gone for this IC
mode too; a du-stop Mdot still carries the usual path spread, so the Newton
finish stays the prescription. Addendum recorded in that memo's section 9.

### 62.7 Verification ledger for the series

| stage | tree | `make check` |
|---|---|---|
| runtime `N` (62.1) | 62.1 only | 5/5 byte-identical |
| + layer scale (62.3) + `load_IC` guard | 62.1+62.3 | 5/5 byte-identical |
| + breakdown fix (62.4) + Lya bottom off (62.5), k_B still old | all but 62.2 | 5/5 byte-identical |
| + `kb_erg` CODATA (62.2) | full series | all 10 comparisons FAIL vs old goldens, differences as in 62.2; goldens re-snapshotted; final `make check` 5/5 against the new goldens |

Files: `parameters.f90`, `input_read.f90`, `write_setup_report.f90`,
`load_IC.f90`, `EXHALE_main.f90`, `excited_hydrogen.f90`, `lya_rt.f90`,
`ionization_equilibrium.f90`, `util_ion_eq.f90`, `steady_newton.f90`,
`lyman_werner.f90`, `python/paper_data.py`. Documentation:
`docs/hd209_metal_stagnation.md` section 10,
`docs/newton_scaling_and_base_wall.md` section 9 addendum,
`docs/EXHALE_user_manual.tex`, `README.md`, `README_HOWTO.md`,
`TO_BE_DONE.md` item (E).

## 63. The rest of the physical constants, and why iron cannot move to Badnell (2026-08-19)

Two items that had been sitting in the "recorded, blocks nothing" list, taken
together because the first needs a golden refresh and the second turned out to
need none.

### 63.1 The constants pass

Section 62.2 updated `kb_erg` to the CODATA value and listed the constants
deliberately left alone. They are now done. Measured against CODATA 2018 and
IAU 2015:

| constant | was | now | relative change |
|---|---|---|---|
| `kb_eV` | 8.6167e-5 | 8.617333262e-5 | 7.4e-5 |
| `mu` | 1.673e-24 | 1.67353284e-24 | 3.2e-4 |
| `Gc` | 6.67259e-8 | 6.67430e-8 | 2.6e-4 |
| `hp_erg` | 6.62620e-27 | 6.62607015e-27 | 2.0e-5 |
| `hp_eV` | 4.1357e-15 | 4.135667696e-15 | 7.8e-6 |

`pi`, `erg2eV`, `c_light`, `parsec`, `AU` and `RJ` were already exact to the
printed digits and are unchanged. `mu` is the mass of the hydrogen ATOM
(m_p + m_e - 13.6 eV/c^2), which is what the density normalization means; it
is not the proton mass, and the comment now says so.

The same values appear in four other places, all synced: `m_h2` in
`lyman_werner.f90` (2 mu), `Gcgs` in `lower_column.f90`, and the Wind-AE
constants `wae_G`/`wae_K`/`wae_MH` (`wae_params.f90`), `wae_K` (`wae_soe.f90`)
and the hydrogen atomic mass in `wae_exhale_input.f90`. The Wind-AE files are a
port of an external code, so the original `defs.h` values are written into the
comment there: restoring `G = 6.67259d-8`, `K = 1.380658d-16`,
`MH = 1.6733d-24` reproduces the bit-exact port gates in
`backup/regression/windae_oracle` (working copy only -- `backup/` is not in the
git remote). Against those three the shifts are 2.6e-4, 6.5e-6 and 1.4e-4, so
the port now agrees with the reference at that level instead of bit-exactly;
`docs/wind_ae_solver.tex` says so where it lists the gates. The helium atomic
mass in `wae_exhale_input.f90` was left at Wind-AE's 6.6464790722e-24, which is
the CODATA He atom mass to 1e-6.

**Measured effect.** On the regression matrix the fields move by 1e-3 to 5e-3
relative and `log10 Mdot` does not move in the printed digits (`wasp_full`
13.23, counts 14065 -> 14059 and 14040 -> 14036). One outlier: in
`wasp_he23off` the photoheating of the single cell at r = 1.6056 -- the last
interior cell, next to the r = 1.610 boundary -- changes by 97%, from 2.617e-6
to 5.158e-6, because the old value was an isolated dip below its neighbour
(5.065e-6 at r = 1.6098) and the new one is not. Everything else in that file
agrees to 0.2%, median 8e-4. The goldens and the parse corpus were
re-snapshotted (the parse dump carries `n0`, `q0` and `b0`, which are derived
from `kb_erg` and `Gc`), and `make check` is 5/5 against them.

The four paper models were re-converged on the new constants: all `info = 0`,
and `log10 Mdot` is unchanged at 9.47 / 9.14 / 11.70 / 13.20. The LaRT fields
were NOT re-run: their inputs (n_HI, T, v) move by ~1e-3, an order of magnitude
below the Monte Carlo noise of a 2e6-photon run.

### 63.2 Iron recombination: Phase 6 of the Huang plan is closed, not deferred

`Cool_coeff.f90` carried a note that switching iron from the Huang et al.
(2023) Eqs. (5)-(6) analytic fits to Badnell RR+DR was "deferred to a later
pass", and the Huang reproduction plan listed that switch as its last
outstanding phase. It cannot be done, for a reason that has nothing to do with
this code: **the data do not exist.**

Fe I is produced by recombining Fe II, which is Mn-like (25 electrons), and
Fe II by Fe III, which is Cr-like (24). The Badnell dielectronic-recombination
project publishes by isoelectronic sequence, and as of its latest instalment
(paper XVI, the phosphorus sequence, 2022) it reaches 15 electrons. Neither
iron stage is covered, and neither is the K-like sequence that the Ca I
comment in the same file already records as missing. The other nine metals are
on Badnell RR+DR through `alpha_rec_metal`, which is what the phase asked for;
iron has one available source and uses it.

This is also what the plan itself prescribed --- its Phase 1c says to use
Badnell "falling back to Huang's explicit Fe fits Eqs (5)-(6) where Badnell
coverage is thin" --- so the code was following the plan rather than lagging
it. The misleading "deferred" comment is corrected at the code site and the
phase is closed in `docs/Huang_update_plan.md` with this reason. What remains
unbuilt is only the bookkeeping that phase's gate wanted: a runtime key
selecting one rate source for all metals, and a Huang-vs-Badnell delta table.
Neither applies to iron.

Files: `src/modules/init/parameters.f90`,
`src/modules/lower_atmosphere/{lyman_werner,lower_column}.f90`,
`src/modules/wind_ae/{wae_params,wae_soe,wae_exhale_input}.f90`,
`src/modules/radiation/Cool_coeff.f90` (comment only).
Documentation: `docs/Huang_update_plan.md` Phase 6, `TO_BE_DONE.md`.

### 63.3 Documentation pass on the md/tex corpus (2026-08-19)

The sweep that accompanied section 62 was run before 63.1/63.2 existed, so a
second pass propagated them and checked the reference documents against the code
again. What it found was not only the expected propagation:

- **The manual and the README claimed the iron swap was pending.** Both now
  carry the isoelectronic reason, and the manual note names the Ca I parallel
  (K-like, 19 electrons). `docs/atomic_data_EXHALE_vs_MoCHII.{md,tex}` and the
  Phase 1c entries in this file and in `docs/Huang_update_plan.md` gained the
  same pointer, so a reader landing on the old "deferred to Phase 6" wording is
  told immediately that it is closed.
- **The README described the line trapping as covering "the two fine-structure
  coolants" solved in a "two-level" scheme.** It has been eight lines across
  C I, C II, N II and O I (`n_fsline = 8`) in exact ground-term statistical
  equilibrium since that set was expanded; `cooling_formulas.tex` and the manual
  were already correct, the README and the comparison memo were not. The two
  dated memos that legitimately say "two" (`resonance_line_trapping.md`,
  `EXHALE_BC_and_IC.tex`) now say so of their own date.
- **`EXHALE_BC_and_IC.tex` quoted `log10 Mdot = 9.31` for HD 209458 b as the
  value "under the present code".** That was the 2026-08-11 number; the
  production value is 9.47 at `||R|| = 1.8e-4`, and the sentence now carries
  both with their dates.
- **Four line-number citations in `atomic_data_EXHALE_vs_MoCHII.tex` pointed at
  the wrong lines** (`rr_badnell` is at 551, the memo said 506, and so on).
  Line numbers in prose go stale on every edit above them, so they were replaced
  by the symbol names.
- **The resolving powers the transit tool actually uses were undocumented**, and
  Ca II / Na I default to *whatever `Instr_res_Ha` is*, so overriding the
  H-alpha resolution moves the two metal doublets with it. `README_HOWTO.md` now
  tabulates the seven defaults and states the trap.
- **The Wind-AE port gates.** `docs/wind_ae_solver.tex` said the base-BC layer
  matched the C/Python reference bit-exactly. After 63.1 it matches at
  1e-6--3e-4; the paragraph now says which three lines to restore to recover the
  bit-exact comparison, and names the driver directory correctly
  (`backup/regression/windae_oracle`, working copy only) here and at the code
  site.

Two mechanical cross-checks were also run in both directions: all **56** input
keys parsed by `input_read.f90` appear in the manual, and every key-like entry
in the manual's tables is still parsed (the 15 apparent orphans are keys the
manual spells with their unit suffix, e.g. `Planet mass [M_J]`, which the parser
matches on the prefix). No key is missing or obsolete.

PDFs rebuilt and pages inspected: user manual 54, Update_EXHALE 118,
EXHALE_BC_and_IC 23, wind_ae_solver 7, atomic_data 7. Zero errors and zero
undefined references in all five; the 6 and 38 overfull boxes in the manual and
this log are pre-existing long table rows, none of them on an edited page.

### 63.4 What the constants pass had still missed, and the golden refresh it needed (2026-08-19)

Section 63.1 swept `.f90` files for the constants named in `parameters.f90`. That
was the wrong search: the same constants also sit in local copies with other
names, in Python tools that reproduce the same physics, and -- in one case --
inside a cooling coefficient. Rescanning for the *values* rather than the names
turned up the following.

**In the solution path (this is what needed the golden refresh).**
`Cool_coeff.f90` computed He II recombination cooling as `1.38e-16*T*alpha_B` in
two places (`rec_cool_HeII` and `rec_cool_HeII_func`). The literal is the
Boltzmann constant, not a fit coefficient: the `k_B T alpha_B` form *is* the
Hui & Gnedin (1997) He II expression, while H II and He III use their fitted
case-B rates. These two copies were missed when section 62.2 replaced the
truncated ATES literal `1.38e-16` with the CODATA `kb_erg` everywhere else, so
the coefficient was 4.7e-4 low. Now `kb_erg`.

**Report-only or generated-file paths.**

- `lower_column.f90`: `mh_g` was the proton mass while multiplying a
  dimensionless mean molecular weight, so the Tier-1 analytic scale height was
  5.6e-4 too small. It is the hydrogen atom now. This affects no golden:
  `lower_column_solve` is called report-only when `input.inp` carries
  `Lower column:`, and no `input.inp` in the repository does.
- `wae_ic_writer.f90`: `ICW_mH` likewise, the unit the Wind-AE IC is written in.
- `species_diffusion.f90`: `m_amu_g` was labeled "amu in g" but held the proton
  mass. The species masses in that file are on the H = 1 scale
  (`m_H_amu = 1`, `m_He_amu = 4`), so the multiplying mass is the hydrogen atom
  -- neither the proton nor the atomic mass unit u. Value and comment fixed.

**Python tools.** Fifteen constants across ten files, all now equal to their
`parameters.f90` counterparts: `src/utils/glob.py` (the GUI's mirror of the
global block -- `Gc`, `mu`, `kb`, `hp_eV` were all stale),
`exhale_transit_lib.py` (`G` was the three-digit 6.67e-11; it enters only the
tidally-locked rotation period, so the induced change is -1/2 dG/G = 3.2e-4 on
the broadening velocity), `examples/exhale_io.py` (the loader every notebook
uses), `EXHALE_plots.py`, `eta_approx.py`, `src/utils/run_lower.py`,
`src/utils/vulcan_to_base.py`, `src/utils/vulcan_driver.py`,
`src/utils/windae_to_exhale_ic.py`, and the two `docs/compare_*.py` scripts.
`cooling_data/{make_cooling_doc_figures,build_cno_notebook}.py` keep
`KB_CODE = 1.38e-16` deliberately -- that is the rounded value the CHIANTI fits
were made with, so the scripts must keep it to reproduce the published numbers;
only the comment claiming it is "the code-wide kb_erg" was corrected.

`run_lower.py` integrates the same column as `lower_column.f90`, so its
correction can be measured directly: `r_base` moves 1.47226 -> 1.47216 R_J
(HD 209458 b), 1.17353 -> 1.17350, 2.00977 -> 2.00953 and 1.43716 -> 1.43700 --
the fifth decimal, and the three-decimal table in
`examples/13_lower_atmosphere/README.md` does not change. The generated
`base_{iso,guillot}.inp` files were **not** regenerated: the examples' recorded
results were produced from them, and refreshing the inputs without re-running
the examples would leave the two inconsistent. The pinned
`backup/regression/mol_*/base.inp` are frozen with their goldens and are never
regenerated.

**Measured effect and the golden refresh.** Staged so the two kinds of change
are separable. First the report-only fixes alone: `make check` **5/5
byte-identical**, which is the measurement proving `mh_g` never reaches a
solution. Then the `kb_erg` fix in the He II cooling: all five cases FAIL, as
they must. Against the old goldens the fields move by

| case | max rel. change | median non-zero | step count |
|---|---|---|---|
| `wasp_full` | 3.7e-6 (`cool`), 7.6e-6 (`Ion_species`) | 3.7e-7 / 1.5e-7 | 14059, unchanged |
| `wasp_he23off` | 3.8e-6 (`cool`), 7.3e-6 | 3.7e-7 / 1.4e-7 | 14036, unchanged |
| `mol_base_handoff` | 1.5e-3 (`v`), 1.1e-5 | 1.6e-8 / 7.9e-9 | 12000 (pinned) |
| `mol_metals` | 9.6e-5 (`cool`), 2.0e-4 | 1.4e-8 / 3.8e-9 | 12000 (pinned) |
| `mol_lyman_werner` | 1.0e-4 (`cool`), 1.5e-5 | 1.6e-8 / 6.0e-9 | 12000 (pinned) |

The column that moves most is `cool`, which is the term that was edited. The
molecular cases have almost no He II, so their median change is two decades
smaller than the atomic ones; their larger *maxima* are single cells of a
relaxation snapshot pinned at 12000 steps, where a 1e-8 perturbation is still
being amplified by the transient. `log10 Mdot` is 13.23 for both `wasp` cases,
unchanged. Goldens re-snapshotted, `make check` back to 5/5.

The four paper models were **not** re-converged. The bound above -- 3.7e-6 on
the hottest, most strongly ionized case, which is the WASP-121 b-like one -- is
three decades below the printed precision of every number in `paper/`, and two
to three decades below the constants pass of section 63.1, which itself moved
exactly one printed transit digit. Re-running them would reproduce the same
printed values at a cost of hours.

Files: `src/modules/radiation/Cool_coeff.f90`,
`src/modules/lower_atmosphere/lower_column.f90`,
`src/modules/wind_ae/wae_ic_writer.f90`,
`src/modules/functions/species_diffusion.f90`, and the ten Python tools listed
above. Goldens re-snapshotted (`make check` 5/5); the parse corpus is unchanged
and was **not** re-snapshotted -- 125 cases OK -- because none of these constants
enters a derived setup value.


## 64. Resolved-configuration output, and the transit tool reading it (2026-08-22)

`EXHALE_setup.out` records the resolved configuration in prose, but
`EXHALE_transit.py` re-parsed `input.inp` only — so a `base.inp` handoff that
overrides the planet radius, base temperature, or He/H reached the wind and
not the transit geometry.  Measured on `benchmarks/wasp52`: the wind used
1.40492 R_J / 1181.2 K / He/H 0.0959 (from `base.inp`) while the transit read
1.27 R_J / 1304 K / 0.0204 from `input.inp`.  The radius enters chord
lengths, absorbing areas, and the rotation velocity, so this was a real
transit-side error for every handoff run.

Fix, per `docs/lhs1140b_lower_atmosphere_plan_new.md` Phase B:

- **`EXHALE_resolved.out`** (new, machine-readable): written next to
  `EXHALE_setup.out` at startup by `write_resolved_config`
  (`write_setup_report.f90`) after `input.inp` + `base.inp` resolution.
  Key-value lines: `planet_radius_RJ`, `planet_mass_MJ`,
  `equilibrium_temperature_K`, `HeH_number_ratio`, `orbital_distance_AU`,
  `star_mass_Msun`, `base_inp_present`.
- `exhale_transit_lib.read_input_params` prefers that file when it sits next
  to `input.inp` (returns `resolved=True`, plus the resolved `HeH`); absent
  file falls back to `input.inp` exactly as before.  A present-but-unreadable
  file is an error, not a silent fallback.  `EXHALE_transit.py` prints which
  source it used.
- **`he_line_metrics.py`** (new, repo root): the He 10830 line metrics of
  Cherubim et al. (2026) — three Gaussians with shared width and shift;
  blended-red depth, blue depth, red/blue amplitude ratio, FWHM of the
  blended feature — with a synthetic-profile selftest
  (`python3 he_line_metrics.py --selftest`).  `EXHALE_transit.py` now fits
  the instrument-convolved He curve (air frame) and writes
  `tpm_He10830_metrics.txt`; a fit failure warns and does not stop the run.

Verified: a scratch run with a `base.inp` overriding radius/temperature/He
shows the transit consuming 1.50 R_J / 900 K (the override), and the fallback
path still returns the `input.inp` values; `make check` byte-identity, since
the new file is additive and the compared outputs are untouched.

## 65. The metastable helium triplet is on by default (2026-08-22)

`Include He23S` was a mandatory key initialized to `.false.`, so every run
carried whatever the file it was copied from happened to say.  Across the
tree that produced 28 runs with the triplet off, almost none of them off on
purpose: the setting was inherited when the case was created for something
else entirely (a solver, an IC, a domain).  A helium-bearing wind has the
metastable level whether or not the input file remembers to ask for it, and
the level is not a diagnostic add-on — it carries its own photoionization,
Penning ionization with H and H2, recombination, and collisional-excitation
cooling, so leaving it out changes the electron budget, the energy budget,
and the 10830 A observable.  Off is the special case, not on.

- `parameters.f90`: `thereis_HeITR = .true.`.
- `input_read.f90`: the key became **optional** (it was read through `req`,
  so the initializer never actually applied to a valid run).  A file that
  omits the line gets the triplet; `Include He23S? False` is the deliberate
  opt-out and is what the HeITR-off regression branch uses.  `thereis_He`
  false still forces it off, as before.
- The Tk interface (`src/utils/EXHALE_interface_functions.py`) follows: the
  checkbutton variable starts at 1, the reset button restores it to 1
  instead of clearing it, and `load_input` treats the key as optional with
  the same semantics as the Fortran reader — on unless the file says
  `False`.  Since a file may now omit the line, the sequential reader
  checks the label before consuming the line, so an older file without it
  no longer shifts `Load IC`, `Do only PP`, and `Force start` by one.

Verified: a one-step run with the key absent prints "Including helium
triplet chemistry" in the setup report, and `Include He23S? False` does not;
the Tk reader parses the tutorial input with and without the line to the
same `Load IC` / `Do only PP` values.  `make check` is unaffected — every
`backup/regression/*/input.inp` sets the key explicitly, including the two
`wasp_he23off*` cases that exist to exercise the off branch.

The `examples/01`--`12` ladder keeps `Include He23S? False` on purpose:
each folder is its base plus exactly one line, and atomic helium is the
baseline that leaves `06_he23s` a one-line difference against `03_newton`
instead of a folder that differs in nothing.  Outside the ladder the triplet
is on -- including `examples/13_lower_atmosphere/`, which had been running two
of its four planets with the triplet and two without.  That split mattered:
the He 2^3S + H2 Penning channel is active only when a molecular run also
tracks the triplet, so one example was doing different physics on different
planets.  All four now carry `True`; none of them has stored output.

The runs whose stored output predates the flip are listed in
`docs/he23s_default_recompute_list.md`; nothing there has been recomputed.

## 66. Primitive conversion of unfilled ghost cells (2026-08-24)

A runtime-checked JFNK run of the LHS 1140 b He/H = 1000 case (Phase C
audit, `LHS1140b/exhale/audit_summary.md`) ended with gfortran reporting a
signalling `IEEE_INVALID_FLAG`.  Rebuilt with `-ffpe-trap=invalid,zero,
overflow`, the solver died in its first residual evaluation at
`UW_conversions.f90:58`, `W(2,:) = U(2,:)/U(1,:)`: `U_to_W` converts the
whole array, ghost cells included, and `Apply_BC` called it on a state whose
ghosts were still zero -- the comment two lines above the call in
`steady_newton.f90` had already noted that the ghosts are zero at that
point.  The `0/0` produced NaNs in the ghost primitives that `Apply_BC_W`
overwrote on the next line, so no result was ever affected; but every
steady solve tripped an invalid operation, which is what had kept the
trapping variant of the `-fcheck` regression off the steady-solver path.
With that site fixed the trap moved to the two line-search trial
conversions in `steady_newton.f90` (`call U_to_W(utry, Wtry)` in the hybrid
and JFNK searches), where `unpack_U` fills only the interior of `utry` and
only `Wtry(1,1:N)` and `Wtry(3,1:N)` are then read.

- `UW_conversions.f90`: new `U_to_W_interior`, converting cells `1..N`
  through `U_to_W_comp` and zeroing the ghost columns -- one definition
  for every site that converts a state before its ghosts are set.
- `Apply_BC.f90`: `U_to_W` -> `U_to_W_interior`; the ghosts are outputs of
  `Apply_BC_W`, never inputs.
- `steady_newton.f90`: both line-search trial conversions ->
  `U_to_W_interior`.

Verified: under `-O1 -g -fcheck=bounds,do,mem -ffpe-trap=invalid,zero,
overflow` the He/H = 1000 JFNK now completes, `info=0 ||R|| = 7.929E-04`
(identical to the untrapped run) with no signalling exception, and so does
the H-rich `solar_gj699` case (`info=0 ||R|| = 6.433E-04`).  `make check`:
5/5 byte-identical (all five cases), as the discarded ghost values require.

## 67. Element diffusion on the direct-steady route, the element drift metric, and the restart of a diffused state (2026-08-25)

Three defects on the `He_diffusion` path, all of them independent of the
transport formulation and all of them blocking the Phase-D milestones of
`docs/binary_diffusion_design.md` (sections 7.3, 7.4).

**The steady solve and the diffusion relaxation are now one routine.**
The steady residual (`steady_newton.f90`) carries no diffusion: the He/H
field is moved by the operator-split `he_diffusion_step`, so a single
steady solve freezes the composition at the state it was handed.  The
marching hand-off answered that with an outer Picard loop (JFNK -> 500
relaxation steps -> `ioniz_eq` -> drift test, at most five passes), but
the direct-steady route -- `EXHALE_PTC=1`, the route the LHS 1140 b cases
take -- ran a bare `solve_steady_jfnk` / `solve_steady_ptc`, so a diffused
wind could not be brought to a steady state there at all.  That loop is
now `steady_wind_with_element_diffusion` (contained in `EXHALE_main`),
called from both routes; the solver, its iteration budget and its initial
pseudo-time are its arguments, so the direct-steady route reaches it with
JFNK or with pseudo-transient continuation.  `dt_loc`, which the
relaxation step needs, is built by the `eval_dt` that route already
called.  With `He_diffusion` off the loop runs exactly once and the body
is the single solve plus the state refresh both sites performed before.

**The drift metric counts elements.**  It read
`(f_sp(:,3)+f_sp(:,4)+f_sp(:,5))/(f_sp(:,1)+f_sp(:,2))`, which is HeI +
HeII + HeIII over HI + HII: it omitted the He 2^3S triplet and every
molecular species, and so did not measure the same quantity that
`species_diffusion.f90` transports (its own `sumHe` includes the triplet).
Replaced by `element_ratio_HeH` in `composition.f90`, which counts nuclei
with the `bsp_nH` / `bsp_nHe` weights of `species_table.f90` over all ten
base species -- H2 and H2+ two H nuclei, H3+ three, HeH+ one of each.  It
is the same element count `load_IC` applies to a restart file, and it is
correct in the molecular region without a second definition.

**A diffused state can be restarted.**  `load_IC.f90` rescaled every cell
to the global `HeH` whenever the file's composition differed from it by
more than 1e-6, on the ground that the input file is the authority on
composition.  With `He_diffusion` on that erases the answer: the
cell-by-cell split is what the diffusion operator produced.  The rescale
is now confined to the base cells (`1-Ng` to `1`), where the operator
itself holds a Dirichlet reservoir composition; every cell above keeps the
helium fraction it was written with.  With the flag off the whole column
is rescaled as before, bit for bit.  The HeH+ refusal is unchanged --
molecular chemistry and `He_diffusion` remain an excluded pair until the
molecular closure is validated (memo section 5, milestone M4).

No diffusion operator changed; `species_diffusion.f90` is untouched.

Verified.  `make check`: `==> REGRESSION PASS (all cases byte-identical)`,
5/5, which is the evidence that the `He_diffusion`-off path is unmoved.

The restart was measured against the stored `examples/14_diffusion` output
under `EXHALE_DUMP_IC=1`, which writes the state as loaded, before the
first ionization-equilibrium solve.  The loaded He/H profile reproduces the
file's to 3.3e-16 relative, worst cell, over all 504 cells; the same file
through the pre-change binary comes back flattened to the input ratio
(r = 1.05 R_p: file 1.616e-02, loaded 8.333e-02, a factor 5.2).  The
direct-steady route with `EXHALE_PTC=1 EXHALE_PTC_JFNK=1` now runs the
outer loop -- five passes, `info=0` on every one, drift falling 9.40e-03,
7.74e-03, 5.79e-03, 4.64e-03, 4.03e-03.  It reaches the pass limit rather
than the 1e-3 drift test, which is what a five-pass Picard budget gives
from a state that was not converged to begin with (the stored output is the
unconverged JFNK `info=2` one); the same test from a converged state is
part of the baseline capture below.

The pre-change binary is kept as
`backup/phase_d_baseline/EXHALE_prephaseD.x`, with the reference runs it
produces, as the memo's T2 comparison.

## 68. Binary H/He element diffusion replaces the trace-helium kernel (2026-08-25)

The composition operator of `He_diffusion` is rewritten from the
trace-helium-in-hydrogen kernel to the two-component (binary) formulation
of `docs/binary_diffusion_design.md`, sections 2-4.  New module
`src/modules/functions/binary_element_diffusion.f90`, routine
`element_diffusion_step`; the old module is kept for the duration of
milestone M2 with its routine renamed `he_diffusion_step_trace` and its
body untouched, reachable only through `EXHALE_DIFFUSION_KERNEL=trace`,
so the two can be compared (test T2a).  Both call sites in
`EXHALE_main.f90` go through one dispatcher, `element_composition_step`.

**What changed in the physics.**

*The transported variable.*  The trace kernel evolved the ratio
`f = n_He/n_H`, which diverges as hydrogen becomes the minor element --
exactly the limit a helium-rich atmosphere needs -- and it wrote the
diffusive flux as `n_H D d(n_He/n_H)/dr`, which divides by a vanishing
background.  The new operator transports the helium MASS FRACTION
`X = rho_He/rho` of two components, hydrogen (with the trace metals
slaved to it) and helium, and is bounded at both ends of the composition
axis.  The gradient coefficient of the memo, `A = (m_1 m_He n/mbar) D_12
(dx/dX)`, collapses to `A = rho D_12`, so nothing in the operator divides
by an element density.  The settling flux carries the prefactor
`rho X (1-X)`, which switches settling off as either element is
exhausted, while the gradient flux does not vanish there -- a helium-free
cell next to a helium-bearing one still receives helium.

*The cap is gone.*  The trace kernel clamped `f <= HeH` over the whole
domain, and in the stored `examples/14_diffusion` output the profile sat
exactly on that clamp between 1.2 and 2 R_p.  The clamp existed because
the trace flux has no reason to vanish when hydrogen runs out; the binary
flux does.  It is removed, together with the base pile-up limiter.  A
pile-up, if one appears, is physics for the base boundary condition to
answer, not something to clamp.  What remains is a round-off assertion
(`0 <= X <= 1`), and the bounds it asserts are a property of the
discretization: the backward-Euler matrix is a strictly diagonally
dominant M-matrix, with the Peclet hybrid on the settling term chosen so
that it never flips the sign of an off-diagonal.

*The ambipolar field is computed, not assumed.*  The previous effective
settling mass `3 - 0.5 (Zbar_He - Zbar_H)` is the hydrogen-plasma value
`eE = m_H g/2` written as a constant; in a fully ionized helium plasma the
massless-electron balance gives `eE = 4/3 m_H g` instead, and the relative
settling mass 5/3, not 2.5.  The field is now evaluated from the solved
electron pressure of the cell, `eE = -(1/n_e) d(n_e k T)/dr`, in its
logarithmic form `-k T dln(n_e T)/dr` (exact for a barometric
stratification, the same second-order central difference elsewhere).  The
three analytic limits are reproduced to round-off (test T9): neutral 3,
H+ plasma 2.5 -- so the validated helium-trace behaviour is unchanged --
and He++ plasma 5/3.

*Advection.*  The operator receives only the new `rho`, `v`, `T` and `dt`,
not the hydro's Runge-Kutta face mass fluxes, so a second conservative
advection of `rho X` could not be made consistent with the hydro's own
continuity step.  `X` is therefore transported in the advective form,
`rho dX/dt + rho v dX/dr = -div(r^2 J)/r^2`, with the velocity divergence
cancelled analytically in the row assembly rather than formed and
subtracted (that cancellation is catastrophic whenever the step is long
compared with the advective crossing time).  A uniform `X` is then
preserved exactly whatever the hydro does, and at a steady state, where
`rho v r^2` is constant, the advective and conservative forms are the same
equation.  In a transient the helium mass is conserved only to the order
of the hydro's own truncation; this is the stated limitation of the phase.

*Write-back.*  The species vector is projected onto the new element totals
by a nonnegative projection: helium-only species scaled by
`r_He`, hydrogen-only species (H2, H2+, H3+ included) and the slaved
metals by `r_H`, `HeH+` by `min(r_H, r_He)`, and the shortfall of the
element with the larger factor deposited into its neutral ground species.
Every operation is a multiplication by a nonnegative factor or an
addition, so no negative intermediate can arise, and both element totals
are met exactly.  The old write-back used one factor per element and
inferred the metal mass as `1 - m_H f_H - m_He f_He`, which counts
molecular mass as metal.

*Trace metals* keep the previous element-by-element treatment behind
`He_metal_diffusion` (a trace element cannot feed back on the background
it moves through); only the background definition is now the nucleus count
over all species rather than `HI + HII`.

**User-visible change.**  `He_Kzz` now defaults to `0.0` instead of
`1.0e9 cm^2/s`.  With the old default every diffusion run silently carried
a constant eddy term; an eddy coefficient is a property of the atmosphere
being modelled, so it is stated rather than inherited.  An input that
relied on the default must now say `He_Kzz: 1.0e9`; `examples/14_diffusion`
already does.  The four stored run directories that carried
`He_diffusion: True` without the key -- `benchmarks/hd209`,
`docs/version_compare/v2_diff`, `docs/version_compare/wasp_v2_diff`,
`scratch_fs/hd189` -- now state `He_Kzz: 1.0e9` explicitly, so each input
records the eddy coefficient its stored output was produced with; no output
was regenerated.  `docs/input_schema.md` K19 and `README_HOWTO.md` are
updated.  `examples/exhale_io.py` gains a derived `Run.heh_profile`, the
He/H element ratio with the same nucleus weights the operator uses, so the
analysis reads one definition.

**Verification.**  A standalone driver, `src/tests/diffusion_tests.f90`,
builds with `make diffusion_tests` and runs the acceptance tests of the
memo's section 6 on synthetic columns and on a stored wind; all 17 checks
pass.

| id | measured |
|---|---|
| T1a | closed isothermal column relaxes to the barometric separation `dx/dr = -x(1-x)(m_He-m_1)g/kT`: max error of `logit x` 4.56e-06 at N = 400; helium mass drift 2.0e-16 at a step of about the cell diffusion time |
| T1b | Dirichlet base, integrated boundary flux accounts for the mass change to 2.2e-11 relative |
| T3 | He/H = 1000 with an outflow: finite, `0 <= X <= 1`, hydrogen rises (X falls below the reservoir value), no cap |
| T4 | closed column, arbitrary initial `X(r)`, 1e4 steps: helium mass drift 2.4e-14 (`eos_metals` on) and 3.7e-14 (off), metals present in both |
| T5 | two-component mass closure `m_1 n_H + m_He n_He = rho` after the write-back: 4.0e-15 relative, worst cell over 200 steps |
| T6 | uniform `X`, `g = 0`, compressing and expanding `v(r)`: 4.4e-14 relative departure after 500 steps |
| T9 | relative settling mass 3, 2.5, 5/3 to 8.3e-13 in the three analytic limits; finite and continuous across a resolved ionization front |
| T10 | second order in `dr` (observed 1.998); the steady state is step-size independent to 1.2e-07 over two decades of `dt` |
| T2a | the operator increment against the trace kernel on a stored LHS 1140 b wind falls off first order in the helium fraction: 2.28e-03, 3.00e-04, 3.04e-05, 3.05e-06 at He/H = 0.0833, 1e-2, 1e-3, 1e-4 -- a factor 10 per decade, and 3.0e-06 at 1e-4 against the memo's 1e-3 gate.  Over 500 relaxation steps the two profiles agree to 4.2e-05, 7.4e-06, 7.8e-07, 7.8e-08 on the same sequence |

T2a is measured with the wind's `rho` and `T` and `v = 0`.  The stored wind
is not steady below r ~ 1.2 -- its mass flux `rho v r^2` varies by orders of
magnitude there and the base cell has `v < 0` -- so with that velocity the
conservative advection of the trace kernel and the advective form differ by
`X div(r^2 rho v)/r^2`, a difference of the STATE rather than of the trace
limit, which does not fall off with He/H (the driver reports that column
too).  The comparison of the transport closure, which is what changed, is
therefore made at `v = 0`.

The measured helium-mass drift of the closed column grows with `dt`
(2.0e-16, 6.2e-12, 3.1e-10, 1.1e-07 at `dt_code` = 1e1, 1e3, 1e5, 1e7 over
1000 steps).  This is the round-off of the tridiagonal solve amplified by
the conditioning of the implicit step, not a leak: the discretization
telescopes exactly, and the drift falls to round-off as soon as the step is
comparable with the cell diffusion time.  The driver prints the table.

T0: `make check` is 5/5 byte-identical, which is the evidence that the
`He_diffusion`-off path is unmoved by any of this.  T7 (the molecular
closure) is milestone M4 and T8 (the moving-wind face flux on a converged
LHS 1140 b run) waits on the baseline capture; neither is claimed here.

---

## 69. Milestone M3: the binary operator on a wind, the LHS 1140 b helium-rich runs, and the removal of the trace kernel (2026-08-25)

Milestone M3 of `docs/binary_diffusion_design.md`: the binary element
diffusion operator of section 68 is run on converged winds -- the HD 209458 b
`examples/14_diffusion` configuration against the captured baseline of the
kernel it replaces (T2b) and the three LHS 1140 b helium-rich cases -- the
elemental face flux is measured on a moving wind (T8), and the trace kernel
is then deleted.

### T2b: the binary operator against the trace kernel on HD 209458 b

Both sides are the same configuration (HD 209458 b, He/H = 0.0833,
`Include He23S? True`, `He_metal_diffusion: True`, PLM+WENO3), Newton
finished, with the same five-pass diffusion outer loop.  The baseline
(`backup/phase_d_baseline/kzz0`, `kzz1e9`) was marched from a cold start with
the pre-Phase-D binary; the binary-operator runs (`new_kzz0`, `new_kzz1e9`)
restart from those converged states on the direct-steady route
(`EXHALE_PTC=1 EXHALE_PTC_JFNK=1`), which is the route section 67 gave a
diffusion loop.  This is a report, not a gate (the memo's section 6): the
baseline profile sits on the `f_He <= HeH` cap over part of the domain, so a
difference there is the cap, not a defect of either side.

| quantity | trace, `Kzz=0` | binary, `Kzz=0` | trace, `Kzz=1e9` | binary, `Kzz=1e9` |
|---|---|---|---|---|
| JFNK `info` / final `\|\|R\|\|` | 0 / 7.737e-4 | 0 / 8.043e-4 | 0 / 8.017e-4 | 0 / 8.308e-4 |
| outer passes / final drift | 5 / 1.22e-2 | 5 / 1.29e-1 | 5 / 3.88e-3 | 5 / 3.19e-2 |
| `log10 Mdot` [g/s] | 9.7882 | 9.7883 | 9.7909 | 9.7910 |
| He 10830 red depth [%] | 9.7970 | 9.8134 | 9.7958 | 9.8118 |
| (He/H)/HeH at 1.05 R_p | 0.62354 | 0.62354 | 0.62806 | 0.62808 |
| at 1.10 R_p | 0.69771 | 0.69744 | 0.70193 | 0.70170 |
| at 1.20 R_p | 0.95514 | 0.95515 | 0.95536 | 0.95536 |
| at 1.50 R_p | 0.97542 | 0.97437 | 0.97555 | 0.97447 |
| at 2.00 R_p | 1.00000 | 1.00111 | 1.00000 | 1.00112 |
| at 3.00 R_p | 0.98628 | 0.98712 | 0.98643 | 0.98725 |
| at 4.00 R_p | 0.85734 | 0.85834 | 0.85806 | 0.85922 |
| cells within 1e-6 of `HeH` | 45 (r = 1.00 .. 2.44) | 3 (the base) | 45 (r = 1.00 .. 2.44) | 3 (the base) |

Three things are worth reading off it.

*The flat stretch was the cap.*  Forty-five cells of the trace profile sat
exactly on `HeH` between 1.0 and 2.44 R_p.  With the binary operator only
the three Dirichlet base cells are on that value, and the profile crosses
*above* it -- 1.0011 at 2 R_p -- which the clamp forbade.  Helium is mildly
enriched there relative to the reservoir, and the enrichment is small.

*The base depletion is not a cap artifact.*  At 1.05 R_p both kernels give
0.6235; the two profiles differ by at most 4e-3 anywhere above 1.05 R_p
after 2500 composition steps of the binary operator.  The depletion of the
first cells (down to 1e-3 of the reservoir ratio in the cell at 1.0004 R_p,
in *both* kernels) is not the limiter and not the lag of the trace form.

**Corrected in section 70 (2026-08-25).**  That last depletion was read here
as settling.  It is not: the barometric He-H separation scale of that base is
24 cell widths, so a 550-fold drop across one cell cannot be barometric, and
the settling flux measured at those faces is 2-20% of the gradient flux.  It
was the face-averaged upwinding of the advective term, which decouples a cell
whose two face velocities straddle zero.  The two kernels agreed because they
shared it.  Section 70 also supersedes the drift and profile numbers of the
tables in this section.

*Neither side is composition-converged.*  Both leave the outer loop at its
five-pass ceiling, and the binary operator's drift is an order of magnitude
larger than the trace kernel's.  That drift is the maximum *relative* change
of He/H over a pass and it is dominated by those near-base cells, where He/H
is ~1e-3 of the reservoir value: over the same five passes the profile above
1.05 R_p moved by less than 0.4%.  The observable moved as little: `Mdot` by
1e-4 dex and the He 10830 red depth by 0.017 percentage points (9.797% ->
9.813%).  The five-pass ceiling is unchanged by this milestone; making the
outer loop reach `drift < 1e-3` needs either many more passes or a drift
metric that is not dominated by a helium-free base cell, and is left open.
(Answered in section 70: neither -- the relaxation was being stepped at the
hydro CFL step, four to six orders below the composition time scale, so it
moved the base by ~1e-6 of the way per step.  With the step taken from the
composition time scale the `Kzz = 1e9` case converges to 1.2e-4 absolute in
ten passes of ~30 steps.)

`load_IC` kept the diffused profile across the restart, which is what
section 67 added: the run reports `top cell He/H = 7.0843E-02` against the
input `8.3333E-02` and pins only the base cells to the reservoir value.

### T8: the elemental face flux on a moving wind

With `EXHALE_DIFFUSION_CHECK=1` the operator writes
`diffusion_faceflux.txt` in the run directory at every step (replaced each
time, so the file holds the state the run ended on): the diffusive face flux
`J`, the advective face flux, their sum as
`F_He = 4 pi r^2 (rho X v + J)`, the total mass flux `4 pi r^2 rho v` and
the relative settling mass, all taken from the same face coefficients the
implicit solve used -- the file reports the discrete fluxes, not a
re-derivation of them.

On `new_kzz0` (HD 209458 b, `Kzz = 0`), with the spread measured exactly as
the code's own `du` is, `(max - min)/|mean|` over a radial window:

| window | spread of `F_He` | spread of `rho v r^2` |
|---|---|---|
| r >= 2.0 R_p (the code's escape window `[j_min:N]`) | 5.13e-3 | 7.19e-3 |
| r >= 1.5 R_p | 1.68e-2 | 7.71e-3 |
| r >= 1.2 R_p | 5.14e-2 | 1.25e-2 |
| r >= 1.1 R_p | 3.84e-1 | 7.77e-2 |

**T8 passes**: in the window the code itself uses to declare the wind
steady, the elemental helium flux is as flat as the mass flux (flatter, in
fact); over the wider windows the two spreads stay within a factor 5, i.e.
the same order of magnitude, and both are dominated by the same unsteady
base region below ~1.1 R_p.  The diffusive part is 12% of the advective part
at most above 1.2 R_p, so this is a real test of the sum and not of
advection alone.

### M3-2: LHS 1140 b with diffusion

Three cases in `LHS1140b/exhale/`, each the seed case's `input.inp` plus
`He_diffusion: True` (`He_Kzz` at its new default 0, `He_ambipolar` at its
default on), restarted from the seed's converged output and finished on the
direct-steady route by `finish_case.sh` (JFNK -> post-processing pass ->
`EXHALE_transit.py`).  All three converged, `info = 0` on every pass.

| | `heh0p55_diff` | `heh1000_diff` | `heh0p06_gj699_diff` |
|---|---|---|---|
| seed | `heh0p55` (GJ 1132 SED) | `heh1000` | `heh0p06_gj699` |
| He/H reservoir | 0.55 | 1000 | 0.06 |
| JFNK `info` / `\|\|R\|\|` | 0 / 4.439e-4 | 0 / 7.929e-4 | 0 / 8.482e-4 |
| outer passes / final drift | 5 / 4.73e-2 | 5 / 5.17e-2 | 5 / 2.20e-2 |
| `log10 Mdot`, seed -> diffused | 7.7642 -> 7.7645 | 7.7442 -> 7.7437 | 8.5885 -> 8.5883 |
| (He/H)/HeH at 1.05 R_p | 0.9997 | 0.9973 | 0.9996 |
| at 1.2 R_p | 0.9996 | 1.0020 | 0.9996 |
| at 2 R_p | 1.0008 | 1.0030 | 1.0010 |
| at 5 R_p | 1.0048 | 1.0074 | 1.0041 |
| at 10 R_p | 1.0085 | 1.0101 | 1.0050 |
| at 20 R_p | 1.0084 | 1.0102 | 1.0034 |
| profile extremes | 0.683 at 29.0 R_p, 1.009 at 14.5 R_p | 0.626 at 29.0 R_p, 1.195 at 1.03 R_p | 0.827 at 29.0 R_p, 1.004 at 6.7 R_p |
| H ionization front (x_HII = 0.5) | 9.43 R_p | 1.03 R_p | 5.66 R_p |
| He 10830 red depth [%], seed -> diffused | 4.2847 -> 4.2875 | 73.389 -> 73.321 | 3.8335 -> 3.8453 |

**The separation is small on these winds.**  Between the base and 20 R_p the
element ratio departs from the reservoir value by ~1% at most; only in the
last few cells below the 30 R_p outer boundary does it fall to 0.63-0.83.
The mass-loss rate and the He 10830 depth move by less than 0.1%.  This is
the expected ordering rather than a null result: these are fast, low-gravity
winds, and the advection time across the domain is short against the
diffusion time, so the wind carries the mixture out before it can separate.
Whether the residual 5e-2 drift would grow the separation if the outer loop
were allowed to run to `drift < 1e-3` is not answered here, and the
five-pass ceiling was not changed for this milestone.

**The ambipolar settling mass on the solved states.**  From the `dmeff`
column of `diffusion_faceflux.txt`, above 1.05 R_p: 2.56-3.00
(`heh0p55_diff`), 2.76-3.05 (`heh0p06_gj699_diff`) and 1.65-5.47
(`heh1000_diff`).  The two hydrogen-rich cases stay near the neutral value 3
because helium is largely neutral there and the H+ term enters through
`Zbar_He - Zbar_1`.  The helium-rich case spans the widest range and touches
1.67 at 1.11 R_p -- the He++ value 5/3 to three digits -- but **it does not
get there as the He++ limit**: `x(He III)` reaches only 0.28 at the top of
the domain in that run, so the 5/3 limit of the memo's (3a) is not sampled
by these winds and remains guarded by the T9 unit test alone.  Below
1.05 R_p `dmeff` is large and of either sign (down to -43 in the base cell
at 1.03 R_p), where the electron pressure gradient is steep and the
diffusive flux itself is negligible.

`LHS1140b/exhale/finish_case.sh` gains one optional variable, `EXHALE_BIN`,
which pins the executable a case is run with (these three were run with a
snapshot binary while the tree was being rebuilt).  With it unset the script
behaves exactly as before.

### The trace kernel is removed

With T2a (section 68) and T2b recorded, `species_diffusion.f90` is deleted,
with it the `EXHALE_DIFFUSION_KERNEL=trace` selector and the
`element_composition_step` dispatcher in `EXHALE_main.f90`: both call sites
now call `element_diffusion_step` directly.  In the test driver T2a can no
longer be measured -- its comparison partner is gone -- so it prints the
M2 measurement as a recorded result and counts no verdict; the driver is
15 checks, all passing.  `make check` stays 5/5 byte-identical (no
regression case turns diffusion on, so this is the statement that the
`He_diffusion`-off path is untouched).

Not done in this milestone: T7 and the molecular closure (milestone M4, the
`error stop` on molecular chemistry + `He_diffusion` stands), and the
regression case with `He_diffusion: True` that section 7.5 of the memo asks
for at the end of the series.

---

## 70. The base helium hole, the metal drain beside it, and a composition relaxation that converges (2026-08-25)

Three defects on the element-diffusion path, all found on the milestone-M3
runs of section 69, all in the same place: the first cells above the base.

### What was measured

On `backup/phase_d_baseline/new_kzz0` (HD 209458 b, `Kzz = 0`, JFNK
`info = 0`) the helium ratio of the first free cell above the base sat at
`(He/H)/HeH = 1.8e-3` between two cells at 1.0 -- a 550-fold drop across one
cell width, where the barometric He-H separation scale of that base is
`H = kT/(dm g) = 0.0045 R_p`, **24 cell widths**. A 550-fold drop across one
cell is not barometric physics. Every metal was emptied from the same cell by
~10^3 (`C I` 2.3e7 cm^-3 against 2.6e10 in its neighbours), and with them the
metal-line cooling of that cell. The dip sat exactly where the breathing base
alternates the sign of the velocity: `v = (-796, +63, +28, -8, -9, +5) cm/s`
over cells 1 to 6.

Rebuilding the operator's face coefficients from the converged state gives
the row terms of the first free cell (cell 2), per unit `X`, in
`g cm^-3 s^-1`:

| cell | `(He/H)/HeH` | diffusive `aa` | diffusive `cc` | advective `aa` | advective `cc` |
|---|---|---|---|---|---|
| 2 | 2.4e-3 | -2.97e-17 | -2.92e-17 | **0** | **0** |
| 3 | 4.5e-1 | -2.78e-17 | -3.22e-17 | -4.53e-15 | 0 |
| 4 | 7.2e-1 | -3.09e-17 | -3.36e-17 | -9.94e-16 | -8.46e-16 |

and the face fluxes it balances, in `g cm^-2 s^-1`:

| face | gradient `J` | settling `J` | total `J` | `rho_f v_f` |
|---|---|---|---|---|
| 1/2 | 1.451e-11 | -3.193e-13 | 1.420e-11 | -7.65e-08 |
| 3/2 | -6.098e-12 | -1.459e-13 | -6.244e-12 | 8.78e-09 |
| 5/2 | -4.114e-12 | -3.575e-13 | -4.472e-12 | 1.93e-09 |
| 7/2 | -2.011e-12 | -4.723e-13 | -2.484e-12 | -1.64e-09 |

The settling flux is 2-20% of the gradient flux and the Peclet number of the
settling hybrid is 0.02 everywhere at the base, so the hybrid is on its
central branch and the settling is fully resolved: **the depletion is not
settling**. Both advective off-diagonals of that row are *identically zero*,
so the cell was held by its diffusive coupling alone -- and that coupling
gives a relaxation time `rho/(|aa|+|cc|) = 3.3e6 s` against the ~2 s hydro
step the relaxation was being driven at. Whatever composition the cell
happened to hold was frozen there.

### (a) The advection was upwinded on face-averaged velocities

`face_coefficients` built the advective coefficients from
`v_f = 0.5 (v_j + v_{j+1})` and chose the donor per face. A cell whose two
face velocities straddle zero (`v_f(j-1) < 0 < v_f(j)`, which is what the
numbers above give: `-366` and `+45` cm/s) then donates at both faces and
receives at neither, and its advective term vanishes identically -- even
though its own `v_j = +63 cm/s` is large and its donor, the Dirichlet base
itself, sits one cell away.

The advective term is now the one-sided upwind difference on the **cell**
velocity,

```
rho_j v_j (X_j - X_{j-1})/(r_j - r_{j-1})   for v_j >= 0
rho_j v_j (X_{j+1} - X_j)/(r_{j+1} - r_j)   for v_j <  0
```

which contributes `+cadv` to the diagonal and `-cadv` to the donor and to
nothing else: a uniform `X` stays exact (T6), the row sum of the advective
part is still zero, the M-matrix argument is unchanged, and the flux-minus-
divergence construction the old form needed is gone. On the cell above, the
coupling it restores is `6.43e-15`, **216 times** the diffusive coupling that
was holding the cell.

New unit test **T11**: a column seeded with the measured 500-fold hole in the
first free cell and a velocity that alternates in sign cell by cell -- so
every face average vanishes while the cell velocities do not -- must refill
the hole. It gives 6.7e-5 relative departure from the reservoir with the new
form and 2.8e-2 (28 times the 1e-3 tolerance) with the old one. The driver is
now 16 checks.

### (b) The trace-metal solver drained the same cell

`solve_trace_element_in_hydrogen` transported the metal DENSITY
conservatively with the same face-averaged velocities. In a conservative form
the configuration above is not a decoupling but a drain: the cell has an
outflow term at both faces and an inflow term at neither. That divergence has
no counterpart in the hydro's own `rho`, whose face fluxes come from the
Riemann solver.

The solver now transports the **mixing ratio** `f_X = n_X/n_H` in the same
advective form as helium, so a uniform `f_X` is preserved under any velocity
field and no spurious divergence can create or destroy the element. The
`f_X <= f_X(base)` cap in the write-back went with it: settling piles an
element up as readily as it depletes one, and the cap is the same limiter the
memo (section 1) rejects for helium. In the runs below no metal exceeds its
reservoir ratio anywhere, so the cap removal changes nothing in them by
itself; what changes the metal profile is the transport.

### (c) The relaxation ran on the CFL step, and on a flow with no mass flux

`steady_wind_with_element_diffusion` called the operator 500 times per pass
with the hydro's local CFL step `dt_loc`. The implicit operator is
unconditionally stable, so nothing tied its step to a sound-crossing time,
while the base composition relaxes on `dr^2/D_12 ~ 10^6 s` against
`dt_loc ~ 2 s`. That is why every case in section 69 left the outer loop at
its five-pass ceiling with the drift falling by only ~0.75 per pass. The
drift metric compounded it: `max |dHe/H| / (He/H)` is a *relative* measure,
and it was dominated by exactly the emptied cells where it means nothing.

The relaxation is now `relax_element_composition` in the diffusion module:

- the step is the composition time scale, `min(dr^2/(D_12 + K_zz),
  dr/max(|v|, w_s))` per cell with `w_s = D_12 |G|`, grown geometrically
  (x1.5, capped at 10^12 times the start) so the late steps are direct steady
  solves, with the coefficients Picard-updated at every step;
- the measure is absolute, `max_j |X_j^new - X_j^old| / X_base`, for the
  inner stop (1e-12 per step) and for the outer pass stop (1e-3, at most 20
  passes);
- the advecting flow is the wind's **steady mass flux**,
  `rho v = mdot_steady/r^2` with `mdot_steady` the median of `r^2 rho v` over
  the escape window `[j_min:N]`. The advective form is the conservative
  equation only where `rho` and `v` satisfy continuity, and below 1.02 R_p
  the converged states do not: the spread of `r^2 rho v` there is 4.8e2
  (`new_kzz0`), 7.4e3 (`new_kzz1e9`) and 7.3e3 (`heh0p55_diff`) times its own
  median, against 1.7e-2 to 4.6e-1 over `1.1-2 R_p`. Relaxed on that field
  the operator converges to the composition of a flow that neither conserves
  mass nor exists. Solving the operator's steady state offline on the M3
  `new_kzz0` wind: with the cell velocities the base takes a five-fold cliff
  at the fourth cell and a plateau at 0.10 of the reservoir ratio; with the
  steady flux the first free cells are flat to 1% and the profile leaves the
  base on the barometric scale.

The marching call keeps `dt_loc` and the cell velocities: there the wind is
genuinely transient and both are the consistent choice.

`element_ratio_HeH` is no longer used by the outer loop itself; it is kept
as the element count of the diagnostic below, which writes one
profile per pass.

### The runs, re-run from the same restart files

`backup/phase_d_baseline/new_kzz0_fix`, `new_kzz1e9_fix` and
`LHS1140b/exhale/heh0p55_diff_fix`: the section-69 configurations, the same
`output/*_IC.txt` restart files, the same direct-steady route
(`EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0`). The M3 directories
are untouched.

| | `Kzz=0` M3 | `Kzz=0` now | `Kzz=1e9` M3 | `Kzz=1e9` now | LHS `heh0p55` M3 | LHS now |
|---|---|---|---|---|---|---|
| outer passes | 5 (ceiling) | 20 (ceiling) | 5 (ceiling) | **10** | 5 (ceiling) | **6** |
| final drift | 1.29e-1 rel | 1.71e-1 abs | 3.19e-2 rel | **1.18e-4 abs** | 4.73e-2 rel | **1.67e-4 abs** |
| relaxation steps per pass | 500 | 24-33 | 500 | 28-31 | 500 | 28-400 |
| `(He/H)/HeH`, first free cell | 1.51e-4 | 0.9971 | 0.8293 | 1.0000 | 1.0001 | 0.9400 |
| next four cells | 0.381, 0.662, 0.810, 0.871 | 0.9939, 0.9960, 0.9994, 1.0025 | 0.7146, 0.6549, 0.6371, 0.6385 | 1.0000 x4 | 0.9990, 0.9998, 1.0001, 1.0000 | 0.8818, 0.8510, 0.8143, 0.7695 |
| at 1.05 R_p | 0.6235 | 0.9085 | 0.6282 | 0.9823 | 0.9997 | 4.9e-4 |
| at 1.10 / 1.20 R_p | 0.6977 / 0.9551 | 0.8342 / 0.7986 | 0.6967 / 0.9553 | 0.9268 / 0.8962 | 0.9995 / 0.9996 | 9e-5 / 2e-5 |
| at 1.50 / 2.00 R_p | 0.9744 / 1.0011 | 0.7912 / 0.7974 | 0.9745 / 1.0011 | 0.8903 / 0.8959 | 1.0001 / 1.0008 | ~0 |
| at 3.00 / 4.00 R_p | 0.9871 / 0.8583 | 0.7935 / 0.7145 | 0.9871 / 0.8603 | 0.8923 / 0.8279 | 1.0021 / 1.0035 | ~0 |
| `log10 Mdot` [g/s] | 9.7883 | 10.396 | 9.7910 | 9.998 | 7.7645 | 7.6065 |
| He 10830 red depth [%] | 9.8134 | 6.4032 | 9.8118 | 8.2115 | 4.2875 | 1.0e-4 |

The base is what the fix was for and the base is fixed. On both HD 209458 b
cases the first free cells now sit within 0.6% of the reservoir ratio and the
profile leaves them monotonically, instead of the isolated 550-fold hole
followed by a five-cell recovery. On LHS 1140 b the first free cells fall
0.940, 0.882, 0.851, 0.814, 0.770 -- ratios of 0.94-0.96 per cell against the
`exp(-dr/H) = 0.93` of that column's barometric scale, i.e. the base now
leaves the reservoir on the barometric scale, which is the criterion this
change was measured against. The metals follow helium: `C/H` in the first
free cell was 0.57 of the reservoir with `Kzz = 1e9` and is 1.0000 now, and
above 1.2 R_p it settles to ~0.56 where the removed cap pinned it at 1.0000.

Three things this exposes rather than fixes, and none is claimed here:

- **The separation is much larger than M3 reported** -- 0.80 instead of 0.98
  of the reservoir ratio at 1.2-3 R_p for `Kzz = 0`, and on LHS 1140 b
  helium is gone above 1.05 R_p instead of flat at the reservoir. That is
  what a relaxation that reaches its steady state gives; the M3 profiles were
  the restart files barely moved. For the LHS 1140 b column it is also the
  expected physics of `Kzz = 0`: at 226 K and `g = 1.8e3 cm/s^2` the
  homopause of a hydrogen-helium mixture with no eddy mixing sits *at* the
  base, and the He 10830 line goes with it (4.29% -> 1e-4%). Those cases
  need a `He_Kzz` argued from the planet, not the default 0.
- **The `Kzz = 0` HD 209458 b case does not co-converge.** JFNK returns
  `info = 0` on every one of the twenty passes, but the pass drift settles
  into a limit cycle instead of falling, so its `log10 Mdot` of 10.396 is not
  a converged number -- it is one branch of a period-2 cycle whose other
  branch is 9.867. The subsection after next shows the cycle is in the wind,
  not in the composition. The `Kzz = 1e9` case converges (1.18e-4 at pass 10)
  and moved by +0.21 dex against M3; the next subsection attributes that.
- The wind itself is better converged than M3's: the spread of `r^2 rho v`
  over 1.1-2 R_p is 1.15e-2 (`Kzz = 0`) against M3's 1.81e-1.

### Where the +0.21 dex comes from: metal transport, tested directly

`backup/phase_d_baseline/new_kzz1e9_nometdiff` is `new_kzz1e9_fix` with one
key changed, `He_metal_diffusion: False`, from the same restart files and on
the same route. It is the control that separates the two halves of the fix:
helium transport and relaxation on one side, metal transport on the other.
Note what the flag does and does not do -- with metal diffusion off the
metals are not reset to the reservoir ratio, they are FROZEN at whatever the
restart file carried, which here is the M3 profile including its base hole.
So the control is "the new helium side, the old metal profile".

| | M3 (`new_kzz1e9`) | metal diffusion off | full fix |
|---|---|---|---|
| outer passes / final drift | 5 (ceiling) / 3.19e-2 rel | 3 / 1.23e-4 abs | 10 / 1.18e-4 abs |
| `log10 Mdot` [g/s] | 9.7910 | **9.8160** | 9.9976 |
| He 10830 red depth [%] | 9.8118 | **9.1081** | 8.2115 |

`C/H` relative to the base cell (`O/H` is the same to four digits in the
first eight rows), first eight rows and then interpolated:

| | row 1-3 | 4 | 5 | 6 | 7 | 8 | 1.05 | 1.10 | 1.20 | 1.50 | 2.00 | 3.00 | 4.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M3 | 1.0000 | 0.5694 | 0.6266 | 0.6213 | 0.6255 | 0.6379 | 0.628 | 0.696 | 0.953 | 0.975 | 1.000 | 0.944 | 0.573 |
| metal diffusion off | 1.0000 | 0.5694 | 0.6266 | 0.6213 | 0.6255 | 0.6379 | 0.628 | 0.696 | 0.953 | 0.975 | 1.000 | 0.944 | 0.573 |
| full fix | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.9999 | 0.9999 | 0.909 | 0.668 | 0.569 | 0.551 | 0.569 | 0.565 | 0.426 |

The metal-diffusion-off column is identical to M3's, which is the check that
the flag froze the restart profile as intended, and the `(He/H)/HeH` base
cells of that run are 1.0000 across, which is the check that the helium side
of the fix is fully active in it.

**The attribution holds.** Switching the metal transport off returns
`log10 Mdot` to 9.8160 and the He 10830 depth to 9.11%, i.e. essentially back
to M3's 9.7910 and 9.81% and away from the full fix's 9.9976 and 8.21%. Of
the +0.207 dex, +0.182 is the metal transport and +0.025 is the helium side.
(The decision rule as posed had this the other way round: if the metals were
*not* responsible, switching them off would have left `Mdot` near 9.998. It
did not.)

One correction to the wording above: the dominant part of that +0.182 dex is
not the base hole. This control cannot split the metal change into its two
pieces -- the base drain and the settling above the base that the removed cap
previously forbade -- but the base hole spans five cells at 1.0004-1.0010 R_p
while the settling changes `C/H` from ~0.95-1.00 to ~0.55 over 1.2-4 R_p,
where nearly all of the metal-line cooling column sits. On the size of the
regions the settling has to be the larger term. It is a reasoned inference,
not a measurement.

### The `Kzz = 0` non-convergence is the wind, not the composition

The outer Picard loop -- wind at fixed composition, composition at fixed wind
-- has no damping of its own, so the composition update is now under-relaxed,

```
X <- X_old + omega (X_relaxed - X_old)
```

with `omega` starting at 0.5 and halved (floor 0.125) on any pass whose drift
failed to fall. The blend is projected back onto the species vector by the
same nonnegative projection `element_diffusion_step` ends on, so both element
totals are still met exactly. The drift that the convergence test and the
halving rule read is the UNDAMPED distance `max|X_relaxed - X_old|/X_base`,
so a small `omega` cannot buy a false convergence by making the applied step
small. `EXHALE_DIFF_OMEGA=<val>` pins the starting value for experiment, and
`EXHALE_DIFFUSION_CHECK=1` now also writes
`output/diffusion_pass_profiles.txt`, one composition and mass-flux profile
per outer pass -- the question a non-converging loop raises is about the
profiles, not about the scalar drift.

**It does not converge, and the reason is not the composition.** Two runs from
the M3 restart files, both `Kzz = 0`:

| | start `omega = 0.5` (`new_kzz0_fix`) | start `omega = 1` (`new_kzz0_undamped`) |
|---|---|---|
| drift by pass | 9.97e-1, 4.99e-1, then the pass-3 JFNK stagnates (`info = 2`, line search collapses to `lam = 9.5e-7` at r = 1.009, energy row) and the loop exits on its solver-failure guard | falls to `omega = 0.125` by pass 11 and then holds a period-2 cycle: 7.62e-2, 6.84e-2, 8.10e-2, 7.31e-2, 8.43e-2, 7.60e-2, 8.95e-2, 7.97e-2 (passes 13-20) |
| final state | not converged, `\|\|R\|\| = 2.4e-3`, `log10 Mdot = 9.838` | 20 passes, every one `info = 0` |

The profiles written each pass say what is cycling. Over passes 3-20 the median
`r^2 rho v` in the escape window alternates between **2.1-2.2e-7** on odd
passes and **5.9-6.8e-7** on even ones -- a factor 2.8 in mass flux -- while
each branch is individually flat to 1e-3 and every JFNK solve returns
`info = 0` with `||R||` between 9.0e-5 and 9.6e-4, all inside the 1e-3
target. The composition, meanwhile, barely moves between the two:

| r [R_p] | odd (pass 19) | even (pass 20) |
|---|---|---|
| 1.010 | 0.9968 | 0.9974 |
| 1.05 | 0.9645 | 0.9674 |
| 1.10 | 0.9065 | 0.9123 |
| 1.20 | 0.8816 | 0.8893 |
| 1.50 | 0.8768 | 0.8849 |
| 2.00 | 0.8815 | 0.8894 |
| 3.00 | 0.8785 | 0.8865 |
| 4.00 | 0.8219 | 0.8332 |
| `log10 Mdot` | **9.867** | **10.376** |

The largest odd-even difference in `(He/H)/HeH` above 1.05 R_p is 0.0117, and
each parity is nearly stationary from one occurrence to the next (0.0028 odd
to odd, 0.0024 even to even). So this is not a composition chasing a wind: it
is **two wind solutions, 2.8x apart in mass flux, that the residual test
cannot tell apart, sitting on compositions that differ by ~1%.** Damping the
composition cannot remove it, and the measurement above is the evidence --
`omega` reached its 0.125 floor and the cycle held at full amplitude.

Two things follow, neither of them settled here. The `Mdot = 10.396` quoted
for this case in the table above is one branch of that cycle (the other is
9.867, close to the 9.998 of the converged `Kzz = 1e9` case and to the 9.816
of the metal-diffusion-off control). And a residual norm that accepts two
states a factor 2.8 apart in `Mdot` is a weak convergence test for this
configuration, which is the standing caveat behind the flux-based
convergence decision of section 19 (`||R||` is reported for reference only)
rather than anything the diffusion introduced.

### Verification

`diffusion_tests.x`: 16/16 (T11 new; the other 15 unchanged). `make check`:
5/5 byte-identical -- no regression case turns `He_diffusion` on, and the
`he_diffusion` guard at the top of the operator means the flag-off path never
enters any of this.

## 71. The friction the ionized wind actually has: stage-resolved binary diffusion (2026-08-25)

Element diffusion used one coefficient over the whole domain, the neutral
Banks & Kockarts hard sphere. That is the wrong friction above the
ionization front and it is wrong in the direction that matters: a helium or
metal *ion* moving through a proton gas is held by Coulomb collisions, whose
momentum-transfer cross section at `T ~ 1e4 K` is `~1e-12 cm^2` against the
`~1e-15 cm^2` of a hard sphere. With a neutral coefficient everywhere the
solver let helium and every metal keep settling through the ionized wind,
and once the `f_X <= f_X(base)` limiter went away (section 70) nothing
stopped them: on the HD 209458 b `Kzz = 1e9` wind iron and calcium reached
*zero* by 1.1 R_p and magnesium 2.5% of its base ratio by 4 R_p.

Koskinen et al. (2013, section 2.1) solve the same Chapman & Cowling
equation with collision terms that "account for neutral-neutral, resonant
and non-resonant ion-neutral, and Coulomb collisions", and report (their
section 3.2.2) that Coulomb collisions "are much more efficient in
preventing diffusive separation than collisions with neutral H."

This is a correction to the physics, not an option: there is no input key
for it.

### What the coefficient is now

Each element is resolved into its ionization stages and the friction is
built pair by pair, every pair being the Chapman-Enskog first approximation
`D = 3 k T/(16 n mu Omega^(1,1))` evaluated with the cross section that pair
actually has (`binary_element_diffusion.f90`):

| pair | routine | coefficient | scaling |
|---|---|---|---|
| neutral-neutral | `hard_sphere_pair_diffusion` | `1.52e18 (1/A_s+1/A_t)^(1/2) T^(1/2)/n` | `T^(1/2)/n` |
| ion-neutral, non-resonant | `ion_neutral_pair_diffusion` | `1/D = 1/D_pol + 1/D_hs`, `D_pol = k T/(2.21 pi e n (alpha_n mu)^(1/2))` | `T/n` then `T^(1/2)/n` |
| ion-ion | `coulomb_pair_diffusion` | `3 (kT)^(5/2)/(4 (2 pi mu)^(1/2) n (Z_s Z_t e^2)^2 ln Lambda)` | `T^(5/2)/n` |

and the element-element coefficient is the stage-fraction harmonic sum
(frictions add; Blanc's law across stages)

```
1/D_eff = sum_{s in element A} sum_{t in element B} y_s y_t / D_st .
```

Resonant charge exchange `H+ + H` is deliberately absent: both partners
carry the same element, and a binary diffusion coefficient is driven by the
friction *between* the two elements, so every ion-neutral pair that appears
is non-resonant. `HeH+` carries both elements and stays out of both carrier
lists, as it already was out of the mean charges.

**The ion-neutral pair carries two channels.** The polarization (Langevin)
coefficient is the long-range limit: it keeps the induced-dipole attraction
and no repulsive core, so its friction falls as `T^-1` where a rigid core
would hold it at `T^-1/2`. Taken alone it would make an ion-neutral pair
*more* mobile than a neutral-neutral one of the same masses above ~1.5e3 K,
which is unphysical -- opening a second channel cannot make a pair more
mobile. The two are momentum-transfer cross sections of the same encounter,
so to first order their `Q^(1)` add, hence their collision integrals add,
hence their frictions add; and since `D = 3kT/(16 n mu Omega^(1,1))` is
linear in `1/Omega`, adding frictions is adding inverse coefficients,
`1/D_in = 1/D_pol + 1/D_hs`. That is the same rule the stage mixture uses.
Both limits come out exactly: `D_in -> D_pol` where the polarization channel
dominates and `-> D_hs` where the core does, crossing at **1504 K** for an
H/He pair. At the 1e4 K of the transition layer `D_pol` is 2.6-4.7 times
`D_hs` and `D_in` sits at 0.82 of `D_hs`, i.e. 0.18-0.28 of the
polarization value alone.

**Sources.** The hard-sphere form is the one already in the code
(Banks & Kockarts 1973). The polarization constant `2.21 pi` is the one
behind the non-resonant ion-neutral collision frequency of Schunk & Nagy
(*Ionospheres*, eq. 4.88) and the Coulomb coefficient is the standard
Chapman-Enskog result (Paquette et al. 1986; Schunk & Nagy eq. 4.142 give
the equivalent collision frequency). Neither textbook is in `references/`,
so **both numerical constants were re-derived from the Chapman-Enskog
collision integral rather than copied**, and the derivation is written into
the comment beside each routine. The check that the chain is consistent is
that the rigid-sphere collision integral put through the same
Chapman-Enskog expression returns the Banks & Kockarts prefactor exactly,
with a collision diameter of 2.7 Angstrom.

Polarizabilities are the recommended static dipole values of
Schwerdtfeger & Nagle (2019), in atomic units, with `alpha(H) = 4.5 a_0^3`
the exact nonrelativistic ground-state value; H and He come out at 0.667 and
0.205 Angstrom^3, the values Schunk & Nagy tabulate for aeronomy.

**Validity caveat, recorded at the coefficient.** The polarization limit
keeps the induced-dipole attraction and no repulsive core, so its friction
falls as `T^-1` where a rigid core would hold it at `T^-1/2`. It is the
smaller of the two coefficients (stronger friction) below ~1.5e3 K and the
larger above -- by 2.6 to 4.7 at 1e4 K, where the true coefficient is the
smaller because the core contribution adds to the collision integral. The
consequence is confined to the partially ionized transition layer; both
asymptotic limits are exact.

### Two defects fixed with it

- **The trace-metal loop assumed the hydrogen-plasma ambipolar field.** It
  built its settling coefficient from `dm_eff = (A_X - 1) - 0.5 (Zbar_X -
  Zbar_1)`, the constant `eE = m_H g/2` of a pure proton plasma -- the same
  assumption decision D3 removed for helium at M2, still standing here. The
  loop now reads the field `settling_coefficient` computes from the electron
  pressure gradient, so there is one ambipolar field in the module and no
  second, cruder copy. With `He_ambipolar` off the field is zero and the
  loop reduces to the pure-gravity form as before.
- The metal loop also used a single hard-sphere `X`-in-H coefficient; it now
  runs the same stage mixture, metal stages against hydrogen carriers.

### Tests

`src/tests/diffusion_tests.f90` gains **T12**, the three limits of `D_eff`
on a prescribed isothermal column at 1e4 K, each against a closed form
written out independently in the test rather than against the module's own
pair routines: all neutral reproduces Banks & Kockarts to `1.2e-16`
relative, fully ionized reproduces the He++/H+ Coulomb coefficient to
round-off, half ionized reproduces the four-term Blanc average (whose two
ion-neutral terms are themselves the combined form) to round-off. T12d adds
the two temperature limits of that combined form: within `2.5e-2` of the
polarization channel at 1 K, within `1.2e-3` of the hard sphere at 1e9 K,
and never above either. `./diffusion_tests.x`: **21 passed, 0 failed**.
`make check`: **5/5 byte-identical** -- the operator is entered only with
`He_diffusion` on and no regression case sets it.

`EXHALE_DIFFUSION_CHECK=1` adds `D_eff`, the neutral coefficient of the same
cell, their ratio and the dominant stage pair to
`diffusion_faceflux.txt`.

### Measured effect

The full tables are `docs/binary_diffusion_design.md` section 9.3. In short,
on two winds restarted from the same converged state and Newton-finished:

- `D_eff/D_neutral` is 1.00 in the neutral base, 0.16 at the front
  (1.2 R_p), `9.6e-3` at 1.5 R_p and `7.6e-5` at 4 R_p on HD 209458 b -- the
  suppression *grows* outward because the Coulomb coefficient carries
  `T^(5/2)` against the hard sphere's `T^(1/2)`.
- On HD 209458 b every metal mixing ratio becomes **flat above the front**,
  which is what a steady wind must do once diffusion loses to advection.
  Iron goes from 0.000 to 0.107 of its base ratio and calcium from 0.000 to
  0.019 -- the difference between having and not having Ca II H&K and Fe II
  in the transmission spectrum. `log10 Mdot` moves 10.0106 -> 10.0561
  (+0.046 dex) and the He 10830 red depth 8.20% -> 7.85%. The control is a
  run of the same configuration from the same converged state with a binary
  built from this commit's parent, so the friction is the only difference.
- The ion-neutral channel decides how much of that survives. Between 1.0 and
  1.2 R_p hydrogen and helium are still neutral (measured `x(H+)` = 0.000,
  0.004, 0.072 at 1.05, 1.10, 1.20 R_p), so the combined form barely touches
  the *helium* profile -- but the low-ionization-potential metals are already
  ionized there (Fe 0.999, Mg 0.943, Ca 0.269, C 0.329, Na 0.053, O 0.003
  ionized at 1.10 R_p), and going from the polarization channel alone to the
  combined form raises their frozen plateaus in exactly that order: Fe x3.8,
  Ca x2.3, Mg x1.23, C x1.07, Na x1.00, O x0.97.
- On LHS 1140 b at `He/H = 0.55` nothing moves: `log10 Mdot` 7.6065 both,
  the temperature identical to 9e-16 relative, He 10830 within 0.2%, and the
  composition relaxation exits on its first pass. That wind separates its helium inside the first cell above the base,
  where the gas is neutral and `D_eff = D_neutral` by construction, and
  above that there is no helium left for the Coulomb suppression to hold.

---

## 72. The advection correction is gated on the slowest species, not on hydrogen (2026-08-25)

`post_process_adv.f90` replaces the local ionization balance of a cell by the
steady advection-ionization ODE, and skips cells where that replacement can
carry no information. One of the three skip conditions is a Damkohler number,

```
Da = (dr/v) * (P_HI + alpha_HII n_e) > 100   ->   keep the equilibrium state
```

built from the **hydrogen** ionization/recombination rate. The correction it
gates, however, carries the whole solved species vector -- H I/H II,
He I/He II/He III and, when `Include He23S?` is on, the He 2^3S metastable --
so one species' clock decided the fate of all of them. The metastable relaxes
orders of magnitude more slowly than hydrogen (`A31 = 1.272e-4 s^-1` alone
sets a 2.2-hour floor, against sub-second photoionization for H in an
irradiated wind), so wherever hydrogen was equilibrated and the metastable was
not, the gate froze `n(2^3S)` at its equilibrium value in cells where it is in
fact advected.

That is physically wrong independently of how often it bites, and it bit
visibly on the LHS 1140 b run `heh10_diff_kzz1e8` (He/H = 10), where the
equilibrium solution has a sharp hydrogen ionization front: 115 of 503 cells
above the front kept the equilibrium ionization, and the `_adv` metastable
profile stepped from 5.15 to 223 cm^-3 across the single cell boundary
r = 6.083 -> 6.172 R_p (1.64 dex). The transit synthesis reads the `_adv`
profiles, so the run reported a red-pair depth of 83.7% and EW = 23.4 %A -- a
number set by the discontinuity. The limitation was recorded in
`LHS1140b/kzz_decision.md` section 6.2 and the run was excluded from that
memo's crossing fit.

### What the gate is now

The advection systems (`adv_implicit_H`, `adv_implicit_HeH`,
`adv_implicit_HeH_TR`) solve the whole H/He vector at once, so there is no
gate to be had one row at a time: the honest statement is that a cell may be pinned to
equilibrium only when **every** population it solves is equilibrated. The
Damkohler number is therefore formed with the **slowest** relaxation rate of
the solved set, each rate being the total rate at which that population is
destroyed and re-formed:

| population | relaxation rate `nu` |
|---|---|
| H I / H II | `P_HI + (a_ion_HI + alpha_HII) n_e` |
| He I / He II | `P_HeI + (a_ion_HeI + alpha_HeII + alpha_HeI23S) n_e` |
| He II / He III | `P_HeII + (a_ion_HeII + alpha_HeIII) n_e` |
| He(2^3S) | `A31 + P_HeITR + (q31a + q31b + a_ion_HeITR) n_e + Q31 n_HI` |

`Da = (dr/v) * min(nu)`, compared with the same threshold of 100. The
He(2^3S) row is term by term the loss side of `fvec(4)` of
`adv_implicit_HeH_TR` -- radiative decay, photoionization, collisional
de-excitation to 2^1S and 2^1P, electron-impact ionization, and Penning
ionization on neutral H -- and every coefficient is read from the arrays
`HeITR_coeffs` and `eval_cool` already fill in the same routine. No rate is
redefined. The He I row now also carries `rcheiTR`, because with the triplet
on `HeITR_coeffs` overwrites `rcheiiB` with the singlet channel alone and the
He I row of the residuals adds the two. The two other skip conditions --
inflow `v <= 0` on either face, and `x_HII,eq < 1e-6` -- are untouched, so the
breathing-base guard behaves exactly as before.

The hydrogen rate now includes collisional ionization `a_ion_HI n_e`, which
the old expression omitted, for consistency with the three rows beside it. It
can only raise `Da_H`, and `Da_H` is no longer the quantity that decides
anything on its own.

### Measured

Diagnostic build dumping `Da` per cell, LHS 1140 b `heh10_diff_kzz1e8`,
post-processing pass 10:

| r [R_p] | `Da_H`, H-only gate | `Da_H`, new binary | `Da` of the slowest species |
|---|---|---|---|
| 6.083 | 6.1e-4 | 6.1e-4 | 6.1e-4 |
| 6.172 | 3.2e+2 | 6.0e-4 | 2.8e-3 (H-only run), 6.0e-4 (new) |
| 6.946 | 7.8e+5 | 4.8e-4 | 2.3e-3 / 4.8e-4 |
| 7.602 | 4.6e+7 | 1.1e+2 | 2.1e-3 / 1.7e-3 |
| 7.956 | 1.7e+8 | 7.2e+2 | 1.9e-3 / 1.6e-3 |

The metastable Damkohler number is `~2e-3` throughout that region: the
He 2^3S population there is advection-dominated by five orders of magnitude,
in exactly the cells the hydrogen gate pinned to equilibrium. The two `Da_H`
columns differ because the old gate is self-reinforcing -- a pinned cell holds
the fully ionized equilibrium state, whose `n_e` and `P_HI` then keep `Da_H`
above the threshold for the cell above it.

Effect on that run (both binaries run from the same converged state with
`Do only PP: True`, in scratch copies; the case directory itself was not
re-run):

| quantity | H-only gate | slowest-species gate |
|---|---|---|
| cells kept at equilibrium | 115 of 503 | 15 of 503 |
| outermost such cell | 30.0 R_p | 1.0034 R_p |
| max `\|dlog10 n(2^3S)\|` between adjacent cells, r > 1.05 | 1.636 (5.15 -> 223 at 6.13) | 0.386 (4.99 -> 12.1 at 7.43) |
| He 10830 red depth [%] | 83.69 | 21.52 |
| red EW [%A] | 23.39 | 5.69 |
| `log10 Mdot` [g/s] | 7.76 | 7.76 (the wind solution is not re-run) |

The remaining 15 pinned cells all sit at `r <= 1.0034 R_p`, the shielded base
where conditions (i) and (iii) are meant to act. The largest remaining
variation in the metastable profile, 0.386 dex over r = 7.38 -> 7.49 R_p, is
not a switch: it is where the advected hydrogen finally ionizes (`x_HI` falls
from 0.80 to 5e-4 over three cells as `Da_H` crosses 1) and the Penning loss
`Q31 n_HI` on the metastable switches off with it. The `_adv` profiles are
continuous across it, and `hybrd1` returns `info = 1` in every cell of the
region, including those with `Da_H ~ 1e5`.

The corrected run's EW of 5.69 %A is still above what the He/H = 2.7-4.0
sequence of `kzz_decision.md` section 6.1 extrapolates (~2.8 %A at He/H = 10
on that sequence's local log-log slope of 0.73), so the He-rich end appears to
be leaving the regime that sequence describes for its own reasons; what the
fix removes is the discontinuity, not the run's distance from the trend.

### Effect on existing results

Three reference cases were re-post-processed with both binaries in scratch
copies -- `LHS1140b/exhale/heh0p55` (diffusion off), `heh2p13_diff_kzz1e9`
(one of the crossing brackets), and the H-rich standard
`backup/phase_d_baseline/new_kzz1e9_d3b` (HD 209458 b, `Kzz = 1e9`, metals
on). In all three, `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` come out
**byte-identical** to the H-only-gate binary, and the number of pinned cells
is unchanged (42, 30 and 153 respectively): in those runs the pinned cells are
held by the inflow or `x_HII,eq` condition, which the change does not touch,
not by the Damkohler condition. He 10830 red depth and EW therefore move by
0.00%, and no stored result needs regenerating.

`make check` passes 5/5 -- as it must, since the goldens compare
`Hydro_ioniz.txt` and `Ion_species.txt`, which the post-process never writes;
running it is the check that the change stays inside the post-process.

### Re-post-processing the stored LHS 1140 b scans (2026-08-26)

The three reference cases above were the wrong sample for the He-rich end.
All eighteen diffusion-off cases of `LHS1140b/exhale/` -- the eleven of the
GJ 1132 composition scan and the seven `_gj699` companions -- were therefore
re-post-processed in place, each from its own untouched JFNK solution of
record (`output/*_IC.txt`) through one `Do only PP: True` pass and both
transit syntheses; the pre-fix products are kept beside each case in
`pre_gate72_adv/`. Five cases changed, all helium-dominated: `heh10`,
`heh100`, `heh1000`, `heh10000` and `heh1000_gj699` lose 89-93% of their
red-pair depth and equivalent width (`heh10` 75.97% and 21.53 %A to 7.42%
and 2.270 %A; `heh1000_gj699` 105.1% and 35.87 %A to 9.53% and 3.479 %A,
the fitted depth above 100% being the pre-fix saturated trough). The count
of `_adv` cells still holding exactly the equilibrium metastable density
falls from 173-437 of 504 to 13-22 in those five. The remaining thirteen
move by at most 0.06% and keep their cell counts exactly, so the earlier
sample's conclusion holds wherever the hydrogen front is not sharp; the
0.06% is a drift in the equilibrium solve between the binary that wrote
those outputs (2026-08-24) and this one, at most `4e-4` relative, not the
gate.

Nothing the memo `docs/lhs1140b_exhale_vs_pwinds.tex` concludes moves: the
equivalent-width crossing stays at He/H = 0.550 (0.541 with turbulence) on
the GJ 1132 proxy and 0.060 on the GJ 699 one, because both are
interpolated inside He/H <= 1, and so does the matched broadening kernel,
22.30 km/s. What moves is the He-rich end of the two composition tables,
the bump table's He/H = 1000 rows, and the figures built from them.
`heh10_diff_kzz1e6`, the He-rich bracket of the `K_zz` crossing table, was
checked in a scratch copy and is unchanged, so `tab:kzzcross` stands as
published; `heh10_diff_kzz1e8` was already excluded from it.

One artifact is not the gate's and survives it: `heh100` still steps 3.8
dex in `n(2^3S)` across a single cell at r = 1.059 R_p, on the hydrogen
ionization front, where the base velocity is still sign-changing. The shell
is thin enough (0.008 R_p against a 6.8 R_p stellar radius) that it can
account for at most 0.04% of that case's 5.13% red depth.

## 73. Milestone M4: element diffusion through the molecular region, and the end of the molecular / `He_diffusion` exclusion (2026-08-26)

`input_read.f90` refused `Molecular chemistry: True` together with
`He_diffusion: True` -- an `error stop` that
`docs/binary_diffusion_design.md` section 5 put there deliberately, to be
lifted by a validated closure and not by deleting the check. The closure is
now built and the check is gone.

### What the molecular region changes, and what it does not

Equation (1) of the memo transports **elements**, and the chemistry moves
hydrogen between H, H+, H2, H2+ and H3+ without changing the hydrogen
element's mass, so the transport equation is the same above and below the
molecular front. What the species set changes is the closure: helium's
collision partners are then a mixture of carriers, and gravity, the ambipolar
field and the friction all act on **particles**, not on nuclei. Four things
follow, three of them already anticipated by section 5 and one not.

1. **Friction: Blanc's law over the carriers.** The stage sum
   `1/D_eff = sum_s sum_t y_s y_t/D_st` of section 2.6 already ran over the
   carrier lists `{HI, HII, H2, H2+, H3+}` and `{HeI, HeII, HeIII, HeTR}`,
   each carrier with its own mass, charge and polarizability (H2 enters at
   mass 2 and with `alpha(H2) = 5.315 a_0^3 = 0.787 A^3`, the Schwerdtfeger &
   Nagle 2019 recommended value, alongside the H and He entries that were
   already there). It reduces to `D(He,H)` in the atomic region and to
   `D(He,H2)` in a fully molecular one.
2. **Forces per carrier.** `G` now uses the mean carrier mass
   `m_c1 = m_1 * mcar` and the mean carrier charge per carrier, where `mcar`
   is the hydrogen nuclei bound into one collision partner. In the atomic
   region `mcar = 1` and this is the old nucleus form exactly; in an H2 base
   the relative settling mass is `4 - 2 = 2` instead of `4 - 1 = 3`. The
   trace-metal loop settles a metal against the same mean carrier mass.
3. **The carrier density in the pair coefficients.** `n` in the
   Chapman-Enskog `D_st = 3kT/(16 n mu Omega^(1,1))` is a density of
   *colliding particles*. It was the nucleus density; it is now the carrier
   density, which is the same number wherever the hydrogen is atomic and
   smaller below a molecular front.
4. **The mole-fraction driver the chemistry carries** -- the term section 5
   did not state, and the one that makes the other three a force rather than
   bookkeeping. The transported variable is the mass fraction `X`; the
   Chapman-Cowling driving force is the gradient of the **mole** fraction
   `x` among the carriers. With `psi = n_1c/n_H` (1 atomic, 1/2 fully H2),

   ```
   logit(x) = logit(X) + ln( m_1/(m_He psi) ),
   d logit(x)/dr = d logit(X)/dr - dln(psi)/dr
   ```

   The gradient coefficient is untouched -- `A = (m_1 m_He n/mbar) D (dx/dX)`
   collapses to `rho D` with the carriers exactly as it does with the nuclei,
   for any `psi`, the two carrier factors cancelling -- so the whole content
   of the term is one addition to the settling coefficient,
   `G -> G - dln(psi)/dr`. The discretization, the Peclet hybrid and the
   M-matrix argument are unchanged, and the term is identically zero wherever
   the hydrogen is atomic.

   It is a real driver: where hydrogen turns molecular going down, each
   nucleus is spread over fewer collision partners, helium's mole fraction
   rises, and helium diffuses down that gradient -- outward across the
   molecular front -- even at a uniform mass fraction. Its size is set by the
   sharpness of the front, and across a front where `psi` runs 1/2 -> 1 it is
   a factor 2 in the He/H nucleus ratio (test T7d below).

The eddy term carries none of this: eddy mixing transports the mixture as a
whole and has no preferred species, so it acts on `dX/dr` alone.

### T7, and the gate it closes

Four checks in `src/tests/diffusion_tests.f90`, each against a closed form
typed out in the test rather than read back from the module:

| id | check | result |
|---|---|---|
| T7a | forced atomic column: `D_eff` vs `D(He,H)` | 1.13e-16 relative |
| T7b | forced fully molecular column: `D_eff` vs `D(He,H2)` | 1.57e-16 relative |
| T7c | prescribed molecular column, three constant `K_zz`: the homopause read off the relaxed profile vs the radius where `D_eff = K_zz` | 3.4e-3 R_p, about half a cell, at all three |
| T7d | same column with `g = 0`, where `-dln(psi)/dr` is the only driver left | the **mole** fraction levels out to 6.2e-5 while the He/H nucleus ratio moves by 0.999 |

T7c reads the homopause off the profile without assuming where it is: at the
zero-flux steady state of a motionless column,
`rho (D + K) dX/dr = -rho X(1-X) D G`, so `-[d logit(X)/dr]/G` is exactly
`D/(D+K)` and the radius where it passes 1/2 is the radius where `D = K`. The
three `K_zz` were chosen to put that radius at 1.5, 2.0 and 2.5 R_p, and the
measured profile put it at 1.4993, 2.0013 and 2.5010.

T7d is the sharpest statement of the closure: with no gravity, no field and no
eddy, diffusion levels the **mole** fraction of helium among the collision
partners, not its mass fraction, so a column that is molecular below and
atomic above ends with an He/H nucleus ratio that differs by the factor `psi`
across the front and a mole fraction that is flat.

`make diffusion_tests && ./diffusion_tests.x`: **28 passed, 0 failed** (21
before this change; the four new ones and T7c three times). Rebuilt with
`-O1 -g -fcheck=bounds,do,mem -ffpe-trap=invalid,zero,overflow` -- the new
carrier arrays are new array bounds -- the same 28 pass with no runtime
error.

### The exclusion, and the restart

`input_read.f90` no longer refuses the pair; the paragraph there now says why
it does not. `docs/input_schema.md`, `README_HOWTO.md` and
`docs/oxygen_chemistry_new_plan.md` are corrected where they stated the
exclusion.

`load_IC.f90` refused to restart *any* molecular state whose composition
differed from the input `He/H`, because `HeH+` carries a nucleus of each
element and a single rescaling factor cannot set both counts. With
`He_diffusion` on the only cells rescaled at all are the base and its inner
ghosts (the column above keeps its diffused split, section 7.4 of the memo),
so those cells are now projected onto their two element totals the same way
`binary_element_diffusion` writes back -- `HeH+` by the smaller of the two
factors, the shortfall deposited into the neutral ground species. Without
`He_diffusion` the whole column would have to be rescaled and the refusal
stands.

Found while reading that block and fixed with it: the helium nucleus count
`load_IC` builds to decide whether a restart file's composition matches the
input left **`HeTR` out** (`HeI + HeII + HeIII + HeH+`), while the rescaling a
line below scales `HeTR` like every other helium species. The metastable
triplet carries a helium nucleus -- `bsp_nHe(HeTR) = 1`, and every other
element count in the code includes it -- so a restart of a `He23S` state read
its own He/H slightly low and rescaled by that error. It is now in the count.
No golden moves: no regression case restarts.

### End to end

`backup/regression/mol_base_handoff` (the hot Uranus, `Molecular chemistry`
+ a `base.inp` handoff) re-run in a scratch copy with `He_diffusion: True`
added, `EXHALE_MAXSTEPS=12000` as in the regression:

- runs to completion with no NaN and no negative density;
- the two-component mass closure `m_1 n_H + m_He n_He = rho` holds to
  **4.4e-16** relative over the whole grid, read back out of
  `Ion_species.txt` and `Hydro_ioniz.txt`;
- `log10 Mdot = 10.58`, the same as the diffusion-off case at the same step
  count;
- `(He/H)/HeH` stays within 1.3e-3 of the reservoir everywhere.

The flat composition is the physics of this planet, not an inert operator.
`He_Kzz = 1e9` comes from the case's own `base.inp` and is three decades above
the `D_eff ~ 8e5 cm^2/s` of the molecular base, so the column is eddy-mixed;
but the same run repeated with `He_Kzz = 0` and `Kzz_base 0.0` comes out flat
as well, within 4e-3 of the reservoir, and the reason is in the diagnostic:
at the base face the diffusive helium flux is `-3.1e-13 g cm^-2 s^-1` against
an advective `-3.5e-7`, six orders of magnitude smaller. This is a hot Uranus
losing `10^10.6 g/s`; the wind sweeps the composition along faster than either
diffusion or settling can separate it. What the run demonstrates is that the
operator carries the molecular composition without NaN, without a negative
density and with the element budget closed to round-off -- the separation
physics itself is exercised by T7 and by the atomic-region cases of sections
69-71.

That the molecular carriers are the ones being used is visible in
`diffusion_faceflux.txt` (`EXHALE_DIFFUSION_CHECK=1`) of the same run, which
names the stage pair carrying most of the friction:

| r [R_p] | dominant pair | relative settling mass | `D_eff/D_neutral` |
|---|---|---|---|
| 1.000 | HeI-H2 | 2.0005 | 0.7747 |
| 1.107 | HeI-H2 | 2.062 | 0.7855 |
| 1.169 | HeI-H2 | 2.520 | 0.8775 |
| 1.196 | HeI-HI | 2.923 | 0.9783 |
| 1.212 | HeI-HI | 3.000 | 0.9999 |
| 1.691 | HeII-HII | 3.067 | 2.52e-3 |

The relative settling mass runs from 2 to 3 across the molecular front,
exactly the carrier weight of point 2 above, and `D_eff/D_neutral` from
`sqrt(0.75/1.25) = 0.7746` -- the He-H2 hard sphere against the He-H one --
to 1, before the Coulomb suppression of section 71 takes over above the
ionization front.

### Regression

`make check`: 5/5 byte-identical at the time of the change. The operator is
entered only with `He_diffusion` on, which no case set then, and every
molecular quantity added here (carrier lists, mean carrier mass, carrier
density, `dln(psi)/dr`) collapses onto its nucleus form identically when the
molecular columns are zero.

### The diffusion regression case (same day)

Byte-identity of the old cases only shows that the operator stays out of the
way when it is off; nothing in `make check` ran it at all
(`docs/binary_diffusion_design.md` section 7.5 had been carrying that as an
open item since the design was written). The sixth default case,
`backup/regression/mol_diffusion`, closes it. It is `mol_base_handoff` --
the Tier-2 hot-Uranus gate with molecular chemistry on and its pinned
`base.inp` copied unchanged -- plus `He_diffusion: True` and
`He_Kzz: 1.0e9`, run as the same 12000-step relaxation snapshot. One case
covers the widest path: molecular carriers below the front, the
stage-resolved pairs through it, the projection back onto `f_sp`, and the
Coulomb friction of section 71 above the ionization front, which is why this
gate was chosen over the atomic HD 209458 b run the design had sketched.

The measured effect is the one the face-flux diagnostic above predicts for
this planet -- present but small (the same hot Uranus, whose wind sweeps the
composition along faster than diffusion separates it). Against `mol_base_handoff` at the same step count,
`Hydro_ioniz.txt` and `Ion_species.txt` both differ, while

```
mol_base_handoff  final: count=12000  du= 2.6292E+00  dtu= 5.7467E-03   log10 Mdot = 10.58
mol_diffusion     final: count=12000  du= 2.6290E+00  dtu= 5.7474E-03   log10 Mdot = 10.58
```

so the golden pins the operator's arithmetic, not a large dynamical
signature. `make check` is 6/6 with the five existing goldens untouched.

## 74. What a `base.inp` key is allowed to do: five categories, and elemental reservoirs that reach `melem_ab` (2026-08-26)

`base.inp` carried six keys and no statement of what any of them meant for the
wind. That is P2 of `docs/oxygen_chemistry_new_plan.md`, and the reason it
matters is not tidiness: `q_H2_base` had been described in places as pinning
the base composition, when all it does is set the base particle count, and the
handoff had no way at all to carry an element -- the quantity EXHALE actually
transports.

### The contract

Every key now belongs to one of five categories, and the category states what
the value may do. The table lives in three places that must agree:
`read_base_inp`'s header comment, `docs/input_schema.md` section 2c, and here.

| Key | Category | Reaches |
|---|---|---|
| comments only | provenance | nothing (the machine-readable keys are A1a/P0, not yet built) |
| `T_base`, `r_base`, `p_base`, `q_H2_base` | EOS boundary | `T0`, `R0`, `p_base_bar`, `q_h2_base` |
| `HeH_base`, `<El>_H_base` | elemental reservoir | `HeH`, `X_<El>` -> `melem_ab` |
| (none today) | initial guess | - |
| `Kzz_base` | boundary constraint | `he_kzz` |
| `q_H2O`, `q_CO`, ... | diagnostic | nothing; metadata, not physics |

`q_H2_base` is an **EOS anchor**: `comp_ntot_bc` removes the H nuclei bound
into H2 at the base, and that is the whole of its effect. Nothing holds H2 at
that value, no H2 profile is seeded from it, and the molecular network moves
away from it in the first cell. The A1c species keys stay comments because no
part of the code reads them.

### Elemental reservoirs

`<El>_H_base` is new: `C_H_base`, `N_H_base`, `O_H_base`, `S_H_base` and the
same form for `Mg Si Ca Na K Fe`, each the El/H **nuclei** ratio at the handoff
level. It overrides `metals.inp` for that element, and it turns the metal
system on by itself when it is the only nonzero abundance. Hydrogen and helium
already had this: `HeH_base`.

The change that makes it work is one of order. `thereis_metals`, `melem_ab` and
`thereis_lowIP_metal` were derived immediately after the core `input.inp`
block, roughly seven hundred lines before `read_base_inp` ran, so nothing the
handoff said about an element could have been heard. They are now derived in
the composition block after `read_base_inp`, beside the `thereis_He` / HeITR /
molecular reconciliation -- the one place where every composition flag has its
final value. Nothing between the two points touched any of them, which is why
the move is a pure reordering.

`read_base_inp` also matches its keys the way the `input.inp` reader does:
`lbl_match` (anchored, separator-terminated) instead of a bare `index(line,
key) > 0` substring test anywhere in the line. `lbl_match`/`is_sep` moved out
of `input_read`'s internal scope to module scope for that. One file, one
matching rule.

`src/utils/vulcan_to_base.py` now writes the elemental keys it can support:
C, N, O and S, summed over **every** carrier at the handoff level, since
photochemistry moves nuclei between molecules without creating or destroying
them. Mg/Si/Ca/Na/K/Fe are still `metals.inp`'s job -- metal and alkali
chemistry and condensation are outside VULCAN's networks. The same carrier sum
corrects `HeH_base`, which used to be `q_He/(q_H + 2 q_H2)` and so lost the
hydrogen bound in H2O, CH4 and NH3; on the HD 209458 b run that is 0.09698 ->
0.096915, a 0.07% correction and, at 1 microbar, the whole of the difference.

### The gate: do the reservoirs hold?

`src/utils/element_budget.py <run_dir>` answers it from the output profiles.
For every element it compares, cell by cell,

```
n_El / n_H   against   (El/H)_resolved
```

with `n_H` counting hydrogen nuclei in every carrier (`HI`, `HII`, `2 H2`,
`2 H2+`, `3 H3+`, `HeH+`), `n_He` counting `HeI` (which already includes the
triplet), `HeII`, `HeIII`, `HeH+`, and hydrogen itself closing against the mass
density, `rho = mass_per_H * n_H`. The abundances it compares against are the
run's own resolved values: `EXHALE_resolved.out` now carries `mass_per_H`,
`ntot_bc` and one `abundance_<El>` line per element, all of them written to
full double precision so the check has no round-off floor of its own.

Measured on a copy of the `mol_metals` regression case (hot Uranus, molecular
chemistry + solar trace metals) whose `base.inp` was given `O_H_base 9.80e-4`
-- twice the `metals.inp` value -- and `S_H_base 1.32e-5`, an element
`metals.inp` does not list at all:

```
# element budget: <scratch>/mol_metals_elem
# 504 cells, profiles Ion_species.txt, tolerance 1.0e-08 (relative)
elem        reservoir  max |n_El/n_H     at r[Rp]    verdict
               (El/H)     - res|/res
H      rho/mass_per_H      1.982e-04      4.81412       FAIL
H*      triplet-corr.      3.056e-13      1.11481         ok
       NOTE: rho - mass_per_H*n_H equals 4*n(He 2^3S) to 0.1%; the triplet mass is counted twice in calc_rho (HeI already includes it). The H* row removes it.
He       7.930000e-02      1.269e-12      1.11481         ok
C        2.690000e-04      8.132e-13      1.04248         ok
Ca       2.190000e-06      7.915e-13      1.04248         ok
Fe       3.160000e-05      8.177e-13      1.04308         ok
Mg       3.980000e-05      8.310e-13      1.04248         ok
N        6.760000e-05      8.037e-13      1.04308         ok
Na       1.740000e-06      8.091e-13      1.04074         ok
O        9.800000e-04      1.803e-12      1.04553         ok
S        1.320000e-05      8.220e-13      1.04074         ok

1 element(s) miss their reservoir by more than 1.0e-08
```

The nine elemental ratios close to round-off (at most 1.8e-12 relative over
504 cells). Hydrogen does not, and the reason is not the handoff: the
mismatch is exactly `4 n(He 2^3S)`, the triplet mass counted a second time in
`calc_rho`, because `HeI` already contains the triplet (`util_ion_eq` forms
the singlet as `nhei - nheiTR`). Removing it, hydrogen closes to 3.1e-13 as
well -- the `H*` row. The double count is a pre-existing defect of the mass
and particle bookkeeping, not of P2; it is left standing here because fixing
it moves every `He23S: True` golden, which is a change of its own. Its size
on this run is 2.0e-4 of `rho` at 4.8 R_p and zero in the molecular base,
where the triplet population is.

The `EXHALE_resolved.out` of that run shows `abundance_O = 9.80e-4` and
`abundance_S = 1.32e-5`, i.e. the handoff overrode `metals.inp` for oxygen and
activated sulfur.

### Legacy

A `base.inp` without an element key behaves exactly as before: `make check` is
6/6 byte-identical, goldens untouched.

## 75. The He 2^3S double count: an excited level is not a second species (2026-08-27)

Section 74 measured it and left it standing: `rho` carried the mass of every
metastable helium atom twice. The He I density the code transports, `nhei`
(`f_sp` column `isp_HeI`), is the **total** He I population -- the solver sets
it as `nhe*(1 - x_HeII - x_HeIII)` and the photoionization rates recover the
singlet where they need it as `nheiS = nhei - nheiTR` (`util_ion_eq`). The
2^3S column is a **level population inside that number**, not a species beside
it. Every budget sum nevertheless added it a second time with mass 4 m_H, one
gas particle and one helium nucleus.

### Where the second count was

Seven places, all reading the triplet column as if it were an independent
species:

| Place | What was counted twice |
|---|---|
| `calc_rho` (`utilities.f90`) | 4 m_H of mass -> `rho` |
| `calc_ntot` (`utilities.f90`) | one gas particle -> `n_tot`, hence `p = (n_tot + n_e) T` |
| `element_ratio_HeH` (`composition.f90`) | one He nucleus -> the He/H element ratio |
| `element_nucleus_counts` and the `gotHe` shortfall of `project_elements` (`binary_element_diffusion.f90`) | one He nucleus -> the two element totals the diffusion operator conserves |
| `hecar_isp` (`binary_element_diffusion.f90`) | one collision partner -> `carHe`, the helium friction |
| `nHe_l` / `gotHe_l` of `load_IC` | one He nucleus -> the He/H a restart file is judged to carry |
| `Run._NUC_HE` (`examples/exhale_io.py`) | one He nucleus -> `heh_profile` |

`calc_ne` was never wrong: the triplet is neutral, so its `bsp_charge` is zero
and the term it contributed was zero.

### The fix

`species_table` gains `bsp_is_excited_level(n_bsp)`, `.true.` only for `HeTR`,
with the rule written beside it: a flagged column is a sub-population of
another column, so every budget sum skips it while the transport and the
scaling of the level itself do not. The nucleus- and charge-counting loops
(`composition`, `binary_element_diffusion`, `diffusion_tests`) test the flag;
`hecar_isp` loses its `HeTR` row, whose mass, charge and polarizability were a
copy of `HeI`'s anyway; `load_IC` and `exhale_io` drop the triplet term from
their helium sums.

`calc_rho` and `calc_ntot` no longer take `nheiTR` **as an argument at all**.
Removing the dummy is what makes the defect unrepeatable: there is no longer a
triplet number in scope to add. The `thereis_HeITR` branch of
`energy_semi_implicit`'s `calc_ntot` call disappears with it.

### What moved

`make check`: `wasp_he23off` (He 2^3S off) stays byte-identical, as it must --
nothing else changed. The five `He23S: True` cases move, and the size of the
move is the triplet mass fraction of each run:

| case | max abs. rel. drho | max abs. rel. dT | dMdot/Mdot | max 4 n(2^3S)/rho |
|---|---|---|---|---|
| `wasp_full` | 1.24e-06 | 1.32e-07 | -6.5e-08 | 1.29e-06 |
| `wasp_he23off` | 0 | 0 | 0 | 0 |
| `mol_base_handoff` | 4.01e-04 | 1.59e-04 | +3.7e-04 | 1.95e-04 |
| `mol_metals` | 2.57e-04 | 6.77e-05 | +2.1e-04 | 1.98e-04 |
| `mol_lyman_werner` | 1.39e-04 | 1.27e-04 | +4.4e-05 | 2.92e-04 |
| `mol_diffusion` | 5.35e-04 | 2.73e-04 | +5.4e-04 | 1.94e-04 |

WASP-121 b is hot enough that its 2^3S population is a part in 10^6 of the
mass; the hot-Uranus cases carry a few parts in 10^4. Step counts are
unchanged in all six. The goldens of the five moved cases were re-snapshotted
(`check` -> `golden` -> `check`, 6/6 PASS); `wasp_he23off` was left alone.

### The gate

`src/utils/element_budget.py` on the `mol_metals` outputs, before and after,
same `EXHALE_resolved.out`:

```
before:  H   rho/mass_per_H   1.979e-04  at 4.81412   FAIL
after:   H   rho/mass_per_H   2.975e-13  at 1.11481     ok
```

Hydrogen now closes against the mass density at round-off, like the nine
elemental ratios beside it, with no correction of any kind. The `H*`
triplet-corrected row and its NOTE are gone from the script: there is nothing
left to correct.

### Effect on a converged wind

LHS 1140 b, `heh0p55` (He/H = 0.55, He 2^3S on, metals off), a scratch copy
re-finished with JFNK from the same seed under both binaries, then the same
post-processing and transit pass:

| | before | after |
|---|---|---|
| JFNK | `info=0`, \|\|R\|\| = 4.326e-04 | `info=0`, \|\|R\|\| = 4.277e-04 |
| Mdot [g/s] | 6.47929e+07 | 6.48008e+07 (+1.2e-04) |
| He 10830 red depth [%] | 4.279701 | 4.279674 (-6.3e-06) |

The wind is the same wind: `rho` differs by at most 2.5e-04 (at the base,
r = 1.0006), the triplet mass fraction of this run peaks at 2.6e-04 near
11 R_p, and the observable is unmoved at the fifth digit. The stored run
directories were not regenerated.

## 76. The lower atmosphere as a profile, not six scalars (2026-08-27)

Milestone E1 of `docs/phase_e_flux_closure_design.md`. The lower-atmosphere
handoff has been a set of single-level scalars in `base.inp` since Tier 3, and
that shape is what makes the composition an input rather than an output: `HeH`
is one number, the diffusion operator imposes it as the Dirichlet reservoir at
the base, and nothing in the handoff carries a gradient, so no flux condition
between the two models can even be stated. E1 replaces the shape. It does not
yet close the flux; that is E4.

### The file and the key

One new opt-in key,

```
Lower atmosphere profile: lower_atmosphere_profile.dat
```

names a plain-text file: a block of `# key value` header lines, a
`# columns:` schema line, then a fixed-column table running deep to shallow.
Same conventions as every other EXHALE product, so `examples/exhale_io.py`
reads it (`load_lower_atmosphere_profile`) with the loader everything else
uses. The full schema is `docs/input_schema.md` section 2d; the header of the
shipped example is

```
# EXHALE lower-atmosphere profile (docs/phase_e_flux_closure_design.md section 2)
# solution_id 43902c714fe3462149dfe4e5e259553ff3327e405f2a0039176bb5668ee56770
# source_code analytic
# source_version examples/17_lower_profile/make_example_profile.py
# mechanism none (synthetic H2/H/He column, no reaction network)
# stellar_flux none (no photochemistry was solved)
# p_match_bar 9.99999999999999955E-07
# p_top_bar 1.00000000000000002E-08
# p_deep_bar 1.00000000000000002E-03
# trial_flux_H 5.00000000000000000E+10
# trial_flux_He 0.00000000000000000E+00
# iteration 0
# reached_steady_state T
# notes synthetic example column; T(p) prescribed, no climate solution, no cold trap
# units p[bar] r[R_J] T[K] n_tot[cm^-3] rho[g/cm^3] Kzz[cm^2/s] q_*[-] X_*[El/H nuclei] F_*[g/s]
# columns: p r T n_tot rho Kzz q_H2 q_H X_He X_C X_N X_O F_H F_He q_H2O
```

`p_match_bar` is the producer's statement of where the two models meet, and
EXHALE puts its base there: `T0`, `R0`, `p_base_bar`, `q_h2_base`, `HeH` and
every `X_<El>` are the profile's columns interpolated at that pressure, through
the same doors the scalar keys use (`set_element_abundance` for the elements,
direct assignment for the rest). `n_tot` and `rho` are carried, not imposed:
the base density is still `n0` of `input.inp`, and the base particle count and
mass still follow from `T0`, `HeH` and `q_H2` through `comp_ntot_bc` /
`comp_mass_per_H`. Nothing in `composition.f90` changed.

Interpolation is linear in `log p` — the variable both models solve on; the two
codes' radius scales are equal only if their hydrostatic integrations agree —
and is never extrapolated. A target that falls on a level returns that level's
value bit for bit.

Columns the reader has no consumer for are kept, not dropped, and every column
is found by name: a producer's element list is ordered differently from run to
run.

### K_zz stops being a constant

`he_kzz` was one scalar read at three places inside the diffusion operator: the
diffusive time scale that sets the relaxation step, the binary operator's
gradient-plus-eddy face coefficient, and the trace-metal kernel's `D + K`.
`docs/eddy_diffusion_kzz.tex` and `LHS1140b/kzz_decision.md` had already
recorded why a constant is the wrong shape — the molecular coefficient rises by
nearly five decades between the base and 1.2 R_p, and the homopause, the one
consequential thing K_zz does, is where `K_zz = D`, a property of two profiles.

The operator now reads `kzz_cell(1-Ng:N+Ng)`, filled once by
`eddy_diffusion_on_grid` after `define_grid`: from the profile's `Kzz` column
interpolated onto the grid when one is given, and from the scalar `he_kzz` in
every cell otherwise. The two face sites take `0.5*(kzz_cell(j)+kzz_cell(j+1))`,
matching how `Df` and `Gf` are already face-averaged beside them.

`he_kzz` and its key `He_Kzz` keep their names: the constant is what a run
states when it has no profile. With a profile in use the key is inert and a
warning says so.

### The scalar file is no longer a second source

With a profile in use, every `base.inp` key of the EOS-boundary,
elemental-reservoir and boundary-constraint categories is refused with an
`error stop` naming the key and its category; provenance comments and
diagnostic keys are unchanged. This removes the same-solution problem by
construction for those keys instead of checking it.

What is left to check is the pair as a whole, and it is now a refusal rather
than the message it used to be: a `base.inp` beside a profile must carry a
`# solution_id <hash>` comment, and it must equal the profile's. Either
condition failing stops the run with both ids printed.

### What the run records

`EXHALE_resolved.out` gains `lower_profile_present` and, with a profile in
use, the file name, `solution_id`, source and version, `p_match_bar`,
`p_top_bar`, the trial fluxes, the iteration index and whether the profile is
iterable at all (a hand-written profile with no `iteration` index runs; it
just cannot be driven by the closure loop). It is written on the direct steady
route (`EXHALE_PTC=1`) as well now, which it was not before — `element_budget.py`
and the closure both read it, and every LHS 1140 b case uses that route.

`write_element_flux_profile` gained the hydrogen element,
`F_H = 4 pi r^2 (rho (1-X) v - J)`, on the same faces and from the same `X` and
`J` the step used, so `F_H + F_He` is the face mass flux identically. The file
moved from `./diffusion_faceflux.txt` to `./output/element_flux_profile.txt`,
with the other products, and is written unconditionally when a profile is in
use rather than only under `EXHALE_DIFFUSION_CHECK=1`.

### The overlap window is empty, and that is a finding

Section 3.4 of the design asks for one number per element: the elemental flux
over the overlap window, as a median and a relative radial spread, measured on
the escape window `[j_min:N]` intersected with the profile's coverage. That is
implemented and it works — a run with `Escape radius [R_p]: 1.02` reports

```
lower_profile_flux_state  measured
lower_profile_flux_nface  59
```

— but `j_min` is the first cell with `r >= r_esc`, and `r_esc` is 2 R_p in
every configuration that carries diffusion, while a profile reaching from the
microbar match to 1e-8 bar covers about 0.05 R_p above the base. Two decades of
pressure is a few scale heights, not a radius doubling. On every realistic
configuration the two intervals do not intersect and the run reports

```
lower_profile_flux_state  window_empty
```

which is what section 3.4 requires it to do rather than quote a base-cell
number. The consequence is for E4, not for E1: the window the closure needs is
between the base sound-wave region (~1.02 R_p) and the profile top, and its
lower edge has to be stated as such instead of borrowed from `j_min`. The
design document now carries this as an open E4 item.

### Tests

| id | what was run | result |
|---|---|---|
| T-E1 | `examples/17_lower_profile/lower_atmosphere_profile.dat` read by the Fortran reader (`EXHALE_PARSE_DUMP=1`) and by `exhale_io.load_lower_atmosphere_profile`, then re-emitted | 15 columns carried (including the four the code has no consumer for), re-emission bitwise identical to the file, `max abs` round-trip error 0 over every column. The reader's match values equal the loader's: `T0 = 1.450000000000000E+003`, `R0 = 9.794530965504000E+009` (= 1.401 R_J exactly), `HeH = 8.333333300000000E-002`, `q_h2_base = 1.196136002775885E-001`, `p_base_bar = 1.000000000000000E-006` |
| T-E2 | `make check`, no `Lower atmosphere profile:` key anywhere | **`==> REGRESSION PASS (all cases byte-identical)`** -- 6/6, 12 file comparisons, goldens untouched |
| T-E3 | `mol_diffusion` re-run from a scratch copy with `base.inp` deleted, `He_Kzz` removed, and a profile whose `Kzz` column is uniformly `1.0e9`, `EXHALE_MAXSTEPS=12000 OMP_NUM_THREADS=1` | `Hydro_ioniz.txt` and `Ion_species.txt` **byte-identical** to `golden/mol_diffusion`. `K_zz(r) on the grid: 1.000E+09 to 1.000E+09 cm^2/s`; the parse dump differs from the `base.inp` run in exactly one variable, `he_kzz` (0 with the profile, 1e9 from `Kzz_base`), because the operator now reads `kzz_cell` |
| T-E4 | four refusals and one acceptance, each from a scratch run | `T_base` / `O_H_base` / `Kzz_base` beside a profile stop with the key and its category named; a mismatched `solution_id` stops with both ids printed; a `base.inp` with no `solution_id` stops; a `base.inp` carrying only a matching `solution_id` and diagnostic comments runs |
| (schema) | six malformed profiles | each refused with its own message: no overlap (`p_top_bar >= p_match_bar`), missing `solution_id`, missing required column, table not strictly decreasing in `p`, `p_match_bar` outside the coverage, a short data row |
| (unit) | `make diffusion_tests && ./diffusion_tests.x` | 28 passed, 0 failed -- the T7c homopause tests set `he_kzz` and now refill `kzz_cell` through the same routine the code uses |

`make check` was run on the binary built from the state of the tree at the
time; the edits made after that build are confined to the profile module
behind its `if (len_trim(lap_file) .eq. 0) return`, to
`apply_lower_atmosphere_profile` behind `if (.not. lap_in_use) return`, and to
comments, so none of them is on the path the six cases execute.

### Two things fixed on the way

- The face-flux table's `# columns:` line carried labels with spaces in them
  (`F_He=4pi r^2 (F_adv+J)[g/s]`), so it could not be split into one token per
  column the way every other EXHALE product's schema line can. It is now one
  whitespace-free token per column and the formulas moved to the comment lines
  above it.
- `LHS1140b/make_memo_figures.py` and `LHS1140b/exhale/kzz_scan_table.py` read
  `D_eff` out of that table by fixed position (`usecols=(1, 8)`), which the new
  `F_H` column would have shifted silently. Both now find the file under either
  name and the column by name, falling back to position 8 for a table written
  before the schema line was fixed. Stored run outputs were not regenerated.

---

## 77. Two producers for the profile handoff: the Photochem adapter and its VULCAN cross-check (2026-08-27)

Milestone E2 of `docs/phase_e_flux_closure_design.md`. Section 76 gave EXHALE
a reader for the lower-atmosphere profile and an example file built by hand;
what was missing was a producer. E2 supplies two, behind one command line:

- `src/utils/photochem_to_lower_profile.py` — solves a gas-giant
  photochemistry model to steady state on a prescribed `T(p)`, `K_zz(p)`
  column and writes the handoff. This is the production path.
- `src/utils/vulcan_to_lower_profile.py` — the same schema and the same
  fingerprint from a finished VULCAN `.vul` output. This is the cross-check
  arm, not a second production path: VULCAN has no climate model, so it can
  never be the route away from a prescribed `T(p)`
  (`docs/vulcan_photochem_comparison.md` P1.7).
- `src/utils/lower_profile_schema.py` — everything that is a property of the
  FORMAT rather than of either chemistry code: the shared command line, the
  elemental accounting, the hydrostatic radius, the level insertion that puts
  an exact node at the matching pressure, the `solution_id` and the two files
  that carry it. It exists so the two adapters cannot drift apart; a schema
  written twice is a schema that will be written differently twice.

The scalar generator `src/utils/vulcan_to_base.py` stays. It writes the
single-level `base.inp` special case that four regression cases pin, and the
two must not be read into the same run: with a profile in use EXHALE refuses
every scalar physics key beside it (section 76).

### One way, and the file says so

E2 is the one-way milestone. The trial elemental fluxes are recorded in the
header and in the fingerprint, and they are **not** imposed as an upper
boundary condition on the chemistry; imposing them, measuring the wind's own
elemental flux over the overlap and iterating the pair is E4. The `F_H` and
`F_He` columns therefore carry the stated trial values, and the header `notes`
line says so in the file itself, so that no reader can mistake a boundary
value for a measurement:

```
# notes one-way E2 handoff: F_H and F_He are the STATED trial fluxes, not
        measured; no climate step, T(p) prescribed.
```

### What the fingerprint is over

`solution_id` is a sha256 over the lower model's configuration: the source
code and version, the mechanism and thermodynamic files by their own sha256,
the stellar flux file by its sha256 and the dilution applied, the prescribed
`T(p)` file, the planet mass and reference radius, the zenith angle, the
top-of-atmosphere pressure, the elemental abundance vector, the matching
pressure and the two trial fluxes. The closure *iteration index* is
deliberately outside it: two iterations handed the same configuration and the
same trial fluxes are the same solution, and the fingerprint is what says so.
The same id is written into the paired `base.inp`, which carries **no physics
key at all** — the provenance comment is its whole content, and any scalar
key in it would stop the run by design.

### Elemental accounting, and the trap that makes it a rule

Element columns are El/H **nuclei** ratios summed over every carrier, the
lesson `vulcan_to_base.py` already records: photochemistry moves nuclei
between molecules without creating or destroying them, so the nuclei ratio is
the conserved quantity the wind inherits and the only one `melem_ab` can hold.

How the carriers are counted differs by code, and it has to. VULCAN labels are
formulas once the ionization sign and the excited-state suffix are stripped
(`H3+`, `CH2_1`), so they are parsed. Photochem labels are **not**: `O1D` and
`N2D` are excited states, and a formula parser reads them as an oxygen with a
deuterium attached. The Photochem arm therefore takes the composition from the
model's own matrix (`dat.species_composition`), indexed by name against
`dat.atoms_names` — which is also the trap P1.8 item 4 records, since that
atom order changes from run to run. Condensed carriers are excluded by
default: the gas EXHALE's base inherits is the gas phase, and the difference
between the two sums is the cold trap, which is E3's subject.

The adapter's own conservation check (design section 4.2) is stated at the
deepest level, because that is the only level where the elemental ratios are
still the input ones — transport and the upper boundary move nuclei with
height, which is the whole point of the closure. The design proposed a
tolerance of 1e-10. Measured on HD 209458 b, the departure is **1.6e-7** with
the Zahnle H/He/N/O/C set and **2.7e-5** with the VULCAN NCHO network
converted through `vulcan2yaml`, i.e. chemistry-solver noise and not a
miscount, so the default is 1e-4: loose enough not to refuse a converged
solution, tight enough that a miscounted carrier (an O(1) or percent-level
error) is still caught. The measured worst departure is printed on every run.

### The radius column

`r` is the adapter's own hydrostatic integration from a stated reference level
(`--r-ref` at `--p-ref`, 1 bar by default, the transit-radius convention),
RK4 in `ln p` with radius-dependent gravity, refined 64x per level so the
answer does not depend on how the chemistry code spaced its grid. The mean
molecular weight enters as `mu * m_amu`, the mass of a particle of `mu`
atomic mass units.

### Trap 1 was met and is handled the way P1 handles it

`gasgiants` requires `3*TOA_pressure_avg < P_climate_top`. On the HD 209458 b
column the supplied `T(p)` stops at 2.152e-2 dyn/cm^2 while the requested TOA
was 1e-2, so the requirement fails on the file as given. The adapter keeps the
stated TOA and cuts the climate column back to `3.05*TOA`, which is what
`read_climate` does in `docs/p1_matched_comparison.py`; the alternative --
silently lowering the TOA to fit the file — was tried first and rejected,
because it changes a stated input without saying so. Trap 2 (the
thermodynamic data capping the temperature) is handled by **truncating** the
column from the top down, never by clipping `T`, because a clipped `T(p)` is a
profile nobody solved; the truncation is recorded in `notes`, and if it falls
below `p_match_bar` the run is refused. HD 209458 b needed no truncation: its
5891 K top was accepted.

### Failure policy

No profile is written when the handoff would not be one. `reached_steady_state
= False` or `gave_up = True`, a temperature truncation that reaches below the
match, a climate cut that leaves fewer than ten levels, `p_top_bar >=
p_match_bar`, a matching pressure outside the table, fewer than two levels
after the requested window, and the elemental conservation check are all
refusals that print the reason and write nothing.

### Measurements

HD 209458 b, the P1 configuration (`M_p = 0.720 M_J`, `r(1 bar) = 1.36 R_J`,
`R_star = 1.155 R_sun`, `a = 0.0480 AU`, Gueymard solar spectrum diluted,
zenith 48 deg, TOA 1e-2 dyn/cm^2, VULCAN's own `atm_HD209_Kzz.txt` as the
prescribed `T(p)`, `K_zz(p)`), matching level 1 microbar. Both Photochem arms
reached steady state; the VULCAN arm is its stored HD209 output.

| at `p_match = 1e-6` bar | VULCAN / NCHO | Photochem / NCHO | Photochem / Zahnle H,He,N,O,C | PC/VUL, same network |
|---|---|---|---|---|
| `T` [K] | 2332.25 | 2332.28 | 2332.43 | 1.0000 |
| `r` [R_J] | 1.495029 | 1.494127 | 1.494022 | 0.9994 |
| `n_tot` [cm^-3] | 3.1098e12 | 3.1095e12 | 3.1124e12 | 0.9999 |
| `rho` [g cm^-3] | 9.4403e-12 | 9.6497e-12 | 9.7202e-12 | 1.0222 |
| `Kzz` [cm^2 s^-1] | 5.000e11 | 5.000e11 | 5.000e11 | 1.0000 |
| `q_H2` | 0.42366 | 0.45363 | 0.46266 | 1.0707 |
| `q_H` | 0.44987 | 0.41710 | 0.40726 | 0.9271 |
| `X_He` | 9.6915e-2 | 9.6912e-2 | 9.6912e-2 | 1.0000 |
| `X_C` | 2.7745e-4 | 2.7754e-4 | 2.7754e-4 | 1.0003 |
| `X_N` | 8.1809e-5 | 8.1831e-5 | 8.1828e-5 | 1.0003 |
| `X_O` | 5.0329e-4 | 6.0608e-4 | 6.0608e-4 | 1.2042 |
| `q_H2O` | 2.2073e-4 | 3.4312e-4 | 3.5306e-4 | 1.5545 |

The matched-network hydrogen ratio is `q_H(VULCAN)/q_H(Photochem) = 1.079`,
which reproduces the **1.08** P1 measured for the code factor on this planet's
2331 K base (`docs/vulcan_photochem_comparison.md` P1.7) — the adapters carry
the same disagreement the direct comparison found, and add none of their own.
`T`, `n_tot` and the two element ratios the codes were given identically agree
to 1e-4 or better.

`X_O` and `q_H2O` differ by 1.20x and 1.55x for a reason that is not a code
difference: the stored VULCAN HD209 run was configured with `O/H = 5.03e-4`,
against the Lodders 2009 `6.06e-4` the Photochem arms were given. That is a
difference of input, and it is what the `X_<El>` columns are for — it would
have been invisible in a comparison of molecular mixing ratios alone.

### Fingerprint round trip

| test | result |
|---|---|
| the same input twice, VULCAN arm | same `solution_id` (`4f80f048...`), and the profile files are **byte-identical** |
| the same input twice, Photochem arm (two full chemistry solves) | same `solution_id` (`f6232692...`), and the profile files are **byte-identical**: the chemistry solve reproduces, so the fingerprint and the table agree |
| `--trial-flux-H 5e10` instead of 0 | `0a4d2051...`, different |
| `--p-match 2e-6` instead of 1e-6 | `ceed0b58...`, different |
| `--iteration 3` instead of 0 | unchanged, by design: the iteration index is bookkeeping, not configuration |
| the paired `base.inp` | carries the same id as the profile, and no physics key |

### EXHALE reads what the adapter wrote

A scratch case built the way `examples/17_lower_profile/` is (the same
`input.inp`, the Photochem HD 209458 b profile and its paired `base.inp` in
place of the synthetic ones), run single-threaded with
`EXHALE_PARSE_DUMP=1` and then with `EXHALE_MAXSTEPS=20`:

```
(lower_atmosphere_profile) Reading lower_atmosphere_profile.dat ..
   base.inp: solution_id matches the profile.
   profile: T0 ->    2332.4 K          profile: q_H2(base) ->  0.46266
   profile: R0 ->   1.4940 R_J         profile: He/H ->  0.09691
   profile: p_base ->  1.00E-06 bar    profile: C/H ->  2.775E-04
   profile: source photochem, closure iteration 0
(lower_atmosphere_profile) K_zz(r) on the grid:  5.000E+11 to  5.000E+11 cm^2/s
```

Every echoed value is the profile's own number at the matching level
(`T = 2332.43 K`, `r = 1.494022 R_J`, `q_H2 = 4.626610e-1`,
`X_He = 9.691225e-2`, `X_C = 2.775366e-4`), and `EXHALE_resolved.out` carries
the provenance: `lower_profile_solution_id f6232692...`, `source_code
photochem`, `p_top_bar 7.44e-9`, `iterable T`, `steady T`. The paired
`base.inp` — provenance comments only — was accepted, and the id check passed
rather than being skipped.

The run itself is a 20-step relaxation snapshot and its `Mdot` is not a
result. One thing it does show: `lower_profile_flux_state window_empty`, the
same finding E1 recorded — with `Escape radius [R_p]: 2.00` the escape window
`[j_min:N]` and the profile's coverage do not intersect, so the elemental flux
has no window to be measured on. That is the open item E4 inherits, restated
here on a real photochemical profile rather than a synthetic one.

### What is not done

The LHS 1140 b application and T-E8 (`element_budget.py` closing on a profile
run) are E2's stated gate and neither is done here; the adapters were
exercised on HD 209458 b because that is where P1's measurement exists to
check them against. The climate step, the tropopause and the cold-trap water
are E3, and the flux closure is E4.

---

## 78. The climate step, the cold trap, and LHS 1140 b through the profile (2026-08-27)

Milestone E3 of `docs/phase_e_flux_closure_design.md`, and with it the two
gates E2 left open: the LHS 1140 b application and T-E8.

Until now the handoff's `T(p)` was an input. `--climate` makes it a solution:

```
<photochem 0.8.4 env>/bin/python \
    src/utils/photochem_to_lower_profile.py LHS1140b/lower_profile \
    --mp 0.0176220 --r-ref 0.157692 --p-match 1e-6 --kzz-const 1.0e9 \
    --stellar-flux LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt \
    --flux-at-planet --wavelength-unit A \
    --atoms H,He,N,O,C --abundances He=2.09 \
    --climate --climate-p-deep 20.0 --toa 1.0e-2 --boa-pressure-factor 1.0
```

`src/utils/radiative_convective_column.py` is the new module: Photochem's
`clima` (`AdiabatClimate`) with `solve_for_T_trop` on, so the stratospheric
temperature is the skin temperature of the solution and not a number chosen
for it, and with the background gas picked by abundance rather than by name —
which is what lets one code path cover both a hydrogen atmosphere and the
helium-rich one the LHS 1140 b line implies. The elemental vector is
partitioned into the carriers a cool H/He atmosphere holds in equilibrium (O
in H2O, C in CH4, N in N2, He atomic, the rest H2) with the validity range
written at the code and a refusal above 1000 K, because that partition is the
equilibrium one only below the CO/CH4 and N2/NH3 transitions. Nothing
downstream inherits it: the photochemistry computes its own composition from
the same elemental vector.

### What the codes actually do about condensation, measured

The design proposed enabling the mechanism's H2O particle and setting rainout
to the tropopause. Both halves turned out wrong for a gas giant, in opposite
directions, and both were measured rather than assumed:

- the H2O particle needs no enabling. The Zahnle set restricted to
  H/He/N/O/C already carries `H2Oaer`, and asking for
  `water-condensation: true` is **refused** — *"Either
  "fix-water-in-troposphere" or "water-condensation" is turned on in the
  settings file, but the reaction mechanism already implements H2O
  condensation via a particle."* The cold trap is applied by the chemistry
  solver already; the settings switch must stay `false`.
- rainout is a surface process: `gas-rainout: true` is refused for want of a
  `rainfall-rate`, and a gas-giant model has no surface for rain to fall to.
  It is not how a hydrogen or helium atmosphere traps its water.

Section 5 of the design now records this correction beside the procedure.

### The tropopause is solid, its pressure is not

Scanned over the one free parameter of the climate solve, its deep boundary:

| deep boundary [bar] | `T_deep` [K] | `P_trop` [bar] | `T_trop` [K] | `f_H2O` above the tropopause |
|---|---|---|---|---|
| 1 | 242.4 | 0.488 | 186.0 | 3.33e-7 |
| 3 | 308.5 | 0.765 | 185.8 | 2.06e-7 |
| 10 | 439.9 | 0.974 | 185.4 | 1.49e-7 |
| 20 | 555.9 | 1.031 | 185.0 | 1.32e-7 |
| 30, 50 | refused: the pseudoadiabat leaves the thermodynamic data | | | |

`T_trop` moves by 1 K over a factor of 20 in the deep boundary, against the
skin temperature `T_eq/2^(1/4) = 190 K`. `P_trop` moves by a factor of two,
and the trapped water moves with it, since at fixed `T_trop` the trapped
mixing ratio is `p_sat(T_trop)/P_trop`.

**Against Cherubim et al. (2026).** Their closed-form estimate is `f_H2O ~ 7
ppm` at a 0.1 bar tropopause with `T_skin = 194 K`. The climate solution gives
**0.13–0.33 ppm at a tropopause of 0.5–1.0 bar** — a factor 20–50 less water
at a tropopause 5–10x deeper. The saturation physics agrees; what differs is
where the tropopause is put, and ours is a solved level in a helium-dominated
column (`mu = 3.6` at `He/H = 2.09`) rather than an assumed one. Reported as a
result about the climate model, not tuned toward the published number.

### What the trap does to the handoff

Gas-phase oxygen falls from `O/H = 6.062e-4` at the deep boundary to
**4.957e-7** at the matching level, a factor of 1223, while `C/H`, `N/H` and
`He/H` cross the tropopause unchanged: CH4 and N2 do not condense at 185 K.
The profile now always computes the elemental sums twice, over the gas phase
and over every carrier including the condensates, and carries both in the
header `notes` — at the match the two coincide, because the condensate is left
behind at the tropopause. In this handoff the cold trap shows as a depletion
with height, not as a condensed reservoir sitting at the base.

The 20 bar deep boundary is not a physical choice: it is the shallowest one at
which the photochemical model bottom fits inside the climate column. On an
atmosphere this cold the chemistry never equilibrates within any column
`clima` can reach, so `determine_quench_levels` returns the deepest level of
the grid and Photochem's default `BOA_pressure_factor = 5` then asks for a
bottom five times deeper than the column it was handed.
`--boa-pressure-factor` exposes it and the adapter's refusal text names it.

The prescribed-column path of E2 is unchanged where it matters: re-run on the
HD 209458 b configuration of section 77 it reproduces that section's numbers
exactly (`T = 2332.43 K`, `r = 1.494022 R_J`, `q_H2 = 4.626610e-1`,
`X_He = 9.691225e-2`). Its `solution_id` does change, because
`boa_pressure_factor` is now part of the configuration the fingerprint is
over, which is the fingerprint doing its job: a profile made with a different
model bottom is a different solution. Files already on disk are untouched.

### LHS 1140 b through the profile: it converges, and T-E8 closes

`LHS1140b/lower_profile/lower_atmosphere_profile.dat` (102 levels, 16.73 bar
to 1.10e-8 bar, `solution_id 51302c76...`) drives a direct-steady JFNK run,
seeded from `exhale/heh2_diff_kzz1e9`, to `info = 0` with the element-
diffusion outer loop converged (composition drift 9.1e-4 at pass 8). What the
reader put in agrees digit for digit with what the file says at the match, and
`EXHALE_resolved.out` carries it:

```
planet_radius_RJ           1.623994484710662E-001
equilibrium_temperature_K  1.850492232658783E+002
HeH_number_ratio           2.092134697681555E+000
abundance_C                2.777529772010253E-004
abundance_O                4.956907073555432E-007
abundance_N                8.187442081270519E-005
lower_profile_solution_id  51302c7647a319bdfbd142d3b4bf0a367b25d2dc07e2675c8da3cd430dc63177
```

`lower_profile_notes` is now echoed there too, so the climate record — deep
boundary, tropopause, cold-trap water, both elemental sums — travels with the
run and not only with the file (`lap_notes` widened to 1000 characters to hold
it).

**T-E8.** `src/utils/element_budget.py` closes on every element the profile
carries: `C 1.44e-14`, `N 1.47e-14`, `O 1.45e-14`, against a `1e-8` tolerance,
and on H (`0.0`) and He (`2.1e-16`) at the base. *All element budgets close
within 1.0e-8.*

The H and He rows needed a fix in the checker, not in the run. With
`He_diffusion` on, He/H is a *solved profile* — 2.0921 at the base, 1.206 at
1.15 R_p, 0.167 at 30 R_p — and the resolved `HeH_number_ratio` is the
reservoir the base is held at, not a column invariant. The tool was testing it
cell by cell and calling the separation the operator exists to produce a
failure. `EXHALE_resolved.out` now states `he_diffusion`, and with it on the H
and He rows are checked at the base cell, with the solved separation printed
beside them. The trace elements are unaffected and stay column-wide: nothing
diffuses them.

The Fortran side of E3 is two writes and one string length — `lap_notes`
widened to 1000 characters, echoed into `EXHALE_resolved.out`, and
`he_diffusion` stated there — none of which any solver reads. `make check` is
**6/6 byte-identical**, goldens untouched.

### Profile route against the scalar base, and where the difference lives

Beside `exhale/heh2p13_diff_kzz1e9`, the converged scalar-base case at the
same eddy coefficient. Every row below is a JFNK `info = 0` run finished the
same way, `red` being the He 10830 red-pair depth from
`he_line_metrics.py`:

| case | what it adds | log10 Mdot | red |
|---|---|---|---|
| `heh2p13_diff_kzz1e9` | scalar base, `He/H = 2.13`, `T0 = 226 K` | 7.80 | 4.396 |
| control A | `HeH_base = 2.0921` | 7.81 | 4.393 |
| control B | + `T_base = 185.05 K` | 7.81 | 4.389 |
| control C | + `r_base = 0.162399`, `q_H2_base = 0.19265` | 7.86 | 4.985 |
| control D | + `C/N/O_H_base` from the profile | **7.46** | **1.517** |
| profile route | the profile itself | **7.46** | **1.517** |

Two things fall out. Control D and the profile route agree, which is the
consistency statement the single-source rule of section 2.4 is supposed to buy:
handed the same base state, the two routes are the same run. (*They agreed to
every digit printed here. Re-solved on the section-87 coefficient the agreement
is 0.01 dex in `log10 Mdot` and 1 per cent in the line -- inside the workflow's
own 0.3-1.5 per cent path spread, and the two rungs were solved at different
residual tolerances. Agreement to every printed digit was the tolerance, not the
physics.*) And the difference from the scalar LHS 1140 b case is **not** the
base temperature — dropping `T0` from 226 K to the solved 185.05 K moves
`log10 Mdot` by less than 0.01 dex — but the **elemental carbon, nitrogen and
oxygen the profile brings with it**, which the scalar runs never had: −0.40
dex in `Mdot` and a 3.3x shallower helium line, metal-line cooling of a
185 K base at `C/H = 2.8e-4`. *Re-measured on the Garcia Munoz rate of section
87 (`LHS1140b/exhale/basemetals_gm25/`), that step is `-0.37` dex and a factor
`2.89` in the red-pair depth (`R = 68,000`; `2.88` at the `R = 80,000` these
were read at). Note also that the two numbers here belong to different pairs of
runs: `-0.40` dex is the `C -> D` rung, while the factor quoted as `2.9`
elsewhere in the memo is `ms -> D`, the whole ladder, which is `2.59` on the
new coefficient.* `p_base` is 1 microbar on both sides and is not
a candidate.

### What is not done

`lower_profile_flux_state window_empty` again: with `Escape radius [R_p]:
2.00` the escape window and the profile's coverage do not intersect, so T-E5
still has no window to be measured on. That is E4's, together with the flux
iteration itself. The composition here remains an input (`He/H = 2.09`, the
diffusion-limit crossing of `LHS1140b/kzz_decision.md` section 6), not yet an
output.

### Beside it: a particle is not a hydrogen atom (`vulcan_to_base.py`)

Unrelated to E3 and found while reading the hydrostatic integrations.
`src/utils/vulcan_to_base.py` integrated `dr/dln p = -kT/(mu m_H g)` with the
hydrogen ATOM mass, while its `mu_profile` is a mean molecular weight in
**amu**, built from VULCAN's own species masses (H2 = 2.016, He = 4.003,
H2O = 18.02). The mass of one particle is `mu * m_amu`; using `m_H` makes
every particle 0.8% heavier and the scale height 0.8% shorter. Measured on the
HD 209458 b VULCAN output, `r_base` moves from 1.49388 to **1.49504 R_J** —
0.078% of the planet radius, 0.86% of the column above 1 bar, which is the
0.8% arriving where it should. Fixed; `lower_profile_schema.py` already used
`m_amu` and recorded the discrepancy rather than copying it.

The pinned `base.inp` of the regression cases were deliberately **not**
regenerated: they are frozen inputs, and the generator's fix does not
retroactively change a run that was made with the old file.

`run_lower.py` and `lower_column.f90` multiply the same `m_H` by a mean
molecular weight built from integer mass numbers (`mass = 1 + 4 f_He`), where
`m_H = 1.008 amu` is the better of the two constants for the H2 that dominates
that mixture. They are left alone; the mixed convention is noted here so the
next reader does not "fix" them into being worse.

## 79. The elemental-flux closure: the composition becomes an output (2026-08-27)

Milestone E4 of `docs/phase_e_flux_closure_design.md`, and with it the last
three of Phase E's tests, T-E5 through T-E7.

Until now the lower model was told what to do and never heard back. E2 and E3
carried `trial_flux_H` and `trial_flux_He` in the profile header but only as
statements: nothing imposed them on the chemistry, so the handoff was one-way
and `He/H` was still an input the wind reproduced. E4 closes the loop in both
directions.

### The trial flux now reaches the chemistry

`photochem_to_lower_profile.py` imposes each element's trial flux at the
photochemical model top with `set_upper_bc(species, 'flux', flux=...)`, the
elemental flux converted at the model's own top radius,

```
phi_El [nuclei/cm^2/s] = Phi_El [g/s] / m_El / (4 pi r_top^2)
```

and placed on the dominant carrier: helium on `He`, hydrogen on `H2` at
`phi_H/2`, two nuclei to the molecule. The dominant-carrier step is measured
rather than assumed — the atomic-H share of the hydrogen nuclei at the model
top is 1.0e-4 in the converged LHS 1140 b solution — and the adapter refuses
a solution in which that share has risen past 1 per cent, because the flux
would then have been put on the wrong carrier. The sign convention was read
off Photochem's own right-hand side (`photochem_evoatmosphere_rhs.f90`,
`rhs(k) = ... - var%upper_flux/dz`): positive is a loss at the model top,
which is this code's outward. A zero trial flux does not call `set_upper_bc`
at all, so E2's and E3's one-way profiles stay reproducible; that was checked
by reproducing `solution_id 51302c76...` bitwise.

That the condition takes is checked every iteration instead of trusted. The
adapter measures the model's own top-of-atmosphere elemental fluxes with
`gas_fluxes()` and writes them into the header as `measured_flux_H` and
`measured_flux_He`. On the reference arm the imposed 1.800000e7 g/s comes
back as 1.800162e7 g/s, a departure of 9e-5; the closed-top solution's own
numerical floor is 3.0 g/s in hydrogen and 18.2 g/s in helium, seven decades
below.

`--p-top-bar` states the shallowest level to carry directly in bar, and
refuses a top at or above `p_match_bar` (no overlap) or a simultaneous
`--toa`, which would be the same level said twice in two units.

### The window the closure is measured on, and why the obvious one fails

E1 left an open item: the escape window `[j_min:N]` and the profile's
coverage never intersect, so `lower_profile_flux_state` was `window_empty` on
every configuration. E4 replaced the single window with two, and the choice
between them is a measurement.

**The overlap window** is the interval both models describe. Its upper edge
is `r(p_top)`; its lower edge is now measured — the first face at which the
5-face moving spread of `r^2 rho v` falls below 10 per cent, which is where
the standing base sound wave stops dominating — with a 1.02 R_p fallback that
`lower_profile_flux_r_lo_source` reports when no face qualifies. Both edges,
the face count and the reduction go into `EXHALE_resolved.out`.

**On LHS 1140 b that window cannot carry the closure, and extending the
profile is not the fix.** Where it is empty it misses by 1.6e-4 R_p
(`r_lo = 1.009 9998`, `r_hi = 1.009 840`); where it is not empty it holds
12-18 faces between 1.007 and 1.0098 R_p and the flux across them has a
radial spread of **13 in F_H and 35 in F_He**, with a median of 1.4-1.8e8 g/s
against the 3.16e7 g/s the wind actually removes. The reason is geometric:
from the match to the profile top the column covers 1.96 decades of pressure
in 0.0098 R_p, about 0.005 R_p per decade, while the wind is still 6x above
its far-field flux at 1.05 R_p and does not settle within 10 per cent until
~1.5 R_p. Reaching 1.5 R_p on that scale would ask the photochemical column
for of order 100 further decades of pressure. Mapping `p_top` through
EXHALE's own `r(p)` rather than the profile's — the two hydrostatic scales
differ, the wind being hot where the cold column is not — moves the upper
edge to 1.0734 R_p and makes the window non-empty, but its spread is then 5.6
to 7.3. Both mappings say the same thing: the interval the two models share
lies wholly inside the region where the wind's mass flux is not yet its own.

**The steady-flux window**, `r >= r_esc`, is what the closure uses. At a
steady state the elemental flux `4 pi r^2 (rho X v + J)` does not depend on
radius, so a flux measured there *is* the flux through the matching level;
that identity, and not convenience, is what makes it a statement about the
handoff. On the reference arm it holds 191 faces and gives
`F_H = 1.8143e7 g/s`, `F_He = 1.3204e7 g/s` with radial spreads of 0.44 and
0.46 per cent, flat to the same tolerance as the mass flux itself (0.45 per
cent). That is T-E5, met in the window where it can be met. The driver reads
the overlap window first every iteration and substitutes the steady one only
with the reason written into `closure.log` and the window named in the
history table; it never substitutes silently.

The spread depends on how far the wind was converged, and that was measured
rather than assumed: at `Resid tol: 1.0e-3`, the value `finish_case.sh` uses,
the same solution's steady-window spread is **7.4 per cent** (F_H 7.1, F_He
7.9) and would fail the closure's own validity test; continuing the same
solution to `1.0e-4` brings it to **0.50 per cent** and moves `log10 Mdot`
from 7.47 to 7.50. The 7 per cent was under-convergence, not a floor, and the
closure runs the wind at `1.0e-4`.

### The driver

`src/utils/element_flux_closure.py` owns the loop and nothing else. Damped
Picard on both elements with `omega = 0.5`, halved to a floor of 0.125 on any
iteration whose residual fails to fall, the convergence test reading the
**undamped** residual `|F_measured - Phi_trial| / |F_measured|`, `tol = 0.05`
and `k_max = 8`. Both hydrogen and helium are iterated: holding helium at
zero, as the design first proposed, would leave the H:He partition — the
whole question on this planet — as a diagnostic rather than a condition.

Section 6.2's validity rule is a statement about the *window*, not the
iterate: a window whose spread exceeds `tol` cannot decide a residual at
`tol` whatever that residual is, so `spread > tol` stops the run as
UNRESOLVED with both spreads printed, every iteration and independently of
how far the iterate is.

Everything that belongs to the planet rather than the iteration lives in one
JSON configuration (`--print-config-template`), so the driver carries no
planet constants; `LHS1140b/exhale/flux_closure/lhs1140b_closure.json` is the
LHS 1140 b one. `closure_history.txt` carries one row per iteration with both
trial and measured fluxes, both residuals, `omega`, the window used and its
spreads, `He/H` at the match, `log10 Mdot`, the EXHALE `info`, and the
profile's `solution_id` — the fingerprint chain, so the sequence is auditable
after the fact. `--resume` continues from the last completed iteration and
`--dry-run` exercises the reader against an existing run directory.

Two departures from `finish_case.sh`, both measured. `Resid tol: 1.0e-4`, for
the reason above. And the post-processing pass runs **without** the
`EXHALE_PTC*` variables: `finish_case.sh` sets them inline on the JFNK
command alone, and carried into the second pass they make the binary re-enter
the steady solver and stop before writing the advection-corrected profiles,
so `*_adv.txt` and the mass-loss rate never appear. The first LHS 1140 b
closure run reported `log10 Mdot = nan` and that is what it was.

### LHS 1140 b, closed

Three arms, from `Phi_ref` = (1.8e7, 1.3e7) g/s in (H, He) and from `0.3 x`
and `3 x` that, all seeded from the same converged wind, all with
`He/H = 2.09` as the *starting* reservoir.

| start | k | `F_H` [g/s] | `F_He` [g/s] | `He/H` at the match | `log10 Mdot` | He 10830 red depth [%] | EW [mA] |
|---|---|---|---|---|---|---|---|
| `1.0 x` | 0 | 1.81433e7 | 1.32040e7 | 2.0923516 | 7.500 | 1.699 | 4.743 |
| `0.3 x` | 5 | 1.82031e7 | 1.34266e7 | 2.0923486 | 7.500 | 1.721 | 4.804 |
| `3.0 x` | 6 | 1.82119e7 | 1.34710e7 | 2.0923602 | 7.500 | 1.725 | 4.818 |

The three converged `He/H` agree to 5.5e-6, `log10 Mdot` to better than the
0.005 the log prints, `F_H` to 3.8e-3 and `F_He` to 2.0e-2 — all far inside
the 5 per cent criterion, so **T-E7 passes and the closure is single-valued
on this planet.** `omega` never halved in any arm: the residual fell on every
iteration, so the map is a plain contraction here and the damping rule was
not exercised. The iteration counts are exactly what a contraction at
`omega = 0.5` predicts, the residual halving each step, which also says the
3x arm used six of its eight allowed iterations with no margin: a start
further than an order of magnitude out would need a larger `k_max` or a
secant acceleration.

**The result to read, though, is what the closure returns.** `He/H` at the
match comes back at 2.09235 against the 2.09 the reservoir was started from —
a change of 1.1e-3, and only 2e-4 of that from the escape flux itself. The
composition is now an output, and on this planet the output is the
*well-mixed* value: at `K_zz = 1e9 cm^2/s` eddy mixing homogenizes the column
all the way to 1 microbar, so a 1.8e7 g/s hydrogen escape imposed at the
model top moves the matching-level `He/H` by 2e-4 while moving it by 1.8e-2
at the model top itself. It is the eddy coefficient, not the escape flux,
that sets the composition at the match. That is a result about this
configuration and not a general one, and it makes `LHS1140b/kzz_decision.md`
section 0's adopted `K_zz` the quantity the answer now rests on.

Above the match the wind fractionates strongly, as it did before: the
elemental fluxes carry `F_He/F_H = 0.73` by mass against the base reservoir's
8.31, and `He/H` falls from 2.09 at the base to 0.183 at 30 R_p. The closure
does not weaken that; it says the separation happens in the wind and not
below the match.

**The flux budget.** `F_H + F_He` reproduces the mass flux identically —
exactly, in the two arms whose medians fall on the same face — because the
binary operator's `X` is the helium mass fraction of the H+He mixture and the
trace metals ride inside `rho`. So the measured `F_H` is hydrogen *plus*
metals, overstating elemental hydrogen by the metal mass share, which is
4.79e-4 here (`mass_per_H = 9.3739 amu`, C/N/O at 2.778e-4, 8.188e-5,
4.955e-7). That is two decades below `tol` and is stated rather than
corrected.

### T-E8, and one thing it caught

The element budget closes for helium at the base cell, where the reservoir is
the boundary condition, to 1.5e-15. Hydrogen, carbon, nitrogen and oxygen
miss by 5.6e-8, 1.2e-4, 1.2e-4 and 3.3e-4 against the checker's 1e-8, and the
miss is not noise: **it equals the change in the elemental reservoir between
the seeding solution and the run.** `load_IC` restores every metal density
the restart file carries and rebuilds an element from `melem_ab` only when
the file does not carry it at all (`load_IC.f90:283-334`), so an iteration
that changes a reservoir keeps the seed's value everywhere above the base and
`EXHALE_resolved.out` then reports a reservoir the run did not use above the
base cell. It was seen at full size first: an iteration seeded across a 22 per
cent change in `O/H` reported exactly 17.8 per cent, which is that change
expressed on the new reservoir.

This was **not** fixed at the time, because the fix changes restart semantics,
which is a user-visible behavior and expensive to reverse; it was recorded for
a decision. Its size on the result reported here is small and measured: the
closure moves reservoirs by at most 5e-4 between iterations, and the coolant
that matters here is C I, which carries 94 per cent of the radiated energy and
is stale by 1.2e-4; oxygen, the worst row at 3.3e-4, carries under 1e-4 of the
cooling on this planet. The decision was taken and the fix implemented on
2026-08-27: section 80. The numbers quoted in this section are what the code
gave **before** it, and are left as they were measured.

### T-E9: not verified, and why

**Not verified: VULCAN did not converge on this cold He-rich column (longdy
rising 2.65 -> 2.75 against `yconv_cri = 0.01`, `dt` pinned at 1.85e3 s, the
model clock at 1.42e7 s, limiting cell `nz = 72 and H2O` -- the cold trap --
exiting on "Maximal allowed steps exceeded (10000)"); stopped 2026-08-27.**

The cross-check arm was built on the production profile's own `T(p)` and
`K_zz` and never produced a profile, so there is no VULCAN composition at the
match, no wind run on one, and no cross-code difference to state beside the P1
spread. `LHS1140b/exhale/flux_closure/vulcan_arm/` holds only an `input.inp`
and the seeded initial condition; the binary was not launched. The full record,
including the second run that reached step ~500 of 30000 before being stopped
and therefore says nothing either way, is
`LHS1140b/lower_profile_vulcan/vulcan_work/STOPPED.txt`.

Two things found on the way that outlive the attempt. First, a real defect in
`src/utils/vulcan_to_lower_profile.py`, now fixed: `sch.formula_elements`
strips the phase suffix, so `H2O_l_s` parsed as a second `H2O` and condensed
carriers were counted as gas, which made the cold trap invisible in `X_O`
whatever `--count-condensates` was set to. The gas phase is now selected by the
`.vul`'s own `atm['gas_indx']`, condensate is excluded from `mu` always as
VULCAN itself does, and the flag now selects. The Photochem adapter and
`lower_profile_schema.py` are untouched, so the production arm and its
`solution_id 5130...` are unaffected. Second, a configuration trap worth
stating because it silently produces a wrong atmosphere rather than an error:
`vulcan_cfg.He_H` is NOT read on the FastChem `ini_mix = 'EQ'` path
(`build_atm.ini_fc` builds the element file from `atom_list` minus H), so
helium must be in `atom_list` or the run keeps the solar `He 10.9864` and
integrates at `He/H = 0.097`. The first run did exactly that.

Whether VULCAN can be made to converge here is open. It is not a blocker for
Phase E: Photochem is the production path precisely because it carries the
climate model this planet needs, and T-E5 through T-E7 are closed without it.

## 80. A restart onto a reservoir the handoff has moved (2026-08-27)

Section 79 left T-E8 closing for helium and missing for hydrogen, carbon,
nitrogen and oxygen by 5.6e-8, 1.2e-4, 1.2e-4 and 3.3e-4, and named the cause:
`load_IC` restores every metal density the restart file carries, so a
flux-closure iteration that moves an elemental reservoir keeps the seed's value
everywhere above the base cell while the boundary condition already uses the
new one. The fix was deferred there because it changes restart semantics, which
is user-visible. The decision was taken on 2026-08-27; this section is what was
implemented and what it measures.

### The rule

An element is renormalized as it is loaded if, and only if, the
lower-atmosphere handoff states its reservoir -- a `<El>_H_base` key of
`base.inp`, or an elemental ratio of the file named by
`Lower atmosphere profile:`. For such an element the whole loaded column is
multiplied by one factor,

```
r_El = (stated El/H) / (El/H of the loaded state at the base cell),
```

the same factor for every ionization stage. The normalization moves and
nothing else does: the ionization split of each cell and the radial shape of
the loaded profile are preserved exactly, because a factor common to all
stages and all cells cancels out of every ratio the state is made of. `El/H`
at the base cell is counted in nuclei -- the element summed over its stages
against the hydrogen nuclei carried by H I, H II, H2, H2+, H3+ and HeH+ --
which is the convention of `element_ratio_HeH` and of
`src/utils/element_budget.py`, so the factor is measured in the same units the
budget check reads it back in.

Helium is not part of this rule and keeps its own, older convention: the base
cells are set to the input `He/H`, and with `He_diffusion` on the column above
keeps the diffused split, because that split is the state being restarted
rather than a defect of the file.

What decides "the handoff states it" is not a test on the numbers. Both
readers set an abundance through one routine, `set_element_abundance` in
`input_read.f90`, and passing through that door is what marks the element in
the new `melem_from_handoff(n_melem)` of `global_parameters`. An abundance
that reached the run from `metals.inp` alone never passes through it, so it is
not marked and its column is loaded untouched. There is no tolerance and no
`abs(r - 1) > tiny` guard: the branch is entered or it is not, and when it is
not entered the arithmetic is never performed.

Each rescaled element prints its own line, with the loaded ratio, the stated
ratio and the factor, so the factor can be compared against the reservoir
change it is supposed to equal:

```
 (load_IC) C/H at the base cell: restart file  2.777530E-04, handoff  2.777859E-04; column rescaled by  1.000118E+00
 (load_IC) O/H at the base cell: restart file  4.956907E-07, handoff  4.955292E-07; column rescaled by  9.996741E-01
 (load_IC) N/H at the base cell: restart file  8.187442E-05, handoff  8.188411E-05; column rescaled by  1.000118E+00
```

### What it closes

The test is the last iteration of the LHS 1140 b closure of section 79,
`LHS1140b/exhale/flux_closure/hi/k06`, re-run from the same seed
(`k05/output/`) on the same profile, once with the previous binary and once
with this one. The three factors above are the run's own log. They are the
measurement, not an illustration: `1.000118` and `0.9996741` are the same
numbers as the `1.183e-4` and `3.259e-4` the budget was missing by, which is
what the section 79 diagnosis said they would be.

| element | budget before | budget after |
|---|---|---|
| H (against `rho/mass_per_H`) | 5.633e-08 | 1.543e-15 |
| He | 1.910e-15 | 1.698e-15 |
| C | 1.183e-04 | 3.688e-14 |
| N | 1.183e-04 | 3.724e-14 |
| O | 3.260e-04 | 4.017e-14 |

`element_budget.py` goes from "4 element(s) miss their reservoir by more than
1.0e-08" to "all element budgets close within 1.0e-08". Hydrogen closes with
the rest and for the same reason: it is checked against `rho/mass_per_H`, and
the metal mass a stale reservoir leaves in `rho` is exactly what put it 5.6e-8
out.

The effect on the physical result is small, as section 79 predicted from the
size of the reservoir change. The steady mass flux at the top of the domain
moves by `-2.4e-5` in relative terms (`log10 Mdot = 7.50` either way), and the
He 10830 line, run through `EXHALE_transit.py` on both solutions, moves by
`4.7e-5` in relative terms in the red-component depth (1.725492 -> 1.725411)
and by at most `1.2e-6` in the normalized line profile.

### What is unchanged, checked rather than argued

Two restarts that must not move, both re-run with the two binaries and
compared bitwise on `output/Hydro_ioniz.txt` and `output/Ion_species.txt`:

- `LHS1140b/exhale/heh0p55` reloaded and advanced one step -- a restart with no
  handoff of any kind. Identical.
- `backup/regression/wasp_full` seeded from its own converged output with
  `Load IC? True` and advanced one step -- a restart that carries a full metal
  state, from `metals.inp` and not from a handoff, so it exercises the metals
  block of `load_IC` with every element unmarked. Identical.

The binary the comparison was made against is kept as
`backup/EXHALE_pre_handoff_restart_rescale.x`.

`make check` is unaffected by construction and by measurement: all seven
default cases are cold starts (`Load IC? False`), so none of them reaches this
code at all, and the full matrix -- `wasp_full`, `wasp_he23off`,
`mol_base_handoff`, `mol_metals`, `mol_lyman_werner`, `mol_diffusion`,
`lower_profile` -- comes back byte-identical (`==> REGRESSION PASS (all cases
byte-identical)`, 2026-08-27).

### `r_El = 1` is not an identity, and the reason is worth stating

The gate is the flag, not the value of the factor, and that is deliberate:
`r_El = 1` does **not** leave the state bitwise where it was. Restarting
`backup/regression/wasp_full` from its own converged output twice with the
same binary -- once with the abundances coming from `metals.inp` alone, once
with a `base.inp` restating those same seven numbers as `<El>_H_base` -- prints
`1.000000E+00` for every element and still gives a different state after one
step. The factors are 1 only to about `1e-14`: the seed's base-cell `El/H` is
what the ionization solver converged to, and it reproduces the reservoir to its
own tolerance rather than to the last bit. After one step the departure is
`1.8e-16` in the median and `1.5e-6` at its worst, and the worst cell is
`O III` at `n = 11 cm^-3` near the base, six orders below that ion's own peak
in the same column -- the solver's tolerance seen from the other side, not a
change in the solution.

So the byte-identity claim is exactly the one the flag supports and no wider:
**an element the handoff does not state, and a restart with no handoff at all,
never enter the arithmetic**, and those are what were compared bitwise against
the pre-change binary above. An element the handoff does state is renormalized
whenever it is loaded, which is the user-visible change of restart semantics
this section is about.

---

## 81. Is the wind collisional where it is set? A Knudsen diagnostic, and LHS 1140 b's answer (2026-08-27)

Until now the repository had no way to ask whether a converged solution is a
solution of the physical problem. A hydrodynamic wind is a continuum
solution, valid only where the gas is collisional on the scale over which
the flow varies -- above all through the critical point, which is where the
topology, and with it the mass flux, is fixed. Nothing measured that.
`src/utils/collisional_validity.py` (new, Python post-processing; it reads an
existing run directory and changes nothing, and no Fortran source was
touched) computes it. Definitions, provenance and the full result tables:
`docs/collisional_validity.md`.

**The collision model is not a new one.** The four momentum-transfer limits
are the ones the binary element-diffusion operator already uses, in the same
Chapman-Enskog first approximation: rigid-sphere neutral-neutral (Banks &
Kockarts 1973), non-resonant ion-neutral with the induced-dipole and
rigid-core frictions added, and screened Coulomb with the Spitzer
`ln Lambda` for ion-ion and ion-electron. Python cannot call the Fortran, so
the four routines are transcribed from
`src/modules/functions/binary_element_diffusion.f90` with the source line
cited at each definition and every constant taken from that file rather than
re-chosen; the two are meant to be read side by side. The step from a
diffusion coefficient to a collision frequency is the definition of the
binary coefficient itself,
`nu_st = n_t k T/(mu_st Dhat_st)` with `Dhat = n D`, in which the total
carrier density cancels.

Two channels are deliberately missing, and both make the verdict *more*
conservative by lengthening the mean free path: resonant charge exchange
(H+ + H, He+ + He), which the element-diffusion operator excludes on the
grounds that both partners carry the same element, and electron-neutral
momentum transfer, which no element-transport coefficient needs. The
electron Knudsen number is therefore quoted only with that caveat, and the
bulk is built from the heavy particles alone.

**Definitions, all stated.** `lambda_s = vbar_s/nu_s` with the harmonic sum
over partners; the structure scale `L = min(H_p, L_v, r)` with
`H_p = |dln p/dr|^-1` measured from the solution and `L_v = |v/(dv/dr)|`
admitted only above Mach 0.01 (below that it is the distance to a stagnation
point, and the base sound-wave layer would otherwise drive `Kn` to
infinity); `Kn_s = lambda_s/L`. The **exobase** is `Kn_bulk = 1`; the
**critical point** is the sonic point with the code's own sound speed
(`cs = sqrt(g p/rho)`, `g = 1.666666666667`, `eval_dt.f90:29`,
`parameters.f90:385`), not a re-definition; the **critical region** runs
from peak volumetric heating to the sonic point, and "collisional" means
`max Kn < 0.1` across it. Coupling times `1/nu_s` and the electron-ion
energy-coupling time are reported against the flow and heating times, with
the note that they are a *different* test: a strongly subsonic flow gives
hundreds of collisions per flow time even at `Kn ~ 1`, and it is the Knudsen
number that decides continuum validity, because the terms the closure drops
scale with the gradient over a mean free path.

**The result, measured 2026-08-27 on existing outputs.** Three
representative LHS 1140 b solutions -- `LHS1140b/exhale/heh0p55` (diffusion
off), `heh2p13_diff_kzz1e9` (diffusion at the adopted `K_zz`), and the
converged flux closure `flux_closure/hi/k06` -- against the HD 209458 b
control `backup/phase_d_baseline/new_kzz1e9_d3b`:

| case | critical point | exobase | `Kn = 0.1` at | max Kn, crit. region | verdict |
|---|---|---|---|---|---|
| `heh0p55` | none in domain (max Mach 0.540) | 25.22 | 6.42 | 2.2 | unvalidated |
| `heh2p13_diff_kzz1e9` | none in domain (max Mach 0.565) | 26.66 | 6.79 | 1.9 | unvalidated |
| `flux_closure/hi/k06` | none in domain (max Mach 0.550) | 19.76 | 5.53 | 3.0 | unvalidated |
| HD 209458 b `new_kzz1e9_d3b` | 4.085 | above 4.15 | -- | 0.021 | validated |

Radii in `R_p`; the LHS 1140 b domain ends at 30 `R_p`.

**LHS 1140 b's wind never reaches its critical point inside the domain.** It
is subsonic to 30 `R_p` in all three solutions, so the mass flux is set at
the outer boundary -- where `Kn_bulk` is already 1.3-2.3. Against the
p-winds retrieval, which places its isothermal Parker sonic point at 8-9.5
`R_p`: that radius lies *inside* the exobase computed here (19.8-26.7
`R_p`), but `Kn_bulk` there is already 0.14-0.30. So the answer is neither
"the sonic point is outside the exobase" nor "the critical region is
collisional" -- it is **transitional**. Neutral hydrogen sets the number;
the ions and electrons are Coulomb-held two to five orders of magnitude more
tightly. HD 209458 b is the control at a factor of ~100 lower `Kn`, and it
validates cleanly.

**The kinetic estimate is a scale, not a bound.** Jeans escape through the
computed exobase (Chamberlain & Hunten 1987 eq. 7.2.5) gives 1.0-2.0e7 g/s,
a factor about three below the continuum `Mdot` of 10^7.47-10^7.81. But
`lambda_J(H) = 0.78-1.45` at that exobase: it is barely gravitationally
bound, the Jeans integral is most of the Maxwellian rather than its tail,
and the atmosphere is in blow-off, so the number does not bound anything.
What it does say is that the continuum result is within a factor of a few of
what a collisionless outer atmosphere at the same density and temperature
would lose -- the continuum failure is not hiding an order of magnitude. On
HD 209458 b, `lambda_J = 10.8` and the hydrodynamic rate exceeds Jeans by
1e4, which is the signature of a driven wind rather than an evaporating one.

The numbers above are read from the advection-corrected profiles, the ones
the transit tool consumes. Read from the equilibrium pair instead (`--eq`)
the verdicts are unchanged -- still no critical point in the LHS 1140 b
domain, still `Kn ~ 1` where the outer boundary sets the flux, still
collisional through HD 209458 b's critical point -- but the exobase moves
above 30 `R_p`, because the equilibrium far field is more ionized and
Coulomb collisions there are far stronger than the neutral ones. An exobase
radius quoted from this tool has to say which pair it came from.

This closes the plan's baseline row 11 the way that row demanded: an
uncollisional critical region makes the hydrodynamic `Mdot` **unvalidated,
not overestimated**, and the tool's `validity_statement` says exactly that,
with the numbers it rests on. What remains unsettled is the LHS 1140 b rate
itself, which needs a kinetic or transitional-flow calculation this
repository does not contain. The He 10830 line-forming layer is not
implicated: it sits at 1.1-3 `R_p`, where `Kn <= 0.021` in every case.

## 82. The third body was a mass density: a molecular-network audit in the He-dominated limit (2026-08-27)

Phase F of the LHS 1140 b plan asks for an audit of the molecular network
where it is least tested -- a helium-dominated envelope. R16-R20 and R23,
the shared H-He charge-exchange path, the two He 2^3S Penning channels, and
the electron and H-nucleus closures were checked reaction by reaction
against the publisher PDFs, and the closures were measured on runs spanning
He/H = 0.079 to 1000. The record, with the full reaction table and the
numbers quoted below, is `docs/molecular_chemistry_audit_he_rich.md`.

**No rate coefficient was wrong.** All 23 Koskinen et al. (2022) Table-1
entries are transcribed verbatim and correctly, R16-R20 and R23 included;
so are the Huang et al. (2023) Table-4 B1/B2 charge-exchange pair, the
Taylor et al. (2025) He(2^3S)+H fit and the Garcia Munoz (2025)
He(2^3S)+H2 fit (which reproduces its four tabulated points to 0.13%, remeasured).
Every one of the eight balance rows of `mol_heh_rows` carries each reaction
in the rows it belongs to, with the right weight and sign.

**What was wrong was the density handed to the three-body reactions.**
R12 (H2 + M), R13 (H+ + H2 + M) and R15 (H + H + M) are written as
two-body coefficients times the third-body density, and the interface that
carries it, `ion_cell_state%ntot`, says "total particle density". What
`ionization_equilibrium` actually passed was `n_in_dim` = `rho/m_H`. That
is a *mass* density in units of m_H -- `calc_rho` weights every species by
`bsp_mass`, so an H2 molecule contributes 2 and a helium atom 4 -- and it
over-counts the third bodies by the mean particle mass. At the hot-Uranus
molecular base the factor is 2.3; in a helium-dominated base it approaches
4. The error therefore grows along exactly the axis this audit was asked
to examine. The same `rho/m_H` was also used as the gas pressure of the
chemical-equilibrium seed that restarts a failed molecular cell, where the
EOS's own law is `p = (n_tot + n_e) kB T`. Both now use `calc_ntot` (one
particle per species, electrons excluded) and, for the pressure, that sum
plus `n_e`.

The consequence is visible where it should be. On the `mol_diffusion`
configuration re-run at He/H = 1, 10 and 1000 (only `He/H number ratio` and
the `HeH_base` of the handoff changed, 12000 marching steps,
`OMP_NUM_THREADS=1`), the He/H = 10 arm no longer aborts on the marching
loop's NaN detector and completes all 12000 steps, and the number of cells
whose molecular root leaves the physical simplex and has to be clamped onto
the element budget falls from 37 to 4 at He/H = 1000 and from 181 to 169 at
He/H = 10. The He/H = 1 arm still aborts, later (5538 steps against 3402);
it is not diagnosed here.

**A guard that was not guarding.** The molecular systems close neutral
helium as `n_He(1 - x2 - x3) - n_HeH+`, but the admissibility test
`ionization_fractions_physical`, the distance measure
`element_budget_violation` and the clamp `clamp_fractions_to_element_budget`
all summed the helium budget as `x2 + x3` (plus the metastable) alone. A
root with `x2 + x3 = 1` and any HeH+ was accepted although its neutral
helium is negative. All three now carry the HeH+ nucleus, converted with
this cell's `n_H/n_He`. Nothing in the measured runs reaches that state --
HeH+/He stays below 1e-10 at every He/H tested -- so the change moves no
result; it closes a hole rather than fixing a symptom.

**The closures hold.** With molecules on, charge neutrality
`n_e = n_H+ + n_H2+ + n_H3+ + n_HeH+ + n_He+ + 2 n_He++` closes against the
`n_e` the code writes to 4e-16 at every He/H from 0.079 to 1000, the
H-nucleus/mass closure to 4e-16, and `src/utils/element_budget.py` closes
both elements to 6e-14 at He/H = 1000. Element conservation is not what
the helium-dominated limit breaks. What it does break is the solve: the
molecular `hybrd1` system balances helium in rows of order `n_He^2` and the
molecules in rows of order `n_H^2`, so at He/H = 1000 the two blocks of one
residual vector differ by 1e6 with no row scaling applied, and isolated
cells fall into the atomic basin next to neighbours at `f_H2 = 0.88`. The
correlation with He/H is measured; the cause is not established.

**Four caveats are recorded rather than patched**, because each would
replace a published rate with one from a different compilation. (i) The
identity of M: Koskinen et al. never say what their "n" counts, and a
monatomic third body is generally less efficient than H2, so R12/R13/R15
with helium as M are upper bounds. (ii) HeH+ formation: Table 1 has only
`He+ + H2` (R20), and not `H2+ + He -> HeH+ + H`, which needs no He+ and is
therefore the route that grows with the helium fraction (Garcia Munoz 2025
Table A.6 gives 1.0e-11 at 2000 K rising to 1.5e-10 at 1e4 K); the ~10%
associative branches of both Penning channels end in HeH+ too and are not
resolved. Against the same tables R16 is 3.4-8.6x low, R18 1.2x high, R19
1.4-2.6x low, R20 14x high and R17 up to 1.9e4 high. (iii) The H <-> He
charge-exchange pair does not satisfy detailed balance: `k(H+ + He)/k(He+ + H)`
should be `4 exp(-127500/T)` and the published fits give 1e2 times that at
1e4 K, a discrepancy weighted by `n_He/n_H` in the He+ balance. (iv) The
Taylor two-branch Penning fit steps by 1.57x at its own 4000 K break point,
which is the published fit and not a transcription error. Each is noted at
the code site as well.

**Does this move the LHS 1140 b results?** No, and the answer is checked
rather than assumed: **no LHS 1140 b run in the repository has molecular
chemistry on.** None of the 170 `input.inp` files under `LHS1140b/` sets
`Molecular chemistry`, and none of the 193 `Ion_species.txt` files carries an
H2/H2+/H3+/HeH+ column. Every statement corrected above lives inside a
`thereis_mol` branch, so the He/H = 2.09 solution, the closure ladder and the
`heh2p13_diff_kzz1e9` arm are untouched and no re-convergence is needed. What
*can* be measured on those stored solutions is the size of the factor that
would have been wrong had they been molecular -- the mis-supplied `rho/m_H`
against the EOS particle count `p/(kB T) = n_tot + n_e`, at the base cell:

| stored solution | He/H | (rho/m_H)/(n_tot+n_e) at the base |
|---|---|---|
| `LHS1140b/exhale/flux_closure/hi/k06` | 2.09 | 3.03 |
| `LHS1140b/exhale/heh2p13_diff_kzz1e9` | 2.13 | 3.04 |
| `LHS1140b/exhale/flux_closure/heh9p7/k03` | 9.7 | 3.72 |
| `backup/regression/golden/mol_diffusion` | 0.0793 | 2.27 |

A molecular LHS 1140 b base would therefore have run its three-body reactions
with a third body 3.0-3.7 times too dense, against 2.3 times at the
hot-Uranus base the regression case covers.

**Goldens.** The matrix splits exactly along the `thereis_mol` branch:
`wasp_full`, `wasp_he23off` and `lower_profile` PASS byte-identical, and the
four molecular cases move. Their shift on the 12000-step snapshots, largest
relative difference over the column against the previous golden:

| case | rho | p | T | heat | cool | H2 front (f = 0.5) | log10 Mdot |
|---|---|---|---|---|---|---|---|
| `mol_base_handoff` | 1.4% | 2.6% | 1.2% | 1.6% | 3.0% | 1.1597 -> 1.1597 | 10.58 -> 10.58 |
| `mol_metals` | 3.0% | 5.0% | 2.2% | 3.5% | 6.7% | 1.1617 -> 1.1617 | 10.58 -> 10.58 |
| `mol_lyman_werner` | 0.34% | 0.53% | 0.21% | 0.44% | 0.71% | 1.1304 -> 1.1304 | 10.58 -> 10.58 |
| `mol_diffusion` | 1.4% | 2.6% | 1.2% | 1.6% | 3.0% | 1.1597 -> 1.1597 | 10.58 -> 10.58 |

The velocity column moves by more (2.5-8x) but only in the base sound-wave
layer of these unconverged relaxation snapshots, where `v` is small and
oscillating. The H2 front does not move at all, and the base H2 fraction
shifts in the fourth decimal: these are H-rich hot-Uranus bases already at
`f_H2` ~ 1, and a fraction pinned at unity cannot respond to a 2.3x change in
the third-body density -- it responds in rho, p and T instead. That is the
same statement the He-rich scan makes from the other side: the term becomes
decisive only once helium dilutes the H2. Goldens re-snapshotted for the four
moved cases at the end of the series, and `make check` re-run after.

## 83. The outer region a converged wind inherits: path dependence beyond the exobase (2026-08-27)

The LHS 1140 b closure ladder (`LHS1140b/kzz_decision.md` sec. 8) climbs a
reservoir axis by seeding each arm from the converged wind of the arm below.
One rung was reached by a large jump instead of a step, and it converged to a
different solution -- different only where the continuum equations no longer
hold. This section records what separates the two, why neither convergence
test noticed, and what it means for reading any of these winds.

**The two solutions.** Both are flux-closure arms at the same reservoir
`He/H = 12.0` on the same photochemical lower atmosphere. One
(`LHS1140b/exhale/flux_closure/heh12/`) was seeded from the 8.0 arm, a 1.50x
jump in reservoir, and closed at k = 4. The other
(`.../heh12_ctl/`) was seeded from the 11.1 arm, a 1.08x step, and closed at
k = 2.

| | seeded across a 1.50x jump | seeded one 1.08x step |
|---|---|---|
| `He/H` at the 1 microbar match | 12.013522 | 12.013520 |
| cold-trap `O/H` of the handed-over column | 1.0656e-6 | 1.0660e-6 |
| `exhale_info` | 0 | 0 |
| achieved wind `\|\|R\|\|` | 2.529e-4 | 3.864e-4 |
| steady-window flux spread | 1.74 % | 3.30 % |
| `log10 Mdot` [g/s] | 7.420 | 7.410 |
| He 10830 red depth [%] | 4.630 | 3.856 |
| He 10830 red-pair EW [%A] | 1.3126 | 1.1181 |

Identical base, identical matching-level composition to seven digits, both
converged and both inside every tolerance the driver applies -- and a 17 per
cent difference in the observable the whole study is fitted to. Neither
residual separates them: on the coefficient in force when this was written the
outlier reported the *smaller* of the two, and re-solved on the Garcia Munoz
rate of section 87 the order reverses (jump-seeded `2.533e-4` against the
control's `1.249e-4`, both `info = 0`, both inside the `4.0e-4` target). Which
one reports less is not a property of the pair; that `||R||` does not rank them
is (`LHS1140b/exhale/xuvkzz_gm25/results.txt`).

**Where they differ.** Nowhere the line window or the closure window looks.
Two measures that do not depend on any sampling choice locate it. `r_drop` is
the first radius outside 1.5 `R_p` at which `T` falls below half of
`T(12 R_p)`; the other is the share of the He 2^3S radial column that lies
outside 10 `R_p`. Both are monotonic along the ladder, and both single out
one point:

| reservoir `He/H` | `r_drop` [`R_p`] | `T(12 R_p)` [K] | He 2^3S column outside 10 `R_p` |
|---|---|---|---|
| 2.09 | 30.00 | 893 | 0.0119 |
| 3.0 | 27.63 | 1035 | 0.0116 |
| 5.0 | 23.76 | 1310 | 0.0105 |
| 8.0 | 20.80 | 1632 | 0.0104 |
| 9.7 | 19.79 | 1789 | 0.0108 |
| 10.3 | 19.79 | 1834 | 0.0108 |
| 10.7 | 19.47 | 1863 | 0.0109 |
| 11.1 | 19.47 | 1891 | 0.0109 |
| 12.0, one step | 19.47 | 1951 | 0.0112 |
| **12.0, one jump** | **25.83** | 1957 | **0.0219** |

The wind gets hotter and more compact as the reservoir rises: `r_drop` moves
inward from 30.0 to 19.5 `R_p` without a reversal, and the outer share of the
metastable column stays in 0.0104-0.0119 over a factor 5.7 in composition.
The jump-seeded solution sits at 25.8 `R_p` where its own reservoir gives
19.5, and carries twice the outer column share of any other rung, while its
`T(12 R_p)` matches its control to 0.3 per cent. It is the same wind with a
hotter, more extended outer atmosphere bolted on.

**The mechanism, read from that pattern.** The outer state is inherited from
the seed and the JFNK finish does not re-solve it. The density out there
contributes almost nothing to the residual norm, so `||R||` cannot separate
the two states; the elemental-flux spread the closure tests is measured on a
steady window that is flat in both. He 10830 is optically thick at these
columns, so a hot, extended outer region fills the line wings and lifts the
equivalent width even while holding about one per cent of the metastable
column. Small mass, large opacity, invisible to both convergence tests.

**The physical judgement, stated first.** The 20-30 `R_p` region in which the
two solutions differ is where the collisional-validity diagnostic of section
81 places the exobase, 19.8-26.7 `R_p`, and where `Kn > 1`. *The fluid
solution is not valid in exactly the region that makes the difference.* Two
consequences follow, and neither is comfortable:

- The ground for discarding the jump-seeded solution is **continuity of the
  trend, not a physical criterion**. Nothing in the equations being solved
  prefers one outer state over the other, because the equations do not apply
  there.
- By the same argument, the absolute equivalent widths of **every** rung of
  that ladder carry a systematic error the closure tolerance does not catch.
  The size is small -- the region outside 10 `R_p` is about one per cent of
  the metastable column -- but it is not zero, and it is not bounded by any
  test currently applied.

**What this says about the convergence criteria.** `exhale_info = 0`, a wind
residual under target, and a flat elemental flux across the steady window are
together sufficient to certify the solution *inside* the collisional region
and are silent outside it. That is not a defect of the criteria: a residual
norm weighted by density cannot see a region that carries no density, and a
flux-spread test measured over a window is blind to what the window excludes.
It is a limit on what "converged" means for a subsonic wind whose domain
extends past its own exobase. Where an observable is optically thick and
integrates over the whole column -- He 10830 here, and Ly-alpha more so -- the
certified region and the region that sets the observable are not the same
region.

**What this says about running a scan.** A ladder scan must be climbed one
step at a time. The step size that worked here was 1.08x in reservoir; the
1.50x jump did not. There is no way to detect the failure from the solver's
own output, so the protection has to be procedural: keep the step small, and
check a sampling-free outer diagnostic (`r_drop`, or the outer share of
whatever column the observable integrates) for monotonicity along the scan.

**What was and was not tested.** Path independence was checked directly at
`He/H` = 3.0, from a seed 3.7x away in composition and in the downward
direction the ladder never takes: it returns the same `r_drop` to the digit,
with everything else inside the closure tolerance. The 5.0, 8.0 and 2.09 arms
were not checked this way. Determinism was checked separately at 10.31
(`.../heh10p3_ctl/`, same seed, re-run): `Hydro_ioniz_adv.txt` and
`Ion_species_adv.txt` reproduce bitwise, which says nothing about seed
dependence. **No code was changed by this section**; the ladder table, the
crossing and the figure were re-derived on the control arm, and the record
is `LHS1140b/kzz_decision.md` sec. 8.3 and
`docs/lhs1140b_exhale_vs_pwinds.tex` sec. "The outer region a converged wind
inherits".

## 84. Metals had no way back: element projection when the helium mass fraction reaches one (2026-08-28)

**The judgement first.** The element projection of the binary H/He diffusion
operator destroyed trace-metal nuclei with no sink and no way back, and that
is physically wrong whatever it costs. In a cell whose helium mass fraction
reaches `X = 1` the hydrogen count `n_H^new = (1 - X) rho/m_1` is exactly zero
and the projection multiplied the metals by `r_H = n_H^new/n_H^old = 0`.
Zeroing them *at that instant* is right: in this two-component split the metals
are part of component 1, slaved to hydrogen at a fixed metal/H, so a cell with
no hydrogen carries no metals and component 1 carries no mass. What is wrong is
the step after. Hydrogen returns to that cell through the shortfall deposit at
the end of the same projection; **the metals do not**, and nothing later can
lift them, because every writer to the metal columns of `f_sp` on this path
multiplies what is already there. The result is a contiguous band of cells
holding exactly zero carbon, nitrogen and oxygen between cells at the reservoir
ratio -- a state these equations cannot produce, since they contain no sink for
a metal nucleus. The fix is that the metals return **with** the hydrogen.

### The path, in the code

`element_diffusion_step` solves for the helium mass fraction `X` and hands the
new element totals to `project_elements`
(`src/modules/functions/binary_element_diffusion.f90`). There, per cell,

```fortran
nucH_new  = (1.0d0 - Xhe(j))*msum(j)/m_1
if (nucH_old(j) .gt. 1.0d-30) then
   rH = nucH_new/nucH_old(j)
else
   rH = 0.0d0
endif
...
do im = 1, n_mion                       ! metals slaved to hydrogen
   f_sp(j,mion_fsp(im)) = f_sp(j,mion_fsp(im))*rH
enddo
...
if (nucH_new  .gt. gotH)                                         &
   f_sp(j,isp_HI)  = f_sp(j,isp_HI)  + (nucH_new  - gotH)
```

Step `k`: `X` reaches 1, so `nucH_new = 0`, `rH = 0`, and every hydrogen
species and every metal ion of the cell is multiplied by zero. Step `k+1`:
`nucH_old` is now zero, so `rH = 0` again for want of a ratio to preserve,
while `nucH_new > 0` -- and the deposit on the last line puts the hydrogen
back as neutral `H`. There is no such line for the metals. Their zero is
final.

### The symptom that led here

`LHS1140b/exhale/branch_diagnosis/` records the measurement (no code was
changed by that work). At `K_zz = 1e7` the reservoir scan reported two
converged states at the same reservoir `He/H`, `log10 Mdot` 7.54 against 7.87
and He 10830 red-pair equivalent width 0.4287 against 1.3248 %A, both with
`info = 0`. The hot one carries 172 cells, 1.0333-1.6835 `R_p`, at exactly zero
C, N and O in every stage, with the same elements at `C/H = 2.639e-4`
immediately above and below. On the cool solution the C I lines carry 95.6 % of
the volume-integrated cooling over 1.0-1.2 `R_p`; deleting the carbon drops the
integrated cooling of that band by a factor 23 at 0.91x the heating and lifts
`T max` from 3899 K to 5879 K. Refilling C, N and O at the handoff's own
reservoir ratio and re-solving the wind with the same handoff collapsed the hot
state onto the cool branch (`log10 Mdot` 7.530, `T max` 3890 K, EW 0.4065 %A).
Of the 550 `output*/Ion_species.txt` files in the tree, 130 carry such a band;
128 are in `LHS1140b/exhale/kzz_profile_scan`, the other two are runs already on
record as failed. **No regression golden and no `backup/regression/` case output
carries one**, which is why the repair below moves nothing.

### Why `info = 0` and `||R|| ~ 1e-4` could not see it -- read beside section 83

Both states report `info = 0` and a residual of a few times `1e-4`, and both
are right to. **The steady residual is a residual of the hydro and energy
equations**, and their solution with the metals removed is a perfectly good
solution of the equations as posed: the cooling function simply has fewer terms.
Nothing in that norm asks whether the composition it was evaluated on is a
composition the transport equations can reach. The elemental-flux closure does
test the composition, but on the radial spread of the elemental *fluxes* over a
window at 2 `R_p` and above -- outside the band. The handoff's own elemental
conservation check (`src/utils/lower_profile_schema.py`) is stated on the
*chemistry* profile at its deepest level, not on the wind. **Elemental
conservation was nowhere checked for the wind**, which is the gap this section
closes.

Section 83 is the same shape of failure at a different place: there two
converged solutions differed only outside the exobase, where the residual norm
carries no density and the flux window does not reach, and He 10830 saw the
difference because it is optically thick and integrates the whole column. The
statement the two sections share is that `info = 0` plus a residual under target
plus a flat elemental flux certify the solution **on the quantities and in the
region those tests are evaluated on**, and are silent elsewhere. What is new
here is that the silence is not about a region: it is about a *conserved
quantity that was never formed into a test at all*.

### This is a third failure mode, not either of the first two

- **Section 70** (HD 209458 b, `K_zz = 0`): two wind solutions a factor 2.8
  apart in `Mdot`, both `info = 0`, alternating in a period-2 limit cycle of the
  outer composition/wind Picard loop, on compositions ~1 % apart. Every
  `backup/phase_d_baseline/` run was checked for the band and **none has one**,
  so that bistability is not this defect and stays open. (Caveat: only the final
  state of each run is on disk; the alternating passes of
  `new_kzz0_undamped` were not saved with species densities.)
- **Section 83** (LHS 1140 b, jump-seeded rung): two solutions that agree below
  the exobase to seven digits and differ only above it, where the fluid
  equations do not apply. Here the two states differ from `r = 1.0004 R_p`
  upward -- by 30 % in `T` by 1.0013 `R_p` -- and the difference is made in the
  collisional part of the atmosphere (`Kn_bulk ~ 1e-3` at `T max`).
- **This section**: a state that is not a solution of the equations at all,
  because an element has been deleted from a band of cells by the discretization
  and nothing can put it back.

The three should not be merged.

### The repair, and the abundance it uses

The metals are a *ratio*, not a density: what the projection has to preserve for
them is `n_X/n_H`. Multiplying by `r_H` does that in a cell that had hydrogen,
and that branch is untouched. Where `n_H^old` is zero the ratio is `0/0` -- the
cell carries no metal/H of its own -- and the operator supplies one:

```fortran
if (nucH_old(j) .gt. 1.0d-30) then
   do im = 1, n_mion
      f_sp(j,mion_fsp(im)) = f_sp(j,mion_fsp(im))*rH
   enddo
else if (thereis_metals) then
   do ie = 1, n_melem
      i0e = melem_i0(ie)
      f_sp(j,mion_fsp(i0e)) = melem_ab(ie)*nucH_new
      do k = 1, melem_top(ie)
         f_sp(j,mion_fsp(i0e+k)) = 0.0d0
      enddo
   enddo
endif
```

`melem_ab` is not one candidate among several. It is the metal/H that **this
operator's own mass budget already assumes**: the mass per hydrogen nucleus of
component 1 is
`m_1 = mass_per_H_nucleus_without_He() = m_H + sum_X melem_ab(X) A_X`, so the
two-component closure `m_1 n_H + m_He n_He = rho` that the write-back is
required to meet is stated at that ratio and no other. It is also the reservoir
ratio `set_IC` builds a cold start from, the one `load_IC` rebuilds an element
the restart file does not carry from, and the one section 80 rescales a restart
column onto. The alternatives were considered and rejected: the cell's own
previous ratio is exactly what does not exist there, and a neighbour's ratio
would make the result depend on the grid.

All of it goes into the **neutral ground stage**, for the same reason the
hydrogen and helium shortfalls go into `HI` and `HeI`: a cell that held no
hydrogen held no ionization split either, `ioniz_eq` owns the split within an
element and re-solves it from the element total on the next call, and this is
already the convention of `set_IC`, of `load_IC`, and of the trace-metal
re-seed branch of `element_diffusion_step` itself.

The repair prevents the band from being *created*. It does not reconstruct one
that a restart file already carries: there `n_H^old > 0` in the emptied cells,
so the preserving branch is the one that runs, and it faithfully preserves a
ratio of zero. That is the right behaviour for a restart -- silently rewriting
loaded data would be worse -- and it is why the census below is unconditional.

### The check the pipeline did not have

`metal_hydrogen_ratio_departure` now runs at the end of every diffusion step.
It measures `max |(n_X/n_H)/melem_ab - 1|` over the cells `1..N` and the
elements the reservoir states, and it counts the `(cell, element)` pairs holding
no metal nuclei at all while hydrogen is present. The departure is a magnitude
and not an error -- `he_metal_diffusion` moves the metals relative to hydrogen
on purpose, and a restart may carry a settled column -- so it is only reported,
under the existing `EXHALE_DIFFUSION_CHECK=1` switch, together with the rest of
the step diagnostic. The count is different in kind: it is zero on any state
these equations can reach, so it raises a warning on stdout the first time it is
nonzero **whether or not the diagnostic is on**. No new input key was added.

The step diagnostic also now reports how far the solve left `[0, 1]` *before*
the clamp, which is the subject of the next part.

### The clamp is a limiter, and the comment that said otherwise was wrong

The projection was preceded by

```fortran
! Round-off assertion (NOT a limiter): the M-matrix argument in the
! header gives 0 <= X <= 1 per step; only round-off can leave the range.
where (Xhe .lt. 0.0d0) Xhe = 0.0d0
where (Xhe .gt. 1.0d0) Xhe = 1.0d0
```

That claim had never been measured, because the diagnostic printed the
*clamped* `X`, which cannot tell an assertion from a limiter -- both print 1.0.
Measured now, on the reproduction run of the next part: `X` exceeds 1 before the
clamp in **140 of 7962 steps, by up to 0.520** (i.e. `X = 1.52`), and 105 of
those excursions exceed 0.1. The lower end holds: `min X` stayed above 0 in
every step of both runs. **The clamp is a limiter on the upper side.**

The reason is in the discretization, not in round-off. The settling flux
`-rho D_12 G X (1-X)` is linearized with the `(1-X)` factor lagged,
`W_f = rho_f D_f G_f (1 - X_f^lag)`, and entered as `-W_f X`. In all three
branches of the Peclet hybrid the face coefficients satisfy
`PdL(f) + PdR(f) = -W_f`, so the row sum of the backward-Euler matrix is

```
rho_j/dt  -  K_j ( r_{j+1/2}^2 W_j  -  r_{j-1/2}^2 W_{j-1} )
```

and **not** `rho_j/dt`. A uniform `X = 1` is therefore not a fixed point of the
discrete row unless the lagged settling coefficient is divergence-free, and the
comparison principle that the header's M-matrix argument invokes does not apply
-- at either end. What separates the two ends here is measurement, not proof:
the lower one held in every step of both runs, the upper one did not. The header
comment and the corresponding sentence of
`docs/binary_diffusion_design.md` section 3 ("a clip at round-off (1e-15) is
kept only as an assertion, not as a limiter") are corrected to state the
measurement.

**What is not done here.** Making the upper bound hold by construction means
changing the linearization of the settling term -- for instance solving the
`X(1-X)` product with both factors at the new level, or symmetrizing the Picard
sweep -- and that changes every diffusion-on result, including two regression
goldens. It is a formulation decision, it is recorded here, and it is left open.
What this section does is remove the *consequence* that made the clamp
destructive: with the metals returning with the hydrogen, a cell that touches
`X = 1` recovers instead of being permanently emptied. The `X = 1` excursions
themselves still happen -- 140 of them in the repaired run, against 142 in the
unrepaired one.

### What was measured

The step `kzz1e7/heh3p8` (clean) `-> kzz1e7/heh4p2`, re-taken with the same
binary in both columns, the projection branch being the only difference. This is
the step `LHS1140b/exhale/branch_diagnosis/metal_dropout_reproduce.py` takes,
run with `EXHALE_DIFFUSION_CHECK=1` and with the He I 10830 synthesis added;
both runs were made outside the repository, so no recorded output in it was
overwritten.

| | before | after |
|---|---|---|
| zero C/N/O cells, `Ion_species.txt` | 157 (1.0429-1.6716 `R_p`) | **0** |
| zero C/N/O cells, `Ion_species_adv.txt` | 157 | **0** |
| `exhale_info` | 0 | 0 |
| `log10 Mdot` | 7.52 | 7.52 |
| `T max` | 4182.6 K at 1.358 `R_p` | 4185.6 K at 1.358 `R_p` |
| He 10830 red-pair EW | 0.5306 %A | 0.5249 %A |
| diffusion steps | 7916 | 7962 |
| steps whose `X` reached 1 | 142 | 140 |
| max `X` before the clamp | 1.521 | 1.520 |
| steps with a vanished (cell, element) | 7541, up to 471 pairs | **0** |
| max metal/H departure from `melem_ab` | 1.0 (i.e. an element gone) | 8.1e-13 |

The unrepaired column reproduces the recorded diagnosis exactly -- 157 cells
over the same 1.0429-1.6716 `R_p`, 7916 diffusion steps, 142 of them touching
`X = 1` -- so the added diagnostics are numerically inert and the projection
branch is the sole cause of the difference. The repaired solution matches the
cool branch continued to this reservoir by refilling and re-solving
(`branch_diagnosis/cool_branch_continuation.txt`: `log10 Mdot` 7.520,
`T max` 4184.4 K, EW 0.5287 %A).

Note what the `T max` and EW columns do *not* say: in this particular step the
band happened to cost only 3 K and 1 % of the equivalent width. The factor-23
cooling loss and the 3899 -> 5879 K jump quoted above belong to the recorded
scan arms, where the band sits where the C I cooling is decisive. The band is a
defect wherever it appears; how much it moves the answer depends on where it
lands.

### Verification

- `make check`: 7/7 byte-identical. Nothing moves, and the reason is structural
  rather than lucky: the preserving branch `n_H^old > 0` is bitwise the previous
  code, and no golden state ever empties a cell of hydrogen.
- `make diffusion_tests`: 31/31 (28 as before, plus 3 assertions of the new
  T13). **T13 fails on the pre-fix projection** -- 21 of 21 band cells hold
  hydrogen and no metals, metal/H departure 1.0 -- which is the point of adding
  it. It seeds a closed metal-bearing column with a band of pure helium, so
  `n_H = 0` and `X = 1` exactly, lets diffusion refill the band, and requires
  that no cell holding hydrogen holds zero metals and that `n_X/n_H` equals
  `melem_ab` to 1e-12 (measured 5.3e-15).
- `run_fcheck.sh` (`-fcheck=bounds,do,mem`) clean. It had to be repaired first,
  for two reasons that predate this change and are unrelated to it: its repo
  root was still one level up from `backup/regression/` (so every `make` in it
  failed with "No rule to make target 'distclean'"), and it copies
  `HD209458b/input.inp` into an empty scratch directory without turning off the
  `Load IC? True` that file now carries, so `load_IC` stopped on a missing
  restart and the script reported that `error stop` backtrace as a runtime check
  firing. With both fixed the HD 209458 b bounded run is clean
  (`log10 Mdot = 9.60`).
- That case runs with element diffusion **off**, so it does not reach the code
  this section changes. Two runs were added for that: `diffusion_tests.x` built
  with the same flags (31/31, and T13 is the test that enters the new branch),
  and the `lower_profile` regression case -- element diffusion together with
  metals -- bounded to 400 steps. Both clean, no trap.

**Not done, and left to a decision.** The 128 contaminated arms of
`LHS1140b/exhale/kzz_profile_scan` are unchanged on disk. The repair prevents
the band from forming; it does not repair a stored state, and re-running that
scan is a data decision, not a code one.

**The open item above is closed by section 86** (2026-08-28): the two factors
of the drift product `X(1-X)` are now taken from opposite sides of the face
and both at the new time level, so the flux shuts off at either end of the
composition axis as it does in the continuum, and the `X > 1` excursions are
gone. The clamp is an assertion at both ends.

---

## 85. A structure scale that is a scale of the state, not of the pressure (2026-08-28)

`src/utils/collisional_validity.py` divides a mean free path by a local
structure scale `L` to get the Knudsen number that decides whether a wind
solution is a continuum solution. Until now `L = min(H_p, L_v, r)`, and
below Mach 0.01 the velocity term switched itself off, so in a subsonic
region -- LHS 1140 b runs at Mach 1.6e-4 near the base -- `L` was the
pressure scale height alone.

**The pressure is the one field that can go flat while the gas does not.**
`p = n k T`; across a heating front `rho` falls and `T` rises with nearly the
same logarithmic slope, so `d ln p/dr` passes through a broad near-zero
while the state changes fast. Measured in `LHS1140b/exhale/heh0p55`: `H_p`
runs 9.13e6 cm at 1.0244 `R_p`, 2.55e8 cm at 1.0453 `R_p` (its maximum,
3.05e8, is at 1.0486) and 1.41e8 cm at 1.0648 `R_p`, while at 1.0453 `R_p`
`H_T` is 1.80e7 and `H_rho` 1.68e7 cm -- fifteen times shorter -- and
`Kn(H I)` showed a hole of 5.4x (3.67e-5, 6.82e-6, 3.56e-5 at
those three radii) in a quantity that is monotonic on either side.

**The fix is the state vector.** The continuum closure expands about a local
Maxwellian, which is fixed by `(n, T, u)`; those logarithmic gradients plus
the spherical divergence are the complete set, and the change of the state
over a distance is a vector whose length is the sum in quadrature:

```
1/L^2 = (dln rho/dr)^2 + (dln T/dr)^2 + (|dv/dr|/c_s)^2 + (1/r)^2
L     = max(1/sqrt(that), dr_cell)
```

Pressure is not in it: it is a derived field, and adding it would count the
density and temperature gradients twice wherever they do not cancel while
contributing nothing exactly where they do. The quadrature sum, not a
minimum over the same lengths, because a minimum jumps discontinuously to
the next length whenever the field it is taken from flattens -- the failure
above in general form -- whereas in quadrature a flat field contributes
zero and `L` still cannot exceed the shortest individual scale. The
velocity enters as `|dv/dr|/c_s` and not `|(1/v) dv/dr|` because the term
the closure drops is the viscous stress, whose size relative to the
pressure is `~ lambda |dv/dr| / vbar`: a bulk-velocity change distorts the
distribution function in proportion to the thermal speed, not to the local
bulk speed. That removes the Mach-number floor -- `MACH_FLOOR` is gone --
which existed only because normalizing by `v` diverges wherever the base
sound-wave layer crosses `v = 0`. The floor at the local cell width says a
gradient claiming structure below one cell is measuring the mesh; it binds
in 3 of 500 cells of one of the four runs, all at `r < 1.002 R_p`, and in
no critical region.

**No conclusion moves.** Same runs, same outputs, nothing rebuilt or
re-run. HD 209458 b stays validated through its critical point (sonic point
4.085 `R_p`, which `L` does not enter; max Kn in the critical region 0.0214
-> 0.0186), all three LHS 1140 b solutions stay unvalidated with no critical
point in the domain (max Kn 1.9-3.0 -> 2.0-3.1). What moves is at the
10-30 % level and in both directions: the exobase 19.8-26.7 -> 18.8-25.4
`R_p`, `Kn_bulk` at the outer boundary 1.3-2.3 -> 1.4-2.4, `Kn_bulk` at the
8-9.5 `R_p` where the p-winds retrieval puts its sonic point 0.14-0.30 ->
0.13-0.29, and the He 10830 line-forming layer at 1.1-3 `R_p` `Kn <= 0.021`
-> `<= 0.019`. The one large change is inside 1.1 `R_p`, where the old
`H_p` sat on its isobaric plateau: `Kn_bulk` at 1.1 `R_p` goes from
4e-5 - 1.2e-4 to 2.5e-4 - 4.2e-4. Mach maxima, `log10 Mdot` and the
coupling times are independent of `L` and unchanged. `docs/collisional_validity.md`
carries the recomputed tables and a row-by-row account of what moved; the
numbers quoted in section 81 are the ones from the old scale and are
superseded by it. Figure `docs/figures/lhs1140b_knudsen.pdf` regenerated
(that block of `LHS1140b/make_memo_figures.py` only). Post-processing
diagnostic, no Fortran touched, goldens untouched.

## 86. The drift flux has to vanish where the composition runs out: the discretization of X(1-X) and the bounds on the helium mass fraction (2026-08-28)

**The judgement first.** The drift (settling) term of the binary H/He
transport equation carries the factor `X(1-X)` and vanishes at both ends of
the composition axis, because a cell with no helium has none to send and a
cell with no hydrogen has nothing to send it in exchange. The discretization
lagged the `(1-X)` factor at the previous iterate, which throws that shutoff
away: the discrete flux into a cell already at `X = 1` was not zero, and the
step overfilled it. Measured, `X` reached 1.52. That is not a small
inaccuracy, it is a discrete operator that does not have the property the
continuum flux has, and the clamp that hid it was doing the work of a
limiter. The drift product is now taken with its two factors from **opposite
sides of the face** -- helium from the cell it leaves, hydrogen from the cell
it enters -- and both at the new time level, and the resulting nonlinear step
is solved by Newton. The excursion outside `[0,1]` is gone: over every case
measured the solve now stays strictly inside the range.

Section 84 fixed the **consequence** (a cell driven to `X = 1` lost its metals
permanently) and left the cause open in as many words: "Making the upper bound
hold by construction means changing the linearization ... it is left open."
This section closes it.

### Where the old discretization lost the bound

The flux is (module header of
`src/modules/functions/binary_element_diffusion.f90`)

```
J = -rho (D_12 + K_zz) dX/dr  -  rho D_12 G X (1-X) ,
```

and the old face form entered the drift as one lagged coefficient times one
implicit factor,

```
W_f = rho_f D_f G_f (1 - X_f^lag) ,      J_drift(f) = -W_f X(f) ,
```

with `X(f)` central or upwind by the Peclet switch. Two things follow. The
row sum of the backward-Euler matrix is `rho/dt - K (r_R^2 W_j - r_L^2
W_{j-1})` and not `rho/dt`, so no comparison principle applies (section 84
established this much). More to the point, and this is the physical
statement: **a face between two cells both at `X = 1` still carries a drift
flux**, `-W_f`, unless the lag already knows that `X_f = 1`. The endpoint
states of the composition axis are not stationary states of the discrete
operator, while they are exact stationary states of the continuum flux. A
cell driven toward pure helium therefore kept receiving helium after it had
run out of hydrogen to exchange, and the clamp turned the overshoot into an
exact 1 -- which is where section 84's metal band came from.

### The discretization that keeps it

The drift is a **counter-flow**: the helium mass flux `-B X(1-X)` (with
`B = rho D_12 G`) is matched by an equal and opposite hydrogen flux, since
the two components close, `J_1 = -J_He`. A donor-cell rule for a counter-flow
has to take each element's mass fraction from the cell **that element**
leaves, and the two elements leave from opposite sides. With `B >= 0` helium
drifts inward, so helium is donated by cell `j+1` and hydrogen by cell `j`:

```
J_drift(f) = -B_f X(j+1) (1 - X(j))        B_f >= 0   (helium inward)
J_drift(f) = -B_f X(j)   (1 - X(j+1))      B_f <  0   (helium outward)
```

Both factors are at the new time level. A donor at `X = 0` sends nothing; an
**acceptor** at `X = 1` receives nothing. Those are the two statements that
bound `X`, and they are the discrete image of the two limits the continuum
factor `X(1-X)` expresses.

The Peclet hybrid stays, with its switch now on the drift coefficient itself,
`|B_f| dr <= 2 (A_f + E_f)`, and the central branch evaluates the product at
the face average, `X_f(1-X_f)`. That branch is inside the bound as well, and
for a reason worth writing down: with `X(j) = 1`,

```
X_f (1 - X_f) = (1 + X(j+1)) (1 - X(j+1))/4  <=  (1 - X(j+1))/2 ,
```

so the Peclet condition `|B_f| <= 2 A_f/dr` bounds the drift flux by the
gradient flux `A_f (1 - X(j+1))/dr` that is carrying helium *out* of that
cell at the same face. The central branch cannot reverse the sign of the
total face flux at a cell sitting on the bound. (The switch is stricter than
the old one by the factor `1-X`, so slightly more faces are upwind.)

**The bounds.** Let `X` solve the implicit step and suppose it first touches
1 in cell `m`, every other cell still inside `[0,1]`. Then the time term
`rho (X_m - X_m^old)/dt >= 0`; the gradient flux leaves `m` at both faces
because `X_m` is the maximum; the upwind advection contributes
`rho |v| (X_m - X_donor)/dr >= 0`; and the drift flux is an outflow or
exactly zero at both faces by the rule above. Every term of the row has the
same sign, and their sum is the row residual, which is zero -- so no term can
be strictly positive and `X_m > 1` is impossible. `X >= 0` is the same
statement, because the scheme is **exactly symmetric** under `X -> 1 - X`,
`B -> -B`: the donor and the acceptor exchange roles and the central product
is even about `X_f = 1/2`. That symmetry is the binary symmetry of the
transport equation itself -- there is one independent diffusive flux and the
two components carry it in opposite directions -- and the old lagged form
broke it, which is why one end held in practice and the other did not.

The row sum is still not `rho/dt` and does not need to be: the argument is
made cell by cell on the residual of the nonlinear step, not from a
comparison principle for a linear system.

### Newton, and why not another Picard sweep

The step is now nonlinear in `X^new`, so it is solved by Newton:
`composition_residual` builds the residual and the two face-flux slopes, and
`solve_mass_fraction` assembles the tridiagonal Jacobian and takes the
correction. The Jacobian is an M-matrix for any iterate in `[0,1]` -- the
slopes are `dJ/dX_left >= 0` and `dJ/dX_right <= 0` in all three branches --
so the correction is well posed, and the nonlinearity is quadratic, so the
iteration converges in one to three passes. A pass that does not reduce the
residual is halved.

A Picard sweep on a lagged factor was the alternative and is what was there.
It does not restore the bound, whatever the sweep count: at every sweep the
flux at a face between two cells at `X = 1` is `-W_f != 0`, and only the
*converged* fixed point has the shutoff. The two-sweep version in use was not
converged -- measured below, the LHS 1140 b step needed one to three Newton
passes to reach its residual floor -- and it is the unconverged iterate that
was clamped. Reducing `dt` was not considered a fix: the bound is a property
the operator should have at any step, the implicit operator is otherwise
unconditionally stable, and the relaxation deliberately runs at `1e5` cell
diffusion times to reach the composition's steady state at all (section 68).

**Convergence measure.** The Newton residual is measured *relative to the
size of the terms of its own row*. At a step long against the cell diffusion
time the time term and the flux divergence each exceed their difference by
`~dt/t_diff ~ 1e5`, so an absolute residual bottoms out at the round-off of
the larger term (measured: `~1e-11` in units of `X`) and cannot be a
convergence test. The relative residual reaches `~7e-12` on LHS 1140 b in one
to three passes; below that the round-off of the tridiagonal solve, amplified
by the same conditioning, is all that is left, so an already-small residual
that no longer halves is taken as converged. The smallness test is not one
number: the reachable floor is set by the conditioning of the tridiagonal
solve, which is `1e-12` in a wind and `~1e-5` in the `K_zz = 2e12` homopause
column of test T7, where the time term is negligible against the eddy term.
It is therefore either an absolute `1e-8` or a drop of `1e-6` from the
residual the step started at, whichever is reached first. That floor test is
**not**
applied while the residual is still large: the bounds belong to the converged
step, and stopping early loses them -- with the floor test applied
unconditionally, T14 fails at an excursion of `7e-3`.

### What is measured

The excursion outside `[0,1]` *before* the clip, which section 84 added to the
step diagnostic and which is now also exposed as the module variables
`he_fraction_over_one` and `he_fraction_under_zero` so a test can read the
solve instead of its clip:

| | before | after |
|---|---|---|
| T14, 400 steps of a closed column at 20 Jeans parameters, `dt = 1e6` | max `X` = **1.514** (excursion `+5.136e-1`) | max `X` = 0.99908 (excursion `-9.2e-4`) |
| T13, 200 steps from a band seeded at `X = 1` exactly | excursion `-8.584e-3` | excursion `-8.587e-3` |
| LHS 1140 b `heh2p13_diff_kzz1e9`, JFNK finish, 64-66 diffusion steps | excursion `-1.050e-1` | excursion `-1.050e-1` |

**The excursion is negative in every case now measured, i.e. the solve stays
strictly inside `[0,1]` and never reaches the clip at all.** It is *not* the
round-off-sized positive number one would quote if some cell sat exactly on
the bound: in these runs the solution approaches the bound and stops short of
it, which is the behaviour the shutoff produces. The clip lines are therefore
assertions at both ends, and they are kept, with the excursion measured every
step rather than asserted.

The old code fails T14 by 0.514 -- essentially the 0.520 recorded in section
84 from a wholly different run, which is the reproduction that identifies the
mechanism.

### Golden movement

The change moves every diffusion-on result, so it moves two of the seven
regression cases. The other five -- `wasp_full`, `wasp_he23off`,
`mol_base_handoff`, `mol_metals`, `mol_lyman_werner` -- run with
`He_diffusion` off, never enter the module and stay **byte-identical**.

Largest relative difference over the column, against the previous golden, on
the two 12000-step snapshots that do run element diffusion:

| case | rho | p | T | heat | cool | v | worst species | H2 front (f = 0.5) | log10 Mdot |
|---|---|---|---|---|---|---|---|---|---|
| `mol_diffusion` | 1e-11 % | 2e-11 % | 9e-12 % | 1e-11 % | 2e-11 % | 1.4e-8 % | `H3+` 1.9e-4 % | 1.16036 -> 1.16036 | 10.58 -> 10.58 |
| `lower_profile` | 4e-9 % | 5e-9 % | 3e-9 % | 7e-9 % | 9e-9 % | 5.4e-7 % | `HeITR` 1.3e-8 % | -- (no molecular column) | 9.58 -> 9.58 |

**These two cases do not exercise the bound at all**, and that is why they
barely move: `X` stays in `[0.2406, 0.2408]` through all 12000 steps of
`mol_diffusion` and in `[0.2468, 0.2496]` of `lower_profile`, both of which
run at `K_zz = 1e9` where the eddy term dominates the drift, and the largest
excursion the diagnostic reports is `-0.759` and `-0.750` -- i.e. `X` is
nowhere near either end. What moved them is the second-order difference
between the old and new evaluation of the drift product plus the Picard
sweep's non-convergence, and both are at round-off amplified by the step
conditioning here. Goldens re-snapshotted for the two, and `make check`
re-run after: 7/7.

Newton cost, measured on the same runs: `mol_diffusion` converges in one pass
in all 12000 steps (largest relative residual 5.9e-13); `lower_profile` in one
pass in 11807 steps and two in 193 (9.8e-13). The wall time is unchanged
within the noise of a shared machine.

### The physical result this was worried about: LHS 1140 b does not move

The composition of `LHS1140b/exhale/heh2p13_diff_kzz1e9` is the one the
He I 10830 line is read from, so the fix was re-taken on it: the case
re-solved from its own converged state, JFNK finish, post-processing pass and
the He I 10830 synthesis at WINERED HIRES-Y (`R = 68,000`), once with a
binary carrying the old linearization and once with this one, both outside the
repository so no recorded output was overwritten.

| | old linearization | this section |
|---|---|---|
| `exhale_info` / `\|\|R\|\|` | 0 / 4.285e-3 | 0 / 4.285e-3 |
| `log10 Mdot` | 7.80 | 7.80 |
| He 10830 red-pair EW (10832.60-10834.20 A vacuum) | 1.1404754 %A | 1.1404754 %A |
| red depth | 4.186959 % | 4.186959 % |
| FWHM | 0.2557731 A | 0.2557731 A |
| `X` range over the run | 0.6503-0.8950 | 0.6504-0.8950 |
| largest excursion before the clip | -1.050e-1 | -1.050e-1 |

The two solutions agree to `1.2e-9` relative in `rho`, `v`, `p` and `T` over
the whole column (the one exception is `He III` at `1.0103 R_p`, 3.5e-3
relative, on a population that is `6e-13` of the helium there). **This wind
never approaches `X = 1`**, so the clamp never fired on it and the fix has
nothing to remove; what is left is the second-order difference in the drift
term, and at this `K_zz` it is invisible.

On the reservoir crossing: the equivalent width moves by `-4.1e-13` in
relative terms, and with the local slope `d ln EW / d ln (He/H) = 0.3819`
measured over the upper rungs of the closure ladder that is
`d ln (He/H) = -1.1e-12`, i.e. the crossing stays at `He/H = 11.730000`. The
adopted result is untouched.

**What this does not say.** The recorded `tpm_He10830.txt` of that case gives
EW 1.1289 %A against the 1.1405 measured here; the difference is *not* this
change, it is everything the working tree has accumulated since that file was
written (section 84 among it). Both columns of the table above were produced
with binaries differing only in the drift linearization, which is what makes
them a controlled comparison.

### Verification

- `make check`: five cases byte-identical, `mol_diffusion` and
  `lower_profile` moved by the amounts tabled above; goldens re-snapshotted
  for those two and `make check` re-run, 7/7.
- `make diffusion_tests`: 35/35 (31 as before, plus T14's three assertions
  and one added to T13). **T14 fails on the old linearization** at an
  excursion of `+0.514`, which is the point of adding it: it drives a closed
  column hard enough that the base reaches pure helium, requires that the run
  actually reach the boundary (`max X >= 0.99`, else the assertion is
  vacuous), and reads `he_fraction_over_one` / `he_fraction_under_zero`, the
  excursions the operator measures *before* its own clip. Reading the clipped
  `X` instead would have shown 1.0 in both columns and proved nothing.
- `run_fcheck.sh` (`-fcheck=bounds,do,mem`) clean. That case runs with element
  diffusion off, so as in section 84 two further bounded runs cover the code
  this section changes: `diffusion_tests.x` built with the same flags, and the
  `lower_profile` regression case bounded to 400 steps. Both clean, no trap.

### Two more things found in passing, both in the step diagnostic

- The diagnostic printed `max|J_He + J_1|` as
  `max(sumJ, abs(Jf(j) - Jf(j)))` -- an identical zero in the clothes of a
  measurement. A binary mixture has one independent diffusive flux and the
  hydrogen component carries `-J_He` by construction, so there is nothing
  there to measure; the field is removed and the reason is written where it
  used to be printed.
- The trace-metal solver ends on `where (fX .lt. 0.0d0) fX = 0.0d0` under the
  comment "Round-off assertion, not a limiter: the M-matrix keeps fX >= 0" --
  the same unmeasured claim this section removes from the helium equation.
  That one is a different case and is left as it is: the trace equation is
  *linear* in `fX` (the trace limit drops the `1-fX` factor), and an M-matrix
  sign pattern with a nonnegative right-hand side is a real argument for
  `fX >= 0`. But it is an argument, so the excursion is now measured too
  (`trace_ratio_under_zero`, reported with the step diagnostic) instead of
  being asserted.

---

## 87. A rate coefficient that jumped by 1.57 at 4000 K: the He(2^3S) Penning rate, and the step it put in every figure (2026-08-28)

**The judgement first.** A rate coefficient is a continuous function of
temperature. The He(2^3S) + H Penning rate in `Cool_coeff.f90` was not:

```
T <= 4000 K:  1.9e-9*(300/T)**0.07     ->  1.5849e-9 at 4000 K
T >  4000 K:  9.1e-9*(300/T)**0.50     ->  2.4921e-9 at 4000 K
```

a factor **1.5724** across the seam. That is Taylor et al. (2025), ApJ
989:68, Table 2, transcribed correctly -- the code was not wrong about what
the paper says. The published fit is what cannot be used, and "it is the
published value" is not a defence of an expression that is discontinuous in
a state variable. It is replaced by the Garcia Munoz (2025) form, which is
continuous everywhere.

**Checking the fit against its own paper made the case decisive.** Taylor et
al.'s Figure 19 (bottom right) plots the quantity Table 2 is supposed to
represent: their Maxwell-Boltzmann average of the Cohen & Lane (1971) and
Morgner & Niehaus (1979) cross sections. Digitized here, that curve rises to
1.50e-9 near 1400 K and falls to 1.29e-9 at 4000 K. The tabulated low branch
decreases monotonically from 1.83e-9 to 1.58e-9 and, being a power law,
cannot produce a maximum at all:

| T [K] | Fig. 19 (their numerical rate) | Table 2 low branch | ratio | new rate |
|---|---|---|---|---|
| 500  | 1.135e-9 | 1.833e-9 | 1.62 | 1.024e-9 |
| 1000 | 1.467e-9 | 1.746e-9 | 1.19 | 1.207e-9 |
| 1400 | 1.502e-9 | 1.706e-9 | 1.14 | 1.271e-9 |
| 2000 | 1.470e-9 | 1.664e-9 | 1.13 | 1.321e-9 |
| 4000 | 1.291e-9 | 1.585e-9 | 1.23 | 1.356e-9 |

Above 4000 K the high branch is worse: 2.23e-9 at 5000 K, above *every*
cross-section determination collected in Garcia Munoz (2025) Fig. 5, whose
curves lie between about 1.3 and 1.7e-9 there. The appendix of Taylor et al.
says how this happened -- they fitted "an equation already implemented into
our KPP routine for simplicity" with `scipy.optimize.curve_fit`, in two
pieces, with nothing tying the pieces together at the join.

**The replacement.** Garcia Munoz (2025), A&A 698, A199, Fig. 5 publishes a
closed form for the Maxwell-Boltzmann average of the Movre & Meyer (1997)
total ionization cross sections:

```
k = 1e-9 exp(c/T + d1 lnT + d2 (lnT)^2 + d3 (lnT)^3)   [cm^3 s^-1]
c = -8.64804E+01  d1 = -2.86766E-01  d2 = +8.68445E-02  d3 = -5.73001E-03
```

quoted there as accurate to better than 2% from 200 to 10^4 K. It reproduces
that paper's own Table A.5 (1.02, 1.32, 1.35, 1.27 e-9 at 500, 2000, 5000,
10^4 K) to 0.4%, measured here. Note the paper's network file writes the
third term as `exp(c)`; only the `c/T` of the Fig. 5 caption reproduces
Table A.5, and the code comment records that.

This is not a case of trading one paper's fit for another's. Taylor's
Figure 19 and the Garcia Munoz curve agree to 10-20%; it was the *tabulated*
fit that was the outlier. Re-digitizing Movre & Meyer's Fig. 10 cross
sections independently and Maxwell-averaging them reproduces Table A.5 to
5-10%, the accuracy of reading a figure. Movre & Meyer give cross sections
only, over 0.01-1.0 eV, and no rate coefficient -- the 200 K end of the fit
leans on Garcia Munoz's low-energy `E^-0.2` extrapolation (15% of the
integral weight at 200 K) and the 10^4 K end on the high-energy `E^-1` one
(40% of the weight at 10^4 K), which is why the validity range is stated
where the expression is written.

The older prescriptions are outliers for identifiable reasons and are not
candidates: Roberge & Dalgarno (1982)'s temperature-independent 5e-10 (what
EXHALE used before Taylor) sits far below every cross-section calculation,
and Fassia et al. (1998)'s 7.5e-10 (T/300)^(1/2) is Bell (1970)'s single
300 K point stretched with hard-sphere scaling, which the energy dependence
of the measured cross section contradicts.

**The H2 channel and the 0.9:0.1 branching.** Both partners ionize the
metastable through two channels that share one cross section,

```
Penning:      He(2^3S) + X  -> He(1^1S) + X^+ + e^-
associative:  He(2^3S) + X  -> HeH^+ (+ fragment) + e^-
```

and Garcia Munoz splits the total "an average 0.9:0.1". The published
network file applies exactly that split through the amplitudes of two rows
per channel (0.9e-9 / 0.1e-9 for H; 4.86740e-12 / 5.40822e-13 for H2, summing
to 5.408222e-12), which settles what the paper text leaves implicit. EXHALE
had charged 100% of both totals to Penning. It now follows the split, and
the distinction matters in a way a single scale factor does not:

* the **total** removes the metastable -- both channels quench it -- so
  every sink term keeps the sum. Using 0.9 there would have made the
  metastable density 10% too high, and the metastable density is the
  observable;
* `f_penning_HeI23S = 0.9` scales only the terms that leave a lasting proton
  or H2+ behind and that deposit the Penning exothermicity;
* the molecular system carries HeH+ with its own balance row, so the
  remaining 0.1 is now a real HeH+ source there. The atomic systems have no
  HeH+ row and drop it, on the stated ground that in an H2-poor gas HeH+
  dissociatively recombines back to He + H faster than it reacts onward, so
  the branch is a metastable sink but not a lasting proton or electron
  source. Where H2 is abundant that cycle is broken by HeH+ + H2 -> H3+ + He,
  which is why the molecular system does not take the same shortcut.

The H2 rate itself is now the published `5.408222e-12 T^0.675388
exp(-696.275/T)` rather than a fit made here to the four Table A.5 points.
The old fit was accurate (0.13%) but reproduced the *total*, so charging all
of it to Penning overstated H2+ production from this reaction by 1/0.9.

**Naming.** `penning_HeI_23S` and `penning_HeI23S_H2` returned totals, not
Penning rates, and had inconsistent names for the same physics. They are now
`ioniz_HeI23S_H` and `ioniz_HeI23S_H2`, both `elemental` functions, named for
what they compute.

**This step was the discontinuity in the LHS 1140 b memo figures.** The
detector in `LHS1140b/exhale/metastable_step_diagnosis` flags a cell whose
neighbour ratio departs by more than 10% from the geometric mean of the
ratios on either side. Before and after, on the same three runs:

| run | flags before (eq / adv) | flags after |
|---|---|---|
| `heh0p55` | 6 / 6 | **0 / 0** |
| `heh2p13_diff_kzz1e9` | 6 / 6 | **0 / 0** |
| `flux_closure/heh11p1/k01` | 4 / 8 | 0 / 4 |

The four that remain in `k01` are at r = 19.15, 19.47, 28.57 and 29.05 R_p,
were present before with the same values, and are outer-boundary features,
not `Penning T=4000` ones. All three runs do cross 4000 K, so the test is
live: the inner crossing is at 1.1775, 1.1084 and 1.1123 R_p. At the
crossing cell of `heh0p55` the neighbour ratio n(2^3S)[j+1]/n(2^3S)[j] went

```
before:  1.0199  1.0194  0.6651  1.0232  1.0227      (a 33% drop in one cell)
after:   1.0190  1.0185  1.0180  1.0175  1.0171
```

and the other two runs behave the same way (0.6681 -> 1.0242, 0.6968 ->
1.0352).

**What moves.** The new coefficient is *lower* than the old one everywhere
in the wind (ratio new/old 0.69 at 1000 K, 0.86 just below 4000 K, 0.54 just
above, 0.81 at 10^4 K), so the metastable lives longer and the 10830 line
gets stronger. Removing the step is the small part of this; the change in
level is the large part.

Regression matrix, current binary against the goldens as they stand
(`max|new - golden| / max|golden|` over the column, pointwise relative
maximum in parentheses):

| case | rho | v | p | T | heat | cool | log10 Mdot |
|---|---|---|---|---|---|---|---|
| `wasp_full` | 8.6e-05 (3.3e-03) | 7.8e-04 (1.8e-02) | 4.2e-05 (1.5e-03) | 6.4e-04 (1.3e-03) | 2.6e-03 (5.8e-03) | 2.5e-03 (9.3e-03) | 13.23 -> 13.23 |
| `wasp_he23off` | 0 | 0 | 0 | 0 | 0 | 0 | 13.23 -> 13.23 |
| `mol_base_handoff` | 2.2e-06 (2.2e-02) | 1.1e-02 (2.2e+00) | 4.3e-06 (3.9e-02) | 1.5e-02 (1.7e-02) | 6.1e-03 (2.7e-02) | 4.8e-05 (3.5e-02) | 10.58 -> 10.57 |
| `mol_metals` | 1.6e-05 (2.3e-02) | 1.1e-02 (1.4e-01) | 1.3e-06 (4.0e-02) | 1.5e-02 (1.7e-02) | 5.4e-03 (2.7e-02) | 1.6e-04 (5.7e-02) | 10.58 -> 10.58 |
| `mol_lyman_werner` | 6.3e-05 (2.9e-02) | 1.5e-02 (8.4e-01) | 7.8e-06 (4.8e-02) | 1.8e-02 (2.0e-02) | 7.8e-03 (3.3e-02) | 2.0e-04 (1.1e-01) | 10.58 -> 10.58 |
| `mol_diffusion` | 2.2e-06 (2.2e-02) | 1.1e-02 (1.9e+00) | 4.3e-06 (3.9e-02) | 1.5e-02 (1.7e-02) | 6.1e-03 (2.7e-02) | 4.8e-05 (3.5e-02) | 10.58 -> 10.57 |
| `lower_profile` | 3.2e-06 (2.8e-03) | 2.5e-02 (2.6e+00) | 4.6e-06 (3.5e-03) | 3.0e-02 (3.9e-02) | 3.3e-03 (8.8e-02) | 5.8e-03 (5.4e-02) | 9.58 -> 9.58 |

`wasp_he23off` is byte-identical, which is the control: with the metastable
off the coefficient is never evaluated. The large pointwise `v` ratios are
cells where the golden velocity passes through zero, not large absolute
changes -- the scale-relative column is the one to read. The H2 front does
not move (1.0181 R_p before and after in all four molecular cases), and
log10 Mdot moves by at most 0.01 dex.

The species that moves is the one the change is about: peak n(He 2^3S) rises
by 3.1% (`wasp_full`), 6.9% (`lower_profile`) and 22-30% (the molecular cases).
Peak HeH+ rises by a factor 116-710 in the molecular cases (2.6e-1 -> 3.4e+1
cm^-3 in `mol_base_handoff`), which is the newly resolved associative branch
-- previously the only HeH+ source in those runs was He+ + H2, and this one
is larger.

**Goldens are NOT refreshed here.** The matrix is measured, not
re-snapshotted; the refresh belongs at the end of the change series.
`make diffusion_tests` 35/35 and `run_fcheck.sh` CLEAN with the change in.

**The LHS 1140 b numbers move, and the helium-abundance crossings move with
them.** Re-solving the three reference winds (the flux-closure case with the
wind re-converged on a fixed profile, no closure iteration) gives:

| run | red EW [%.A] before | after | FWHM before | after | log10 Mdot |
|---|---|---|---|---|---|
| `heh0p55` | 1.10884 | 1.30991 (+18.1%) | 0.2573 | 0.2583 | 7.76, unchanged |
| `heh2p13_diff_kzz1e9` | 1.12887 | 1.33805 (+18.5%) | 0.2553 | 0.2568 | 7.80, unchanged |
| `flux_closure/heh11p1/k01` | 1.08532 | 1.16816 (+7.6%) | 0.2825 | 0.2840 | 7.40, unchanged |

The mass-loss rate is untouched: this is a trace-level population, not an
energy-budget term. The crossings with the measured equivalent width, solved
on re-computed composition ladders rather than from the local slope (the
+18% shift is well outside the slope's range of validity), move down:

| ladder | before | after | 1 sigma |
|---|---|---|---|
| well mixed, diffusion off | 0.55 | **0.424** | 0.409-0.440 |
| K_zz = 1e9 (adopted) | 2.09 | **1.626** | 1.573-1.685 |
| flux closure | 11.73 | **9.317** | 8.792-10.050 |

The local-slope estimate (slopes 0.709, 0.921, 0.382) gives 0.435, 1.738 and
9.676 for comparison; the direct ladder solve is the number to use. Two
caveats stand on the flux-closure entry: 9.317 is an interpolation across
existing arms, not a re-run closure, and confirming it means solving a
closure arm at a reservoir near 9.3, where the cold-trap composition differs
arm by arm. The stored He/H grid runs elsewhere in `LHS1140b/exhale` are
still on the old coefficient and were not re-solved.

`docs/lhs1140b_exhale_vs_pwinds.tex` is *not* edited here. Every quantitative
statement in it that depends on n(2^3S) is now stale -- the three crossings
above, the Penning share of the metastable sink (89%/60% at He/H = 2.09/10.31)
and the lifetimes in its Table 4, its EW table, and the two places that name
the rate as "the Taylor et al. (2025) fit" (around the `Q_31` description and
the metastable-timescale table). Its figures that show the He 10830 profile
or n(2^3S) need regenerating, and `figures/lhs1140b_ew_vs_heh.pdf` in
particular carries the crossing.

**Files.** `src/modules/radiation/Cool_coeff.f90` (both rate functions
renamed and replaced, `f_penning_HeI23S` added);
`src/modules/nonlinear_system_solver/ion_residual_core.f90` (Penning branch
on the H+ row, total kept in the triplet row);
`src/modules/nonlinear_system_solver/System_HeH_mol.f90` (branch on the H+
and H2+ rows, associative branch into the HeH+ row, total on the H2 and
triplet sinks); `src/modules/nonlinear_system_solver/System_implicit_adv_HeH_TR.f90`;
`src/modules/radiation/util_ion_eq.f90`;
`src/modules/radiation/ionization_equilibrium.f90`;
`src/modules/post_process/post_process_adv.f90` (heating terms).
`docs/molecular_chemistry_audit_he_rich.md` section 4.4 updated, including
the withdrawal of its claim that the step entered the Jacobian -- the
coefficient is constant within a cell, so it appeared only as a jump in the
loss rate between neighbouring cells.
`references/garcia_munoz_2025_network/` holds the published network file the
branched amplitudes come from.

## 88. A population solved as a difference, and a post-process that never read `info` (2026-08-29)

**The judgement first.** Two defects in the advection-corrected post-process,
both numerical, both fixed here.

1. **The ground-singlet helium population was solved as a difference.**
   `System_implicit_adv_HeH_TR` carried the *summed* He I fraction as `x(2)`
   and the metastable as `x(4)`, and formed the singlet as `x(2) - x(4)`.
   Where helium is heavily ionized and the little neutral helium left sits
   mostly in the metastable, that difference loses every significant digit
   and finally evaluates to exactly zero. The summed-He I row then
   degenerates -- it no longer contains the unknown it is supposed to
   determine -- and `hybrd1` returns `info = 4` with a residual it cannot
   reduce. Fixed by solving for the two populations, `n(1^1S)/n_He` and
   `n(2^3S)/n_He`, and forming the *sum* where the total is wanted. A sum has
   no such failure mode.

2. **No return code of the post-process was ever read.** `post_process_adv`
   calls `hybrd1` three times per cell -- the advection ionization system,
   the metal stage re-solve (`pp_metals = 2`), the energy equation -- and
   wrote the returned iterate out in every case, converged or not. Fixed:
   a cell whose solve did not converge keeps the equilibrium state it came
   in with, and the number of such cell solves is reported.

The run that exposed both, `LHS1140b/exhale/heh0p55_diff_ctrl`, is
**not made physical by either fix**, and the reason is a third thing, left
untouched here and written up in section 88.4: its post-process drives
`P_HeI` to `7e15 s^-1`.

### 88.1 A population, not a difference

The advection systems integrate the ionization ODE upwind across one cell.
The triplet system solved for

```
x(1) = n_HI/n_H     x(2) = n_HeI/n_He (summed)     x(3) = n_HeIII/n_He     x(4) = n(2^3S)/n_He
```

and built the singlet, which every rate in rows 2 and 4 needs, as
`xheiS = x(2) - x(4)`. Nothing keeps the two apart. In the He-settled control
run the difference collapses pass by pass until `x(2)` and `x(4)` are the same
double, at which point `xheiS` is exactly zero, row 2 reduces to a balance
between the metastable and He II alone, and the solver stalls. Measured on
that run, before the change: **310 of 4660 advection cell solves returned
`info = 4`** (10 post-process passes, `OMP_NUM_THREADS=1`); the other three
runs checked returned `info = 1` everywhere.

The unknowns are now the two populations themselves, and row 2 is the
ground-singlet balance -- the summed He I row minus the metastable row of
row 4, which is the same system written in different variables:

```
fvec(2) = x_s,old - x_s + c1*(  x_HeII*a_HeII*n_e
                              - x_s*(P_HeI + (b_HeI + q13)*n_e)
                              + x_t*((q31a + q31b)*n_e + A31 + x_HI*Q31*n_H) )
```

Loss of the singlet: photoionization, electron-impact ionization, and
collisional excitation `q13` into the metastable. Gain: recombination of
He II into the singlet ladder (`a_HeII`; the `a_HeI23S` channel feeds the
metastable and is charged to row 4), and every route by which the metastable
comes back to the ground state -- collisional de-excitation `q31a + q31b`,
the `2^3S -> 1^1S` decay `A31`, and the `He(2^3S) + H0` ionizing collisions
`Q31`, whose Penning branch leaves `He(1^1S) + H+ + e-` and whose associative
branch makes HeH+ that dissociatively recombines back to ground-state helium.
The He--H charge-exchange pair stays where it was, on rows 1 and 2, with the
*summed* He I as the reactant density: that is the accounting the
summed-He I form carried and the one the equilibrium systems make, neither of
which gives the metastable a charge-exchange channel of its own.

`post_process_adv` now carries `nheiS` and `nheiTR` as the two profiles it
solves for and forms `nhei = nheiS + nheiTR` wherever a routine wants the
total. The difference survives in exactly one place, the entry point, where
the equilibrium solution hands over a summed He I and a metastable; it is
well conditioned there because the equilibrium metastable is bounded by its
own balance, whose loss side is dominated by `A31`.

**What the change does to the answer.** The reparameterization is a change of
variables, so a converged solution must not move. Measured, PP-only re-runs
at `OMP_NUM_THREADS=1` against a baseline built from the same tree without the
change:

| run | `info /= 1` solves, before -> after | max relative change of the `_adv` species |
|---|---|---|
| `heh2p13_diff_kzz1e9` (K_zz = 1e9) | 0 -> 0 | 2.6e-9 (He II), median 3e-14 |
| `flux_closure/heh11p1/k01` | 0 -> 0 | 2.0e-9 (He II), median 3e-15 |
| `heh0p55` (no diffusion) | 0 -> 0 | 9.7e-4 (He III), median 1.3e-4 |
| `heh0p55_diff_ctrl` (K_zz = 0) | **310 -> 0** | 0.25 (He I and 2^3S, at 2.7 R_p) |

The first two are roundoff: the solve is the same solve. The `heh0p55` entry is
*not* the reparameterization -- it is the new `info` guard of section 88.2 firing
on one cell of the energy equation in passes 5 and 7. The last is the run whose
solution was not a solution.

**Where the cancellation went.** In the control run at `r = 2.53 R_p`, tenth
pass, the same cell before and after:

```
before   x(2) = x(4) = 1.089669e-05   (singlet exactly 0)   info = 4
         fvec = (-9.5e-4, 1.22e-2, 7.7e-6, -3.5e-4)
after    x_s  = 2.53e-22   x_t = 1.098e-05  (x_s/x_t = 2.3e-17)  info = 1
         fvec = ( 4.3e-17, 2.4e-08, -3.3e-18, 9.6e-18)
```

The largest residual falls from `1.2e-2` to `2.4e-8`, and the singlet keeps its
relative accuracy seventeen decades below the metastable, where the difference
form had nothing left. One-cell steps in the `_adv` metastable profile of that
run, counted with the local-trend detector of section 87, go from **9 flagged
cells to 1** -- and that one, at `r = 29.05`, is the outer-boundary feature the
equilibrium profile carries too.

Its neighbours are the more instructive record: cells `334` and `335` returned
`info = 1` with `fvec(2) = 1.17e-2` and `1.31e-2`. `hybrd1`'s `info = 1` is a
test on the *increment* between iterates, not on the residual, so a solver
frozen against a degenerate row reports success. **The `info` check of
section 88.2 would not have caught those cells.** Only the reparameterization
does.

### 88.2 The return code

`post_process_adv` made three `hybrd1` calls and read the return code of none
of them:

| call | what it solves | on failure, before |
|---|---|---|
| `adv_implicit_H` / `adv_implicit_HeH` / `adv_implicit_HeH_TR` | the advection ionization system of the cell | the iterate was written to `nhi..nheiTR` and fed the upwind cascade of the next cell |
| `ion_system_metals_pp` | the metal stage split at fixed H/He (`pp_metals = 2`) | the iterate replaced the equilibrium split |
| `T_equation` | the advected energy equation (legacy MINPACK path) | the iterate became `T_out(j)`, subject only to the metal-mode band reject |

A non-converged iterate satisfies neither the advection balance it was asked
to solve nor the equilibrium balance it started from, so it is not a state of
the gas. The equilibrium solution of that same cell is, and it is already what
the three validity conditions of section 72 fall back to. Each call now keeps the
equilibrium state on `info /= 1` -- the ionization of the cell, the equilibrium
metal split, the converged equilibrium temperature -- and counts it. The counts
are summed over the ten passes, not read off the last one: a cell reverted in an
early pass feeds that pass's upwind cascade whether or not the last pass
converges everywhere.

Measured on four LHS 1140 b runs, the guard is quiet: the advection system
fires nowhere once section 88.1 is in, and the energy equation on one cell in
two of the ten passes of `heh0p55` and of the He 2^3S-off variant of
`heh2p13_diff_kzz1e9` -- worth `1e-4` to `1e-3` in the `_adv` profiles of those
runs. The metal guard is untested by them: `pp_metals` defaults to the frozen
mode (1), so the re-solve block is one those runs never enter.

The guard accepts `info = 1` and rejects everything else. That is a statement
about the increment, not about the residual, and section 88.1 shows the two can
disagree. A residual test would be the stronger guard; it needs a threshold in
the units of the rows (a fraction per cell crossing), which is a choice, not a
measurement, and it is not made here.

**A third, smaller thing in the same block.** With helium off
(`thereis_He = .false.`) the post-process left `nhei`, `nheii`, `nheiii`,
`nheiTR` -- and, with this change, `nheiS` -- undefined until the helium-free
branch of the ionization loop zeroes them at its *end*. The first pass reads
them before that, in `nhe`, in the electron sum and in `eval_cool`. They are
now zeroed where the helium arrays are initialized. Helium-on runs are
bit-identical across the change (checked on `heh2p13_diff_kzz1e9`, all three
`_adv` files).

### 88.3 Why nothing caught this

A run's convergence tests never look at the post-process.

* `du` and the steady residual `||R||` are momentum and energy metrics of the
  *wind*, evaluated in the time loop. The loop is over before
  `post_process_adv` is called; there is no metric downstream of it.
* `ionization_equilibrium` validates the root it accepts -- physical simplex,
  up to three starting points, `ok_rank`, and a one-line sweep summary of
  reseeds, retries, unphysical roots and failures (section 45). The
  advection-corrected post-process, solving a system of the same size with the
  same solver, validated nothing.
* `make check` bit-compares `Hydro_ioniz.txt` and `Ion_species.txt`. Those are
  the *equilibrium* files, written before `post_process_adv` runs. **No golden
  case compares an `_adv` file at all**, so nothing in the regression can move
  when the post-process does -- which is also why the change documented here
  leaves the goldens byte-identical (7/7, 14/14 files), and why that result
  carries no information about the `_adv` profiles. The A/B tables above are
  what stands in for it.

That leaves the `_adv` files checked only by whoever plots them. The two
defects were found that way: an order-of-magnitude one-cell step in a
metastable profile.

### 88.4 What the control run actually shows, and what is not fixed here

With section 88.1 in place the solver converges everywhere in
`heh0p55_diff_ctrl` and still returns a neutral-helium column five orders of
magnitude below the equilibrium one. The singlet is not being lost to
cancellation any more; it is being destroyed by the rate the post-process
hands the solver. Traced pass by pass at `r = 2.53 R_p`:

| pass | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| `P_HeI` [s^-1] | 2.5e-4 | 3.7e-2 | 5.3e0 | 1.1e5 | **6.9e15** |
| `x(1^1S)` | 6.8e-3 | 4.7e-5 | 3.3e-7 | 1.6e-11 | 2.5e-22 |

A photoionization rate of `7e15 s^-1` is not a rate. Two things produce it,
and neither is touched here.

**(a) A volumetric secondary-ionization rate divided by the density of the
species it ionizes.** `util_ion_eq.f90` adds the Shull & van Steenberg (1985)
secondary channel as

```
P_HeI(j) = P_HeI(j) + R_secHeI/max(nheiS(j), 1.0d-99)
```

`R_secHeI` is the energy deposited per unit volume times the SvS85 branching
into He I ionization: a *volumetric* rate, and one carrying SvS85's own fixed
helium abundance. Dividing it by the actual `n(1^1S)` gives the right per-atom
rate only at that abundance, and diverges as the singlet disappears. In this
run helium has settled to `n_He/n_H = 6e-6` at `1.36 R_p`, four decades below
the SvS85 composition, and the ten post-process passes close the loop: a
smaller singlet raises its own destruction rate, which shrinks it further. The
`max(..., 1e-99)` is what stops the divergence, at a rate of order `1e60`.
The physical scaling is the opposite one -- the secondary He I ionization rate
per unit volume should fall with the helium abundance, roughly as
`n_HeI sigma_HeI / sum_a n_a sigma_a`, leaving the per-atom rate near its
solar-composition value.

**(b) In a `Do only PP: True` run the equilibrium solve and the post-process
use different physics.** The SvS85 coupling is staged: `sec_ion_active` starts
false and is switched on once the wind has first converged (section 38).
`EXHALE_main.f90` also sets it unconditionally *after* the time loop, so the
post-process always runs with the full coupling. With `Do only PP` the loop
exits after one iteration, so the single `ionization_equilibrium` call runs
with `sec_on = F` and all twenty-one `PH_heat_HHe` calls of the post-process
run with `sec_on = T` (printed and counted). The equilibrium `P_HeI` at
`1.36 R_p` of the control run is `2.4e-7 s^-1`; the post-process's first pass,
at the same densities, `4.3e-5 s^-1`. The `_adv` files of such a run are not a
correction to the equilibrium files next to them: they are a different model.

**(c) The same difference, in the radiation.** `util_ion_eq.f90` builds its
own singlet as `nheiS = nhei - nheiTR` and uses it for the He I column
density, the He I photoionization integrand and the denominator in (a). The
state vector of the whole code carries (summed He I, metastable) and
reconstructs the singlet wherever it is needed -- `species_table.f90` says so
in as many words. Section 88.1 fixes the one system where the difference was
catastrophic; carrying `n(1^1S)` as the state variable throughout would touch
the global species vector, the output schema and every reader of it, and is
not attempted here.

Six of the 354 stored `LHS1140b` runs carry the signature of (a) --
`n_HeI <= n(2^3S)` in the `_adv` file above `1.5 R_p`:
`heh0p55_diff_ctrl`, `heh0p55_diff_d3`, `heh0p55_diff_fix` (171 cells of 233
each), `heh0p55_diff_kzz1e6` (36), `xuv0p01_heh2p13` (25) and
`xuv0p10_heh2p13` (28). All six are helium-settled: `K_zz <= 1e6`, or a
reduced XUV flux. The `K_zz = 0` arm of Fig. `lhs1140b_kzz_profiles.pdf` and
its quoted peak `n(2^3S) = 0.046 cm^-3` are read from the first of them.

**Files.** `src/modules/nonlinear_system_solver/System_implicit_adv_HeH_TR.f90`
(the two populations as unknowns, the ground-singlet row);
`src/modules/nonlinear_system_solver/System_implicit_adv_HeH.f90` (the same
naming, no arithmetic change: without the metastable all He I is the singlet);
`src/modules/nonlinear_system_solver/ion_cell_state.f90` (`xhei_old` ->
`xheiS_old`); `src/modules/post_process/post_process_adv.f90` (`nheiS` and
`nheiTR` as the profiles solved for, `nhei` formed as their sum, the three
`info` guards, their counters and reports, and `pin_cell_to_equilibrium`,
which the validity conditions and the guards share).
`docs/postprocess_advection_validity.md` carries the same items next to the
2026-08-12 pair. `docs/lhs1140b_exhale_vs_pwinds.tex` is *not* edited.

Found in passing and fixed, outside the post-process:
`LHS1140b/make_memo_figures.py` and `LHS1140b/exhale/kzz_scan_table.py` formed
the elemental helium density as `HeI + HeII + HeIII + HeITR`, counting the
metastable twice -- it is a level of He I and is already inside the `HeI`
column, which `species_table.f90` states in capitals. Both now sum the three
ion stages. Measured on the four `K_zz` arms of the memo figure, the quoted
elemental ratios move by at most `2e-6` relative (`0.318246 -> 0.318244` at
`5 R_p`, `K_zz = 1e10`), so nothing was regenerated.

## 89. A conservation check that could not see conservation: the deepest-level elemental test and the reservoir band it refused (2026-08-29)

> **Read with section 90.** The diagnosis below stands, but its two decisions
> were superseded the same day. Photochem was corrected at the source, so the
> `--abundance-tol` default went `1.0e-3` -> `1.0e-10` and the three `clima`
> failures noted at the end were fixed (two others appeared).

**The judgement first.** The profile adapter's deepest-level elemental check
was refusing converged LHS 1140 b columns on **the residual of Photochem's
chemical-equilibrium solver**, not on any error of the photochemistry, the
transport or this code's carrier summation. At the level the check is stated,
Photochem pins every gas species of the handed-over mixture with
`bc_type='press'` (`gasgiants._initialize_atmosphere`), so that level carries
its equilibrium initialization unchanged and **no photochemical or transport
leak can reach the check at all**. What it reads instead is how far
`composition_at_metallicity` got, and that routine states
`# Do not enforce convergence.` The `--abundance-tol` default is raised from
`1.0e-4` to `1.0e-3` -- an interim value for the solver as it then was; section
90 tightens it to `1.0e-10` once the solver is fixed. Full diagnosis, with
every measurement: `docs/deep_level_elemental_check.md`.

**The symptom.** Reservoir He/H = 8.5, 8.55, 8.6, 8.65 and 8.7 could not be
run in the closure ladder of `LHS1140b/exhale/crossings_gm25/`: the N/H ratio
at the deepest level departed from the input by 1.02e-4, 1.22e-4, 1.43e-4,
1.67e-4 and 1.94e-4 against the 1.0e-4 tolerance. Two closure lineages with
different trial escape fluxes give the same numbers to the three digits the
refusal prints, and a standalone chemistry-only scan reproduces them a third
time -- as section 2 of the diagnosis predicts, the trial flux enters only at
the upper boundary and the refused level is Dirichlet.

**Measured, at He/H = 8.6.** Chemical equilibrium plus the quench overwrite,
the photochemical model's initial state, and the converged steady state agree
at the deepest level to **all 13 digits printed** (`X_N` =
`8.185823483630e-05` in all three; deep `q_NH3` `8.995632484e-06` converged
against `8.995632e-06` initial). Called on its own at the same P and T,
`equilibrate.ChemEquiAnalysis.solve` returns `converged = True` at every
reservoir swept while its own `molfracs_atoms_gas` misses the requested N/H by
up to 1.72e-4.

**Why nitrogen and nothing else.** In all 38 runs of the scan
`dev_N/dev_He = 4072` (4069-4086 over the whole set, 4072 in 32 of them),
against `1/(3 N/H) = 1/(3 x 8.18465e-5) = 4073`. Deep nitrogen is essentially
all NH3 (on the He/H = 8.6 profile the N2 nuclei share of `X_N` is 1.0e-5 and
HCN's is 2.2e-25), so **one excess N rides on three excess H**: He, C and O
carry the common hydrogen error alone, N carries it plus the relative error of
the equilibrium NH3 mole fraction. The check was testing that one mole
fraction 4072 times more tightly than anything else, which nobody chose.

**Not a branch change, and not a special band or planet.** Across
8.45 -> 8.75 the deep H2, H2O, CH4, NH3, N2 and He all move monotonically by
less than 1 % with no jump; the residual jumps by a factor 38 between 8.71
(1.99e-4) and 8.72 (5.2e-6), 1500 between 9.15 and 9.20, and 5800 between
10.25 and 10.30. At the old tolerance the same scan also refuses He/H = 6.0
(1.31e-4) and 10.20-10.25 (1.5-2.0e-4); the tree already held refusals at
He/H = 14.643, 15.090 (**2.12e-4**, the largest recorded) and 25.0
(`LHS1140b/exhale/kzz_profile_scan/`); and at **solar composition** the same
solver misses N/H by 1.50e-4 at 800 K and 1.14e-4 at 1600 K. The two
HD 209458 b departures that fixed the original default (1.6e-7, 2.7e-5;
section 77) were a fortunate sample, not a floor.

**No accounting error.** Recomputed from the written profile, the diagnostic
columns account for **100.00 %** of `X_N` at the deepest level. They account
for 93.9 % at the 1 microbar match, and any earlier reading of a shortfall in
the carrier sum refers to *that* level and to the deliberately short `q_*`
diagnostic column set -- the check reads neither.

**The change.** `src/utils/lower_profile_schema.py`: `--abundance-tol` default
`1.0e-4` -> `1.0e-3`, and both the option help and the comment at the check
now state what the level can and cannot test. 1e-3 leaves a factor 4.7 over
the worst departure on record and gives up next to no detection power: each
element sits in one dominant carrier there, so a real miscount is an O(1)
error, and the smallest conceivable one -- dropping N2 -- is 1.0e-5 of `X_N`
in the band at issue and was already invisible at 1e-4. Over the 38 profiles
the N2 nuclei share runs 4.2e-6 (He/H = 12, deep T = 379 K) to 4.0e-3
(He/H = 1, deep T = 588 K), rising with the deep temperature, so where
dropping it would be a real loss it stays above 1e-3. Nothing else was
changed: not the check, not its location, not the refusal policy. **Superseded
by section 90**: the headroom argument was a statement about how far the 0.8.4
equilibrium solver missed, and once that was fixed the default went to
`1.0e-10`.

**What it means for the LHS 1140 b crossings** (`crossings_gm25/results.txt`).
The flux-closure crossing at He/H = 9.048 is untouched -- the 9.0 and 9.1 arms
both ran and bracket it -- and so is the -1 sigma end 10.77, bracketed by 10.6
and 11.0. The **+1 sigma end 8.49 was interpolated across the un-runnable
band**, a 3.4 % gap in reservoir between the 8.45 and 8.75 arms; with the
tolerance raised those five rungs can be run and that end replaced by
measurement. That was done, on the corrected build: section 90. Because
the refusals at 6.0, 10.20-10.25, 14.643, 15.090 and 25.0 are the same trap at
other reservoir values, and nothing in the pattern predicts where the next one
falls, any future ladder can meet it anywhere.

**Seen in passing, not investigated.** `clima`'s
`surface_temperature_bg_gas` root solve fails outright at three of the 41
reservoir values attempted -- He/H = 0.0969 (with `Failed to compute heat
capacity`), 8.1 and 15 -- so those compositions produce no column at all.
That is a **different** failure from the refusal above (this one happens
before a column exists, the refusal after), and it is recorded as item (M) of
`TO_BE_DONE.md`. All three are fixed by the bounded, scaled climate solves of
section 90; the corrected build refuses two other compositions instead,
He/H = 9.4 and 9.5, and section 91 shows that is the same underlying solver
defect moving address rather than a new one -- 0.8.4 refuses two compositions
of its own on the same grid.

**Files.** `src/utils/lower_profile_schema.py`;
`docs/deep_level_elemental_check.md` (the diagnosis);
`LHS1140b/exhale/nh_refusal_diagnosis/` (the measurement record: README, four
scripts, five tables, 38 adapter runs under `scan/`). No EXHALE run was made
for any of it and no golden was touched; the adapter's chemistry is unchanged
and `make check` does not exercise this path.

---

## 90. The equilibrium residual removed at its source, and the observable that did not move (2026-08-29)

**The judgement first.** The residual section 89 diagnosed was fixed where it
came from -- in Photochem's chemical-equilibrium path, not by relaxing this
code's check -- and **it does not move the observable**. Re-measured end to end
on the corrected build with 11 arms, the LHS 1140 b flux-closure crossing of
the He I 10830 equivalent width is at He/H = **9.0484** (chord) / **9.0461**
(quadratic), the same to every digit the photochem 0.8.4 ladder printed; on
seven like-for-like arms the equivalent width moves by at most `1.2e-6` in
relative terms and `log10 Mdot` is unchanged in its fourth decimal. The one
EXHALE change is that `--abundance-tol` goes back to its design value,
`1.0e-3` -> `1.0e-10`. Records:
`LHS1140b/exhale/crossings_pc090/results.txt`,
`docs/photochem_solver_modification_implementation.md`.

**What was changed, and where.** Nothing in the EXHALE Fortran. The corrected
chemistry is a modified Photochem `v0.9.0` source tree, built as
`photochem 0.9.0` into a dedicated conda environment:

- **Equilibrate** now tests each positive requested elemental abundance against
  its own relative tolerance (`abs(b_0(i_atom))*mass_tol`) instead of a common
  absolute threshold scaled by `MAXVAL(b_0)`, and no longer excludes elements
  at or below `1e-6`. This is the actual cause of section 89.
- **`gasgiants.py`** exposes `equilibrium_mass_tol` (default `1e-12`), assigns
  it to `ChemEquiAnalysis.mass_tol`, verifies the returned elemental vector
  against the requested one, raises rather than silently keeping the last
  result when all five temperature retries fail, and restores `use_prev_guess`
  in a `finally` block.
- **clima** `v0.7.5` gets a bracketed background-pressure solve, a
  `log10(T_surface - T_trop), log10(T_trop)` temperature parameterization, and
  normalized residual acceptance, replacing unconstrained MINPACK solves.

Both dependency changes are reproducible patches applied through Photochem's
CMake configuration. Diagnosis:
`docs/photochem_solver_modification_investigation.md`.

**A blocking defect, found only by running a column.** The first version of the
strict closure test **rejected every composition**, including ones that ran on
0.8.4. It compared the requested elemental vector with `gas.molfracs_atoms`,
the total over gas and condensate, at every level. With
`rainout_condensed_atoms=True` the vector handed to the next level is the
gas-only one, and `gas.molfracs_atoms` started from the previous level through
`use_prev_guess` still carries that level's condensate: measured on the LHS
1140 b column at He/H = 9, the total reads `1.1e+2` against the request where
the same state solved from scratch closes to `5.7e-14`. Identical on 0.8.4 and
on the patched build, so it is Equilibrate's existing reporting behavior, not
something the Fortran patch created. The test now forms the residual **only at
a level whose `molfracs_atoms_condensate` is entirely zero** -- which is every
level deep enough to set an elemental handoff, and where the residual it exists
to catch appears. The reason the defect survived the first round of tests is
worth recording: those tests ran **one level for each of the 38 saved deep
states** and never carried a column, so `use_prev_guess` chaining and
condensing levels were never reached.

**Measured.** Re-solved standalone, the 38 saved deep states of
`LHS1140b/exhale/nh_refusal_diagnosis/scan/` give a largest relative elemental
residual of `3.539391003e-13` and a median of `1.915134717e-14`, against
`2.018197773e-4` for the same measure on 0.8.4 -- all 38 converged and
independently accepted. Over the 68 handoff writes the corrected build made
under `LHS1140b/exhale/` (43 in `scan_pc090/`, 25 in `crossings_pc090/`), the
largest deepest-level departure is `3.28e-13` and the median `1.73e-14`.

**The tolerance.** `src/utils/lower_profile_schema.py`: `--abundance-tol`
default `1.0e-3` -> `1.0e-10`, a factor 305 above the worst handoff write
measured. Detection power is untouched -- each element sits in one dominant
carrier at that level, so a genuine miscount is an O(1) error, and the smallest
conceivable one, dropping N2, is 1.0e-5 to 4.0e-3 of `X_N`. **A 0.8.4 handoff
is refused at this default, and that is the intended signal**: over the 43
stored 0.8.4 handoffs of `crossings_gm25/` the departures run 3.48e-8 (best) to
8.50e-5 (worst), so even the best of them misses `1e-10` by a factor 348.
Reproducing a stored 0.8.4 result needs an explicit `--abundance-tol` on the
command line, which is how a deliberate reproduction states that it accepts
that residual; the scans under `nh_refusal_diagnosis/` already pass
`--abundance-tol 1.0` and are unaffected.

**What it does to the LHS 1140 b crossings.**

| quantity | 0.8.4, 16 arms | corrected build, 11 arms |
|---|---:|---:|
| crossing, chord | 9.0484 | 9.0484 |
| crossing, quadratic | 9.0461 | 9.0461 |
| 1 sigma band | 8.4923 -- 10.7716 | 8.4649 -- 10.7357 |

The crossing does not move. The handoff profile does: at He/H = 9 the
deepest-level N/H departure goes from 2.1e-5 to 1.2e-14 and trace species such
as `q_H2O`, `q_CO` and `q_HCN` move by up to 1e-2 relative -- none of which
reaches He I 10830, a line set by H, He and the temperature. **The 1 sigma band
moves for a different reason than the solver**: its low end used to be
interpolated across the un-runnable band He/H = 8.5--8.7, which now runs
(departures 1.4e-14 to 3.2e-14), so real arms at 8.5 and 8.6 carry that end and
it goes 8.4923 -> 8.4649; the high end goes 10.7716 -> 10.7357 because the 10.6
arm was not seeded the same way in the two ladders, a 1.1e-3 seed difference in
equivalent width, the size of the seed spread the stored ladder already
measured.

**Also opened: solar composition.** He/H = 0.0969 is runnable end to end for
the first time. It was blocked not by the tolerance but by the climate step,
which now solves (deep boundary 720.5 K, against 432.4 K at He/H = 8.1 and
393.1 K at 15 -- the other two recovered compositions).

**Two compositions the corrected build refuses, and what they turned out to
be.** He/H = 9.4 and 9.5 raise `Could not bracket background pressure in
make_profile_bg_gas` on this build and produce columns on 0.8.4. This section
first recorded that as a regression of the new bracketing scan, with the cause
undiagnosed. **Section 91 measured it and that reading was wrong.** The scan is
not the cause -- the root is inside the interval it searches, with one
monotone sign change in the last scan step -- and the same defect is present in
0.8.4, which fails at 9.37 and 9.45 on the same grid where this build fails at
9.40 and 9.50. Both builds are tripped by a forward-difference Jacobian formed
at a step below the flux residual's round-off jitter. **The net accounting of
"three failures removed, one introduced" written here does not hold**: it
compared against 0.8.4 columns that had not been run at the matching reservoir
values, and 9.45 in particular fails on 0.8.4. Neither value is a bracket arm
of the crossing. Diagnosis, repair and status: section 91; item (M) of
`TO_BE_DONE.md`.

**Scope.** No EXHALE Fortran source or binary changed, no golden was refreshed,
and `make check` was not run: this path is not in the regression matrix. Not
tested: condensate-rich and sulfur-network equilibrium sweeps -- note that
after the fix above a condensing level is accepted on the solver's own
`converged` flag, so the closure test does not speak for such a level at all.
The stored 0.8.4 ladder was not re-run, and He/H = 9.5 cannot be repeated on
this build. The 0.8.4 environment is untouched and remains the
reproduction path for every profile whose header reads
`# source_version photochem 0.8.4`.

**Files.** `src/utils/lower_profile_schema.py` (the default);
`docs/photochem_solver_modification_investigation.md` (diagnosis),
`docs/photochem_solver_modification_implementation.md` (what was built and
measured), `docs/deep_level_elemental_check.md` sections 8--10 (the tolerance,
the crossings, the climate solve); `LHS1140b/exhale/crossings_pc090/` (the
11-arm re-measurement) and `LHS1140b/exhale/nh_refusal_diagnosis/scan_pc090/`
(the 45-point reservoir scan on the corrected build). The Photochem changes
live outside this repository, in the local `photochem/` tree.

---

## 91. The clima bracket was not the cause: a Jacobian built from round-off, and a failure that predates the patch (2026-08-29)

**The judgement first.** The two compositions section 90 reported as a
regression are not one. `clima`'s new bracketing scan does not fail to find a
root at He/H = 9.4 and 9.5 -- the root is inside the scanned interval, and the
residual crosses zero once, monotonically, inside the last scan step. The scan
fails because the *outer* `surface_temperature_bg_gas` solve hands it a
surface temperature of `6.5e8` K, where every scan point is thermodynamically
invalid. That step comes from a forward-difference Jacobian whose difference
step is smaller than the residual's own round-off jitter, and **the 0.8.4 build
has the same defect at the same rate** -- it fails at 9.37 and 9.45 where the
corrected build fails at 9.40 and 9.50. The patch changed which compositions
are unlucky, not how many. Record:
`LHS1140b/exhale/clima_bracket_diagnosis/` (README, 12 scripts, six sweep and
variant tables), measured 2026-08-29 against both installed builds with nothing
in `photochem/` modified.

### 91.1 The bracket is not the cause

`make_profile_bg_gas` scans 49 points over `log10 P` spanning 12 decades below
the requested surface pressure, so for the LHS 1140 b ladder's 20 bar deep
boundary the interval is `log10 P = [-4.6990, 7.3010]`. At He/H = 9.4 the 0.8.4
solution sits at `log10(P_bg) = 7.278533`, inside it. So does every neighbor:
over 9.30--9.60 the 0.8.4 background pressure moves only from `7.278297` to
`7.278990`, always inside (`sweep_old_9p30_9p60.txt`). Inside one failing call
the residual is monotone with exactly one sign change, in the last scan step
(`scan_residual.py`) -- a grid of 49 points cannot step over it.

What the scan actually reports is that it had **no valid point at all**. The
trial surface temperature it is called with is `6.5e8` K (`outer_trace.py`), and
at that temperature `make_profile` refuses all 49 points with `Failed to
compute heat capacity`: the thermodynamic polynomials stop at 6000 K
(`T_ceiling.py`). With no valid point there is no sign change to find, and the
message `Could not bracket background pressure in make_profile_bg_gas` is the
only thing the routine can say.

### 91.2 The Jacobian is round-off

`surface_temperature_bg_gas` solves the two-component system (flux closure,
skin temperature) with MINPACK `hybrd1`, which fixes `epsfcn = 0`; `fdjac1`
then forms `eps = sqrt(max(epsfcn, epsmch))` and steps `h = eps*|x|`. At the
ladder's start, `T_surf = 400` K and `T_trop = 120` K, that is `3.6e-8` in the
log variable, a temperature step of `2.4e-5` K. Over that step the normalized
flux residual moves by `3.8e-8`, and its own jitter is `~3e-8`
(`jacobian_noise.py`, `noise_source.py`) -- the difference is the noise, not the
derivative.

The resulting entries have the wrong magnitude and, at these compositions, the
wrong **sign**: at He/H = 9.4 the cross term `df1/dx2` comes out `+52.1` where a
resolved step gives `-0.553` (`first_step.py`). Because that cross term is tens
of times the diagonal `df1/dx1 = -1.21`, an error in the small, noisy diagonal
is amplified into a first step of several decades in `log10 T_surf`, which is
how a 400 K guess becomes `6.5e8` K in one iteration.

### 91.3 The failure predates the patch; only its address changed

The same measurement on 0.8.4 (`first_step_old.py`) shows the identical
mechanism at He/H = 9.37: its noisy Jacobian sends the first trial to
`T_surf = 0.10` K, which is exactly the `T_surf is less than T_trop` refusal
that build reports there. The patch reparameterized the solve as
`log10(T_surf - T_trop)`, which makes that particular refusal unreachable, so
the same bad step now surfaces at the other end as a bracket message. The
failure rates are the same:

| grid | 0.9.0 + patch | 0.8.4 |
|---|---|---|
| He/H 9.30--9.60, step 0.01 (31 points) | 9.40, 9.50 | 9.37, 9.45 |
| He/H 8.00--11.00, step 0.05 (61 points) | 9.40, 9.50, 10.45 | 8.10, 9.45, 9.85, 10.15 |

(`sweep_new_9p30_9p60.txt`, `sweep_old_9p30_9p60.txt`, `sweep_new_wide.txt`,
`sweep_old_wide.txt`.) **The net accounting of "three failures removed, one
introduced" that sections 89, 90 and `deep_level_elemental_check.md` carried
does not hold** -- it counted the corrected build's failures against a 0.8.4
column that had not been run at the same reservoir values. In particular the
9.45 entry those tables record as `not attempted` on 0.8.4 is a **failure** on
0.8.4.

Nothing about the composition marks a failing point. The solution is smooth
across them: the deep temperature falls monotonically from 422.72470 K at
He/H = 9.30 to 420.61128 K at 9.60, with the failures at 9.40 and 9.50 sitting
on that trend. The failure is a discrete event of the solve, not a boundary in
the physics.

### 91.4 The workaround, before the repair was installed

`--climate-t-deep-guess` (default 400.0, in
`src/utils/photochem_to_lower_profile.py`) is the only user-facing control over
where the outer solve starts. At both 9.4 and 9.5, every guess tried other than
exactly 400 K -- 300, 350, 380, 420, 450, 500, 600 -- solves, and each returns
the 0.8.4 answer to five decimals (`tguess_scan.py`). A composition that
refused could therefore be recovered without rebuilding. The installed
environment no longer refuses any composition on either grid (91.5), so the
key is a diagnostic rather than a workaround.

### 91.5 The repair, tested in replication and then installed

`fix_variants.py` drives the patched residual through the same MINPACK `hybrd`
from Python, reproducing the installed build's converged deep temperature to
`5e-5` K wherever it solves, and then varies one thing at a time over the
61-point grid:

| variant | solved of 61 | failed or stalled |
|---|---:|---|
| as patched | 54 | 9.40, 9.50, 9.80, 10.05, 10.25, 10.70, 10.80 |
| initial trust-region factor 1 | 54 | the same seven |
| trial `T_surf` confined to `(T_trop, 6000 K)` | 55 | 8.50, 9.55, 9.80, 9.85, 10.20, 10.35 |
| `epsfcn = 1e-4` | **61** | none |

`epsfcn` is not the step itself: with `eps = sqrt(max(epsfcn, epsmch))` and
`h = eps*|x|`, `epsfcn = 1e-4` is a relative step of `1e-2` in the log variable,
16.2 K at the starting point against `2.4e-5` K at the default. The Jacobian it
produces is smooth and nearly composition-independent,
`[[-1.21, -0.55], [5.0e-4, -3.557]]` across the whole 9.30--9.55 range, against
entries running from `+52` to `-42` at the default step. Applied to 0.8.4's own
parameterization it also solves all 61 points (`old_epsfcn.py`), so the repair
is independent of what the patch changed. The recovered temperatures are
422.00813 K (9.40), 421.30374 K (9.50) and 415.15428 K (10.45) -- on the smooth
trend, and equal to the 0.8.4 values at 9.40 and 9.50 to all printed digits.

Confining the trial temperature to the range of the thermodynamic data helps
but does not fix the solve on its own -- it converts one failure mode into
another, because the tropopause variable stays unbounded.

**Installed, and measured there.** In `clima` the repair is `hybrd` in place
of `hybrd1` in `AdiabatClimate_simple_solver`, the routine both
`surface_temperature_bg_gas` and `surface_temperature` drive, so that `epsfcn`
can be set at all; `epsfcn = 1e-8`, `factor = 100`, every other setting
`hybrd1`'s. `1e-8` and not the `1e-4` of the table above: a relative step of
`1e-4` in the log variable, 0.16 K at the start, is about `1e4` above the
jitter and still reproduces a converged central difference to 0.12 percent on
the flux row, where `1e-4` is a 16 K step whose flux row is 9 percent off the
derivative. Both solve the grid; `1e-8` is the middle of the window in which
neither error source is active. The reasoning is in the patch at the
constant's declaration.

The installed library carries it. Read out of it --
`AdiabatClimate_simple_solver` disassembled, the stack arguments of its
MINPACK call followed to the constants they point at -- it calls
`__minpack_module_MOD_hybrd` with `epsfcn = 1e-08` and `factor = 100.0`
(`LHS1140b/exhale/clima_epsfcn_installed/installed_epsfcn.py`). Measured on
that build, not in replication:

| grid | installed | as patched, before |
|---|---:|---:|
| He/H 8.00--11.00, step 0.05 (61 points) | **61 / 61** | 54 / 61 |
| He/H 9.30--9.60, step 0.01 (31 points) | **31 / 31** | 29 / 31 |

Every composition either build refused now solves -- 8.10, 9.37, 9.40, 9.45,
9.50, 9.80, 9.85, 10.05, 10.15, 10.25, 10.45, 10.70, 10.80 -- and the
temperatures the replication predicted are what the installed build returns:
422.00813 K at 9.40, 421.30375 at 9.50 and 415.15429 at 10.45, against
422.00813, 421.30374 and 415.15428 above, agreeing to the last printed digit.
The deep temperature is monotone in composition across the whole grid,
433.26144 K at 8.00 falling to 411.97203 K at 11.00.

**What it moves.** Nothing in the observable. The two arms of the LHS 1140 b
elemental-flux-closure ladder that bracket the He I 10830 crossing were re-run
end to end on the corrected environment against the stored arms of
`LHS1140b/exhale/crossings_pc090`, which are the same photochem 0.9.0 source
without the repair -- same seed, same parent, same binary, same transit
configuration. Both closure iterations converge at k = 1 on residuals equal to
the stored ones to all five decimals, and the equivalent width moves by
`5e-8` in relative terms:

| arm | EW without the repair | with it | relative | log10 Mdot |
|---|---|---|---:|---|
| 9.0 | 1.1068538680 | 1.1068539227 | +4.94e-08 | 7.39 -> 7.39 |
| 9.1 | 1.1092155021 | 1.1092155563 | +4.88e-08 | 7.39 -> 7.39 |

The crossing moves from 9.0484199298 to 9.0484176256 (chord) and from
9.0460519944 to 9.0460494223 (quadratic), `2.3e-6` and `2.6e-6`, so the
`9.0484` and `9.0461` of record stand. Every quality indicator of that ladder
is unchanged digit for digit. The climate solve itself moves a composition
that already solved by `2.1e-6` K at He/H = 9.0 (424.9534068835 against
424.9534089414 K). Record:
`LHS1140b/exhale/clima_epsfcn_installed/results.txt`.

**The source tree carries the repair; a clone outside this repository does
not.** `photochem/src/dependencies/patches/clima-bounded-scaled-solvers.patch`
holds `real(dp), parameter :: epsfcn = 1.0e-8_dp` (md5 `3c3c40cc`), which is
also what `src/utils/photochem_exhale.patch` writes, and
`src/utils/setup_photochem.sh` checks that line by name and passes. The tree
reproduces the installed library. What is a revision behind is the separate
clone at the workspace level, `ExoAtmosphere/photochem/` (md5 `6734053e`),
where the modification was first developed: it still writes `hybrd1` and has
no `epsfcn`. That clone is outside this repository and nothing here reads it.

### 91.6 What is not established

The `~4e-7` relative jitter in `ISR - OLR` was not traced to a particular part
of the radiative transfer. It is not the inner pressure solve: the background
pressure returned is bit-identical across the steps that produce the difference
(`noise_source.py`). Whether the same trap is reachable at other planets, other
`K_zz` or other deep pressures was not tested -- only LHS 1140 b's
configuration was scanned. No EXHALE Fortran source, binary or golden was
touched, and `make check` does not exercise this path.

**Files.** `LHS1140b/exhale/clima_bracket_diagnosis/` (the diagnosis);
`LHS1140b/exhale/clima_epsfcn_installed/` (the repair measured on the
installed environment); `src/utils/photochem_exhale.patch` (the current Clima
patch, `epsfcn` included); `src/utils/photochem_to_lower_profile.py`
(`--climate-t-deep-guess`); item (M) of `TO_BE_DONE.md`.

---

## 92. A column that came back 2.4 % apart in pressure: the grid the solver stops on, not the temperature it started from (2026-08-29)

**The judgement first.** The handoff-profile sensitivity left open at the end
of section 91 -- a deep boundary temperature 2.1e-6 K apart giving a column
2.4 % apart in pressure, 15 % in one trace species and 0.5 % in the elemental
O/H -- is **the extent of Photochem's altitude grid at the instant the run is
declared steady, and nothing else**. It is numerical, it is not an
amplification of 2.1e-6 K, and it is not the chemistry's convergence tolerance.
**It does not limit the closure chain**: the base state EXHALE reads at the
matching level reproduces to 6.5e-08 in the He/H the closure iterates on,
against `element_flux_closure.py`'s `tol = 0.05`. One defect was found on the
way and corrected: `insert_level` broke `n = p/(k_B T)` at exactly one level of
the table, the matching level. Full diagnosis, with every measurement:
`docs/lower_profile_deep_boundary_sensitivity.md`; record
`LHS1140b/exhale/deep_temperature_sensitivity/`.

**Reproduced inside one build.** The initial guess handed to the climate root
solve is not physics -- it is where MINPACK starts -- so ten guesses over
300-450 K are ten solves of the same problem. They return the deep boundary
temperature over a range of **2.121e-06 K** (424.95340688 to 424.95340876),
the tropopause to 2.5e-08 relative and the water above it to 6.6e-08; re-run in
a fresh process the ladder is bit-identical. Seven full adapter runs at seven
of those guesses reproduce every number section 91's record left open: the
pressure at fixed level index spreads by **2.470e-02**, `q_OH` by 0.58, and the
oxygen handed over at the matching level by **5.207e-03**. `g400a` and `g400b`
are the same command run twice and write byte-identical files, so none of this
is run-to-run scatter.

**Everything is two-valued, and the label is the model top.** Over all eleven
runs the profile top is either `1.0934e-08` bar or `1.12075e-08` bar and never
between; the oxygen at the match is either `8.8461e-07` or `8.8892e-07`. Within
a group the columns agree to 1e-06 or better. The two groups are two grids:
with the inserted node removed so the indices line up, the pressure difference
is **zero at the pinned bottom** (16.73098 bar in both) and grows linearly in
level index to 2.496e-02 at the top -- the signature of a grid laid out
uniformly in ALTITUDE, `dz = (top_atmos - bottom)/nz`
(`vertical_grid`, `photochem/src/photochem_eqns.f90`), with two different
`top_atmos`. Nothing pins `top_atmos` at the exit: `robust_step`
(`photochem/photochem/extensions/gasgiants.py`) re-pins it only every
`freq_update_TOA = 1000` internal steps and then accepts any top pressure
within a **factor 3** of the requested one. The 2.1e-6 K decides which phase of
that cycle the integration exits on, and no more.

**Excluded by measurement.** (i) *An amplification of the temperature*: sorted
by the deep boundary temperature the seven runs interleave -- `g401` at
424.95340782 K carries the lower model top, `g400_1em4` at 424.95340736 K the
higher -- and a pair 3.0e-07 K apart within one group agrees to 4.2e-06 in O/H
while a pair 4.6e-07 K apart across the groups differs by 5.2e-03. (ii) *The
photochemical stopping tolerance*: four runs repeated with `conv_longdy`
1e-2 -> 1e-4, `conv_longdydt` 1e-6 -> 1e-8 and `equilibrium_time` 1e17 -> 1e22 s
give the same two groups and the same spreads (2.466e-02 in pressure,
5.199e-03 in O/H) and reproduce their loose-tolerance twins to 4e-07.
(iii) *A second chemical solution*: within a group the composition reproduces
to 1e-06 across guesses spanning 300-450 K and a factor 1e5 in stopping time.
(iv) *A cold-trap cell switching in the climate solve*: the climate grid is
fixed by its two stated pressures and 60 layers and does not move at all, and
the tropopause reproduces to 2.6e-08. (v) *Round-off amplification*: the code
is bit-deterministic and the outcome is two discrete states.

**The matching pressure never moved.** It is an exact node inserted by
`lower_profile_schema.insert_level` and reads `1.000000000000e-06` bar in every
profile written. The 2.4 % is a comparison at fixed level index; EXHALE reads
the file as a table over pressure.

**What the chain actually inherits.** At the matching level, over the seven
runs: T 2.7e-09, r 4.8e-09, `q_H2` 5.7e-08, **He/H 6.5e-08**, C/H 3.1e-07, N/H
3.2e-07, O/H 5.2e-03. `n_tot` and `rho` are carried and explicitly not imposed
(`input_read.f90`). The only place the 2.4 % reaches EXHALE at all is
`lap_r_top_RJ`, the upper edge of the elemental flux window
(`binary_element_diffusion.f90`), which moves by 5.2e-05 R_p. Against
`tol = 0.05` on `|F - Phi|/|F|`, and against the 0.084-0.085 flux-window spread
the two LHS 1140 b arms carry, the handoff's own reproducibility is six orders
below the tolerance in He/H and does not enter the H residual at all. Worth
stating: the closure re-runs the adapter in every iteration directory, so the
oxygen reservoir carries a 0.5 % step between iterations of one ladder --
nothing in the He I 10830 chain reads it, but a result that read the oxygen
would inherit it.

**The change.** `src/utils/lower_profile_schema.py`: `insert_level`
interpolated every column linearly in log p, including `n_tot` and `rho`, which
are not free -- `n = p/(k_B T)` holds at every node the solution was computed
on. Measured on these columns the identity held at 101 of 102 levels and was
violated at **exactly one, the matching level**, by 7.7e-04 and 1.24e-03. It
now recomputes `n_tot` from the inserted pressure and temperature and `rho`
from the interpolated mean molecular weight. Verified by rerunning one case:
exactly two cells of the file change, the header and `solution_id` are
unchanged, `n_tot/(p/k_B T) - 1` over the whole table goes from 7.706e-04 to
2.2e-16, and the run-to-run spread of `n_tot` at the matching level falls from
4.646e-04 to 4.682e-10. Also corrected: the module docstring of
`src/utils/element_flux_closure.py` stated the convergence test as
`|F - Phi|/max(|F|,|Phi|)` while the function implements -- and its own
docstring documents -- `|F - Phi|/|F|`.

**Not repaired, and proposed instead.** The grid state at the exit is
Photochem's. Re-pinning the grid and re-testing convergence before
`reached_steady_state` is set is the right repair and needs a rebuild of
`photochem/`; doing the same from the adapter after
`photochemical_steady_state` returns (`pc.update_vertical_grid`, then step to
convergence again) needs no rebuild but changes the converged answer of every
future handoff column, against which the stored `crossings_gm25` and
`crossings_pc090` ladders were measured. Recorded as a proposal, not made.

**Scope.** No EXHALE Fortran source or binary changed, no golden was refreshed
and `make check` was not run: the `lower_profile` regression case reads a
stored profile and never runs the adapter, so the schema change cannot reach
it. Not measured: another reservoir value, another planet, another `K_zz`, and
whether the two grid groups are the only two.

**Files.** `docs/lower_profile_deep_boundary_sensitivity.md` (the diagnosis);
`LHS1140b/exhale/deep_temperature_sensitivity/` (the record: README, six
scripts, three tables, thirteen adapter runs);
`src/utils/lower_profile_schema.py` (`insert_level`);
`src/utils/element_flux_closure.py` (the docstring). The observation this
answers: `LHS1140b/exhale/crossings_pc090/results.txt` section 7 and section 91
above.

## 93. Rows of one residual that differed by a million, and a third body corrected only halfway (2026-08-29)

**The judgement first.** Section 82 left two items open in the molecular
network: a NaN abort at He/H = 1 that it did not diagnose, and the suspicion
that the molecular `hybrd1` system fails at high helium because its rows are
not scaled. Both are settled here, and they are **not the same problem**.

* **The rows were the problem they were suspected of being, and scaling them
  fixes it.** Every row of the molecular system is now divided by its own
  turnover rate -- the rate at which the species that row balances is produced
  or destroyed in that cell -- so what `hybrd1` minimizes is the *relative*
  imbalance of each species balance instead of an absolute reaction rate. The
  fraction of solver attempts that end on MINPACK's `info = 4` ("iteration is
  not making good progress") falls from **54.6 % to 0.01 %** at He/H = 1000,
  from 41.3 % to 0.01 % at 100 and from 36.3 % to 0.26 % at 10.
* **The third-body fix of section 82 had reached only one of its three
  reactions.** `ionization_equilibrium` set `ieq_cell%ntot` from `calc_ntot`,
  which is what the R12 term of the residual reads, but still called
  `set_mol_coeffs(T, n_in_dim)` -- and R13 and R15 fold the third body into
  *their coefficient*, so both were still evaluated with `rho/m_H`. Corrected;
  the goldens move a second time, for the same reason they moved the first.
* **The He/H = 1 NaN is in the hydrodynamics, not in the network**, and it is
  now located exactly: face `j = 232`, `r = 1.1398 R_p`, where the reconstructed
  right state carries `p = -3.65e-2` in code units and the first invalid
  operation of the whole run is `aR = sqrt(g*pR/rhoR)`
  (`src/modules/flux/Num_Fluxes.f90:46`). Row scaling does not move it by one
  step. It is recorded as an open item, not patched, because the patch belongs
  outside the molecular branch.

### The measurement

`backup/regression/mol_diffusion` (hot Uranus, molecules + He 2^3S + binary
H/He diffusion + a `base.inp` handoff), copied to scratch with **only**
`He/H number ratio` in `input.inp` and `HeH_base` in `base.inp` changed
together, `OMP_NUM_THREADS=1`, `EXHALE_MAXSTEPS=12000`. Three binaries: the
tree as found, the tree with the third body corrected, and the tree with the
third body corrected **and** the rows scaled. The `info` histogram is a new
run-wide diagnostic (`ieq_n_mol_info`), reported next to the clamped-cell
count; "fail" below is `info = 4` as a fraction of all molecular `hybrd1`
attempts.

| build | He/H | steps | NaN abort | clamped cells | hybrd1 attempts `info=1` | `info=4` | fail | log10 Mdot |
|---|---|---|---|---|---|---|---|---|
| as found | 1    | 12000 | no  | 614 | 6031444 | 56627   | 0.93 % | 11.06 |
| as found | 10   | 12000 | no  | 166 | 5405217 | 2169365 | 28.6 % | 11.04 |
| as found | 100  | 12000 | no  | 212 | 4972511 | 3499546 | 41.3 % | 11.27 |
| as found | 1000 | 12000 | no  | 16  | 4360657 | 5242809 | 54.6 % | 11.11 |
| + third body | 1    | 3406  | yes | 507 | 1677801 | 148735  | 8.1 %  | 8.91 |
| + third body | 10   | 8827  | yes | 274 | 3754138 | 2139333 | 36.3 % | 8.85 |
| + third body | 100  | 12000 | no  | 216 | 5051123 | 3285000 | 39.4 % | 11.26 |
| + third body | 1000 | 12000 | no  | 7   | 4368719 | 5222557 | 54.5 % | 11.11 |
| + row scaling | 1    | 3406  | yes | 337 | 1702668 | 48098 | 2.8 %  | 8.91 |
| + row scaling | 10   | 5474  | yes | 259 | 2758854 | 7098  | 0.26 % | 8.87 |
| + row scaling | 100  | 12000 | no  | 408 | 6048842 | 670   | 0.01 % | 11.26 |
| + row scaling | 1000 | 12000 | no  | 38  | 6048079 | 450   | 0.01 % | 11.09 |

These are relaxation snapshots followed by the steady solver, not a
composition study; the `Mdot` column is a diagnostic of where each run went.
What it does show is that the row scaling is a change of *path*: at He/H = 100
and 1000, where every arm completes, `Mdot` moves by 0.01 and 0.02 dex while
the failed-attempt fraction falls by three and a half orders of magnitude.
A positive diagonal scaling of a residual has exactly the same zeros, so this
is what should happen.

**The clamped-cell count does not follow the `info` count, and that is worth
saying plainly.** At He/H = 100 it rises, 216 to 408, and at 1000, 7 to 38.
The clamp counts cells where *no* starting point produced an admissible root;
with the rows scaled, `hybrd1` now reaches its tolerance on nearly every
attempt, and a converged root that lies on a face of the simplex (a fully
dissociated or fully ionized species) is rejected by the admissibility test
where a non-converged root that happened to sit inside it was accepted before.
The counts are rare either way -- 408 clamps against 6.05 million cell solves,
7e-5 -- and the solution does not move with them. The suspicion recorded in
section 82 was that row scaling drives the simplex failures; the measurement
says it drives the **solver** failures, by three and a half orders of
magnitude, and that the residual clamps are a separate, much smaller
population.

### What the scale is

A row of this system is a production-loss balance in cm^-3 s^-1, and its
magnitude is set by the element that carries it: the helium rows run as
`n_He` times a rate, the hydrogen and molecular rows as `n_H` times a rate,
with a second density inside every bilinear term. At He/H = 1000 the two
blocks of one residual vector differ by about 10^6. MINPACK's `hybrd1` scales
the *variables* (its internal `diag`) and never the rows: the dogleg step
minimizes `||J dx + f||_2`, in which a block 10^6 below the other carries no
weight, so the unknowns that block determines -- the molecular fractions --
are left wherever the helium block puts them.

The scale each row is divided by is **the sum of the magnitudes its terms
reach when every species is set to the whole of its element**: each rate
coefficient of the row times the cell's `n_H`, `n_He`, `n_e` and `n_tot`. It
is an upper bound on the row, strictly positive wherever the row carries any
reaction, and -- this is the part that matters for the solver -- **constant
across the cell's solve**, so the scaled residual is the same smooth function
of `x` that the finite-difference Jacobian of `hybrd1` assumes. The electron
density used is the cell's incoming `n_e` and not the neutrality bound
`n_H + 2 n_He`: in the shielded molecular base the bound overstates the
electron terms by many orders and would scale those rows into insignificance,
which is the failure this is meant to remove.

`set_mol_turnover_rates` (`System_HeH_mol`) fills rows 1-8 and
`set_mol_metal_turnover_rates` (`System_HeH_mol_metals`) the metal rows above
them, so the whole system reaches the solver equilibrated; a metal row runs as
`n_Xe` times a rate and a trace element carries `n_Xe ~ 1e-4 n_H`, so it sits
below the helium block for the same reason the molecular one does. The
identity rows `metal_rows` writes for an absent element, and for the X++ of a
two-stage element, are already dimensionless and keep a scale of 1. The
routine mirrors `mol_heh_rows` term by term and says so at the code site: a
reaction added to a row has to be added to its scale as well, and a missing
term only makes the scale a weaker bound -- it cannot make the root wrong.

### The third body, the other two reactions

Section 82's account is right about what the third body is and wrong about how
far the fix reached. R12 takes `n_M` in the residual, through
`ieq_cell%ntot`, and that path was corrected. R13 and R15 are written in
`mol_rates` as *two-body-equivalent* coefficients, `k = 3.2e-29 n` and
`k = 8e-33 (300/T)^0.6 n`, so the third body enters when the coefficient is
built -- in `set_mol_coeffs`, which `ionization_equilibrium` still called with
`n_in_dim = rho/m_H`. Both now get `n_tot(j)`, the same particle density the
R12 term reads. R13 is the H3+ source at the base and R15 the H2 source, so
this is not a small correction to a small term: it is the pair that builds the
molecular base, and it was running with a third body 2.3x too dense at the
hot-Uranus base and up to 4x in a helium-dominated one.

**The corrected third body is what brings the NaN abort back.** The tree as
found does not abort at any He/H tested -- section 82's abort at step 5538 is
absent from it -- and with R13/R15 corrected the He/H = 1 arm aborts at step
3406 and the He/H = 10 arm at 8827. The physically correct third body is kept;
the abort it exposes is diagnosed below and recorded as an open item.

### Where the He/H = 1 NaN is

Reproduced in a build with `-fcheck=bounds,do,mem` and
`-ffpe-trap=invalid,zero,overflow` (the detector `backup/regression/run_fcheck.sh`
uses, built to `EXHALE_fcheck.x` so it does not overwrite the production
binary). The first invalid operation of the run is at step 3406:

```
Program received signal SIGFPE
#3  __numerical_fluxes_MOD_num_flux  at src/modules/flux/Num_Fluxes.f90:46
#4  __rk_integration_MOD_rk_rhs      at src/modules/time_step/RK_rhs.f90:62
```

Line 46 is `aR = sqrt(g*pR/rhoR)`. A second build, instrumented to print the
face states, names the face and the state that reaches it:

```
face j = 232   r_edg = 1.13983 R_p
WL (rho,v,p) =  2.9894   26.707   +0.23329
WR (rho,v,p) =  2.8916   26.816   -0.036499
cell 231 (rho, rho v, E) = 2.7948   75.246   1013.68
cell 232                 = 2.8916   77.540   1040.21
cell 233                 = 5.8744  152.863   1988.84
```

The reconstructed right state has a **negative pressure**, and HLLC takes the
square root of it. The layer it sits in explains how: below `r = 1.139` the
run has developed a cold hypersonic shell -- `T` about 650 K, `v` about
1e7 cm/s, **Mach 50 to 60** -- riding on a dense hot wall at `r = 1.144`
(`rho` up by 20x, `T` 9100 K). At Mach 60 the thermal pressure is 2e-4 of the
total energy the scheme carries, so the pressure recovered as
`(g-1)(E - rho v^2/2)` is the difference of two numbers that agree to four
digits, and a reconstruction across the jump at the next cell tips it below
zero. Nothing between the reconstruction and the Riemann solver tests the
reconstructed pressure for positivity.

This is a hydrodynamic robustness gap, not a molecular one: it is identical,
to the step and to the face, with and without the row scaling, and the
positivity floor that would close it belongs in the reconstruction, outside
the `thereis_mol` branch this series is confined to. Recorded as
`TO_BE_DONE.md` item (O).

### Containment, and the goldens

Every change is inside `thereis_mol` -- the two `set_mol_*_turnover_rates`
calls and the corrected `set_mol_coeffs` argument sit in the
`if (thereis_mol)` cell-state block of `ioniz_eq`, the scaling is applied only
in the two molecular residual routines, and the `info` histogram is counted
only in the molecular attempt loop. Two measurements say so rather than the
code reading alone:

* `make check` passes byte-identical on `wasp_full`, `wasp_he23off` and
  `lower_profile` -- the three non-molecular cases -- against goldens that
  predate the series.
* The LHS 1140 b arm `LHS1140b/exhale/crossings_gm25/rs_heh0p425`, re-solved in
  a scratch copy from its stored restart with the pre-series binary and with
  the new one, writes **byte-identical** `Hydro_ioniz.txt`, `Ion_species.txt`,
  both `_adv` files, both `_IC` files, `Cooling_breakdown.txt` and
  `Heating_breakdown.txt`. The stored run directory was not touched.

**The four molecular cases move, and the third body is what moves them.** Each
was also re-run with the third body corrected but the rows *unscaled*, so the
two changes can be told apart. Largest relative difference over the column
against the previous golden, with the median in brackets (these are 12000-step
relaxation snapshots, so the maximum sits in the base sound-wave layer or in
the thin outer cells and the median is the better measure of the profile):

| case | rho | p | T | cool | H2 front (f = 0.5) | log10 Mdot |
|---|---|---|---|---|---|---|
| `mol_base_handoff` | 9.1% (0.20%) | 14.9% (0.32%) | 6.1% (0.087%) | 53.5% (1.8%) | 1.1617 -> 1.1577 | 10.57 -> 10.58 |
| `mol_metals` | 6.5% (0.16%) | 9.8% (0.34%) | 6.2% (0.093%) | 68.3% (0.67%) | 1.1637 -> 1.1597 | 10.58 -> 10.58 |
| `mol_lyman_werner` | 27.6% (0.57%) | 45.3% (0.71%) | 17.4% (0.32%) | 59.8% (2.4%) | 1.1321 -> 1.1271 | 10.58 -> 10.49 |
| `mol_diffusion` | 9.1% (0.20%) | 14.9% (0.32%) | 6.1% (0.087%) | 53.5% (1.8%) | 1.1617 -> 1.1577 | 10.57 -> 10.58 |

The row scaling contributes almost none of that. Measured on its own -- the
third-body-corrected build against the third-body-corrected **and** scaled
build, same case, same steps:

| case | rho | p | T | H2 front |
|---|---|---|---|---|
| `mol_base_handoff` | 0.019% | 0.0066% | 0.018% | unchanged |
| `mol_metals` | 0.28% | 0.46% | 0.23% | unchanged |
| `mol_lyman_werner` | 0.013% | 0.020% | 0.0079% | unchanged |
| `mol_diffusion` | 0.017% | 0.0068% | 0.016% | unchanged |

That is the judgement the goldens were refreshed on. The move is the physics
correction, and it is the correction the section-82 audit intended and
described: R13 builds H3+ at the base and R15 builds H2, both linearly in the
third-body density, which was 2.3x too large at this hot-Uranus base. The H2
front shifts by 0.004 R_p and the mass-loss rate by 0.01 dex on three of the
four cases and by 0.09 dex on the Lyman-Werner one, where the photodissociation
front and the three-body formation front sit on top of each other and the base
is therefore the most sensitive. The row scaling, in contrast, leaves the
solution where it was to a few parts in 10^4 while removing the solver failures
-- which is what a diagonal rescaling of a residual should do, and is the
evidence that it is a change of path and not of physics.

Goldens re-snapshotted for the four molecular cases at the end of the series
(`run_check.sh golden mol_base_handoff mol_metals mol_lyman_werner
mol_diffusion`) and `make check` re-run afterwards.

### Files changed

| File | Change |
|---|---|
| `src/modules/nonlinear_system_solver/System_HeH_mol.f90` | `mol_inv_turnover` and `set_mol_turnover_rates`; the scaling applied at the end of `ion_system_HeH_mol` |
| `src/modules/nonlinear_system_solver/System_HeH_mol_metals.f90` | `set_mol_metal_turnover_rates` for the metal block; the scaling applied at the end of `ion_system_HeH_mol_metals` |
| `src/modules/radiation/ionization_equilibrium.f90` | `set_mol_coeffs` given `n_tot(j)` instead of `n_in_dim(j)`; the two turnover-scale calls; `ieq_n_mol_info`, the run-wide `hybrd1` exit-code histogram of the molecular solves |
| `src/EXHALE_main.f90` | the `info` histogram reported at the end of a run, next to the clamped-cell count |
| `TO_BE_DONE.md` | new item (O), the negative reconstructed pressure |
| `docs/molecular_chemistry_audit_he_rich.md` | section 7: both open items of its section 6 closed |

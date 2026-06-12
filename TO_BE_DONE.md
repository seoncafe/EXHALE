# TO BE DONE

Running list of known limitations and planned improvements.

## DONE 2026-06-13 — Metals in the mass/charge budget (`eos_metals`, default ON)

*(Original entry 2026-06-12, implemented 2026-06-13. Logged as
`docs/Update_ATES` §17; the same-day `IC mode: auto` selector is §18.)*

Metals used to be strictly passive trace species: they affected
heating/cooling and the ionization equilibrium, but were excluded from the
mass density, the base density, and the electron/particle densities used
outside the ionization solver. (The only self-consistent place was the
charge balance **inside** the MINPACK systems, `System_HeH_metals.f90`.)

**Implemented** (runtime key `eos_metals 0|1` in `metals.inp`; `1` = new
default, `0` = legacy trace approximation; metals-off runs identical
either way):

- **Base density / composition constants** — `rho_bc` now includes the
  metal mass: `mass_per_H = 1 + 4*HeH + sum(melem_ab*melem_A)` and
  `rho_bc = mass_per_H/(1 + HeH)` (`input_read.f90`); atomic weights in
  `species_table.f90:melem_A`. `n0` keeps its H/He-nuclei meaning, so
  `n_H = n0/(1+HeH)` is unchanged. (Correction to the original entry:
  `v0`, `p0`, `b0` are m_H-based normalizations and never contained the
  mean molecular weight — they were never affected.)
- **Mass density** — `calc_rho` takes an optional `nm` argument and adds
  `melem_A(elem)*nm` (all stages); passed from `ionization_equilibrium`.
- **Electron density** — `calc_ne` takes an optional `nm` and adds
  `stage*nm`; passed everywhere nm is available: `composition`
  (main-loop EOS), `ionization_equilibrium` (opacity/cooling),
  `util_ion_eq` (eval_cool + cooling dump), `energy_semi_implicit`,
  `excited_hydrogen`, `post_process_adv`. The scalar `T_equation` already
  carried metal electrons via `pp_nm_cell`.
- **Particle density** — `calc_ntot` likewise adds the metal nuclei
  (EOS `T = p/(n_tot + n_e)`).
- **Ghost pressure / dp_bc** — the lower-BC ghost pressure is now
  `ntot_bc + dp_bc` with `ntot_bc = (1 + HeH + sum(melem_ab))/(1 + HeH)`
  (so T = T0 holds exactly at the base), and `dp_bc` includes the metal
  electrons. `set_IC`/`load_IC` use `mass_per_H`/`ntot_bc` consistently.

**Validation (WASP-121b regression matrix, 2026-06-13):**
- `eos_metals 0` is **byte-identical** to the pre-change goldens for both
  cases (wasp_full count=7296, wasp_he23off count=7300) — the legacy path
  is untouched.
- `eos_metals 1` (new default): converges cleanly (count 7145/7142, ~2%
  fewer steps); profile shifts are at the expected ~1% mass scale —
  median |dn|~2.4%, |dp|~2.3%, |dT|~0.9% (max |dT| = 132 K at the
  breathing-base region r~1.09, peak rel. 3.8% at r~1.007);
  number flux n*v*r^2 -1.1%, i.e. log10 Mdot 13.31 -> 13.30 g/s — the
  -1.1% number flux offsets the +1.2% mass per particle, so the MASS
  loss rate is essentially unchanged.
- Goldens re-snapshotted under the new default; a repeat run reproduces
  them byte-identically (single-thread determinism intact).

**Still open (small):** the per-cell scalar Brent/Newton `T_equation`
energy balance uses `pp_nm_cell` only in the post-process path; the
time-marching per-cell solve receives ne via `energy_semi_implicit`
(now metal-aware). The `use_2lev_cool` legacy branch remains
AIOLOS-hard-coded (default off).

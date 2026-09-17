The lower atmosphere handed over as a **profile** instead of the six scalars of `base.inp` (`Lower atmosphere profile:`, `docs/phase_e_flux_closure_design.md`, `docs/Update_EXHALE.md` sections 76-79). Copied from `examples/17_lower_profile`: HD 209458 b, `Include He23S? True`, `He_diffusion: True`, `He_metal_diffusion: True`, `Solver: Newton`, with a provenance-only `base.inp` added beside the profile and a 12000-step cap, the same relaxation-snapshot convention `mol_diffusion` uses.

What the case pins, all of it in one run:

- the profile reader and its schema check (`src/modules/files_IO/lower_atmosphere_profile.f90`);
- the matching-level base state taken from the table at `p_match_bar = 1e-6 bar`: `T0 = 1450.0 K`, `R0 = 1.401 R_J`, `p_base`, `q_H2(base) = 0.11961`, `He/H = 0.083333`;
- **the base level itself**: the wind starts at the matching level, so `n0` follows from it, `n0 = p_match/(k_B T0 ntot_bc) = 4.9913e12 cm^-3` with `ntot_bc = 1.000764`. `input.inp` carries no `Log10 lower boundary number density` key: a density key beside the profile would state the level a second time and is refused unless the two agree within 1%. Until 2026-09-05 the case ran at the density key's `n0 = 1e14 cm^-3`, i.e. `p = 2.0035e-5 bar`, 20.03 times deeper than the level whose `T`, `r`, `q_H2` and elemental ratios it was using; the golden of this case therefore moves entirely (log10 Mdot 8.44 -> 7.92 g/s at 12000 steps);
- the elemental reservoirs the profile carries (`C/H = 2.700e-4`, `O/H = 4.900e-4`, `N/H = 6.800e-5`) reaching `melem_ab` through `set_element_abundance`. The case is **metals-on with no `metals.inp`**: the profile is the only source of the metal abundances, so a regression in that route shows up here and nowhere else;
- `K_zz(p)` interpolated onto the grid instead of the scalar `He_Kzz`, which is inert here. The column is not constant (`1.130e9` to `5.623e9 cm^2/s` on this grid), so the interpolation is exercised rather than degenerate;
- the **accepting** branch of the profile / `base.inp` pair enforcement: `base.inp` carries a `# solution_id` comment and no physics key at all, and the run prints `base.inp: solution_id matches the profile.` A profile or a `base.inp` edited without the other stops the run;
- the elemental flux window statistics (`reduce_element_flux_windows`) and `output/element_flux_profile.txt`, which are written unconditionally when a profile is in use.

Molecular chemistry is **off** here: `mol_diffusion` already guards the molecular carriers of the diffusion operator, and this case is about the handoff route, on an atomic base where the profile's own `q_H2` is what sets the base particle count.

Only `Hydro_ioniz.txt` and `Ion_species.txt` are compared bitwise, as in every case; `EXHALE_resolved.out` and `element_flux_profile.txt` are produced but not compared.

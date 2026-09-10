# Ready-made example configurations (HD189733b)

One folder per solver/physics combination. Folders `01`-`11` are the
same planet (HD189733b) so the effect of each option can be isolated;
`12_windae_ic_hd209`, `14_diffusion`, `15_molecular`, `16_molecular_metals`
`17_lower_profile`, `18_oxygen_chemistry` and `19_molecular_ir_bands` use
HD209458b, where
those lower-atmosphere / diffusion options are validated (`11` uses HD189733b
far from the seed to exercise the self-consistent-BC continuation). Each
folder is self-contained: `cd` into it and run the repo-root binary,

```sh
cd examples/03_newton
OMP_NUM_THREADS=8 ../../EXHALE.x
```

Outputs land in the folder's own `output/`. In the table below the **Added /
changed lines** are relative to each example's natural base, not always `01`:
the legacy `01_legacy_marching/input.inp` for the solver examples `02`-`04`
(and the Wind-AE `11`/`12`), and the recommended-solver `03_newton/input.inp`
(two-stage `PLM+WENO3` + `Solver: Newton`) for the physics examples `05`-`10`,
which all inherit that solver and add only their signature option on top. So
`diff 03_newton/input.inp 05_metals/input.inp` shows exactly what the metals
add, while `diff 01_legacy_marching/input.inp 03_newton/input.inp` shows the
solver progression. See the user manual (`docs/EXHALE_user_manual.pdf`)
Sect. "Worked examples" for the corresponding table, and Sect. 2 for what
each option does.

The metastable helium triplet is on by default in the code, but the ladder
keeps `Include He23S? False` on purpose: in both bases (`01_legacy_marching`,
`03_newton`) and in every rung whose signature option is something else
(`02`, `04`, `05`, `09`-`12`). Atomic helium is the ladder's baseline, which
is what makes `06_he23s` a one-line diff against `03_newton` instead of a
configuration that differs in nothing. Everywhere outside the ladder --
`tutorial/`, `13`-`17`, the planet folders -- the triplet is on, and an
input file that simply omits the line gets it.

| Folder | Demonstrates | Added / changed lines |
|---|---|---|
| `01_legacy_marching` | Classic ATES v2 baseline: PLM marching, H/He only, Roche domain, `du`-based stop | (none, baseline) |
| `02_two_stage` | Two-stage PLM -> WENO3 marching | `Reconstruction scheme: PLM+WENO3`, `du_th [PLM,WENO3]: 0.5 1.0e-3` |
| `03_newton` | **Recommended default**: two-stage warm-up + JFNK Newton finish to the true steady state | + `Solver: Newton` |
| `04_newton_from_state` | Newton re-convergence of an existing state (no marching) | `Load IC? True`, `Resid tol: 1.0e-3`; run with `EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 ../../EXHALE.x` |
| `05_metals` | Trace-metal cooling (solar C/N/O/Mg/Ca/Na/Fe) | `metals.inp` present |
| `06_he23s` | He 2^3S metastable level (He I 10830 line): the one rung that leaves the ladder's atomic-helium baseline | `Include He23S? True` |
| `07_balmer_lya` | Non-LTE H(n=2) + Ly-alpha pumping (H-alpha/H-beta) | + `Stellar Teff/radius`, `Deexc heat`, `Jlya escape-prob`, `Stellar Lya flux/halfwidth/boost` |
| `08_full` | Everything on (= the `HD189733b/` planet folder + Newton); feeds `EXHALE_transit.py` | 07 + `metals.inp` |
| `09_spherical` | Spherical domain instead of the default Roche/L1 truncation | `Domain mode: Spherical`, `Outer radius [R_p]: 10.0` |
| `10_warm_seed_ic` | Warm-seed (hot Parker overlay) initial condition | `Hot Parker IC: 10000` |
| `11_windae_ic` | In-process Wind-AE IC for HD189733b: far from the seed, so the **self-consistent-BC continuation** (base-BC re-convergence + molecular-layer turn-off) is exercised; the Wind-AE ramp converges and writes the IC. EXHALE's *own* HD189733b base-breathing instability (separate from the IC) then limits the warm start | `+ IC mode: windae`, `+ Solver: Newton` |
| `12_windae_ic_hd209` | In-process Wind-AE warm-start IC that **works**: HD209458b (not HD189733b), close to the shipped seed, so the ramp converges and EXHALE warm-starts cleanly (spherical 10 Rp) | HD209458b params `+ Domain mode: Spherical`, `IC mode: windae`, `Solver: Newton` |
| `13_lower_atmosphere` | **Lower-atmosphere connection** for four planets: an analytic 1-ubar base column plus a `base.inp` handoff (isothermal or Guillot T(p) generator), see the folder's own README (multi-planet, not the HD189733b baseline) | driver-generated `base_iso.inp` / `base_guillot.inp` (the four `input.inp` files are plain core blocks, all four with the triplet on; add `Lower column: <R_1bar>` by hand for the analytic-column report) |
| `14_diffusion` | **Diffusive separation of He and metals** (HD209458b): the He/H ratio declines with altitude and each trace metal settles independently, reshaping the He 10830 line | HD209458b params `+ Include He23S? True`, `+ He_diffusion: True`, `+ He_metal_diffusion: True`, `+ He_Kzz: 1.0e9`, `+ He_alphaT: 0.0`, `metals.inp` present |
| `15_molecular` | **Full molecular chemistry** (HD209458b): H2/H2+/H3+/HeH+ in the coupled ionization equilibrium; a sharp H2->H front forms above a thin molecular base, the wind above it essentially atomic. Converges to a residual norm of 6.105e-6 in 89 outer iterations from a cold start with the current code (re-run 2026-08-15; the front sits at r = 1.0097 R_p and the molecular layer between the base and the front collapses to ~400 K), see the note below | HD209458b params `+ Molecular chemistry: True`, `+ Solver: Newton 5.0e-2`, `+ Resid tol: 1.0e-5`, `+ Max steps: 150000` (no `metals.inp` here, but metals are allowed, see `16_molecular_metals`; `He_diffusion` is allowed too, the element transport closes over the molecular carriers) |
| `16_molecular_metals` | **Molecular chemistry + trace metals in one system** (HD209458b): the H2/H2+/H3+/HeH+ network and the metal ionization stages share the free electron density, which the metals dominate in the shielded molecular base. Converges only from a warm restart off a freshly converged `15`, with `Low-Mach damping` on, the recipe is in the note below | 15 `+ metals.inp` (solar C/N/O), without 15's three convergence keys |
| `18_oxygen_chemistry` | **Oxygen chemistry** (HD209458b): OH, H2O and CO added to the same coupled system as the molecular network, with the FUV photolysis of H2O and OH and the molecular carriers transported, so the base H2/H partition is computed rather than imported. The four band fluxes are the HD 209458 b flux at the planet integrated over the band edges; the first band is the 912-1201 A Lyman-Werner interval, so `Stellar LW flux: 480.9` (the same value the Lyman-Werner cases use, and the sum 343.0 + 137.9 of the two halves the band was split into until 2026-09-06) is the band flux there. Note that HD 209458 b's 2331 K base is **not** the regime the option exists for: there `H2 + M -> H + H + M` runs the partition and the oxygen cycle adds a few percent. It is the example because it is a molecular configuration that converges; the cool-base case the option is for is HD 189733 b, whose band fluxes are in `README_HOWTO.md` | 16 `+ Oxygen chemistry: True`, `+ Stellar FUV B3/B4 flux`, `+ Stellar Lya flux`, `+ Stellar LW flux` |
| `19_molecular_ir_bands` | **Molecular infrared bands** (HD209458b): `18` plus the coolants an H2 atmosphere carries below the H2 -> H front and the atomic channels do not -- the H2 quadrupole and magnetic dipole line spectrum (Roueff et al. 2019) and the H2O and CO vibration-rotation bands (HITEMP). Each is the NET exchange with the diluted `B_nu(T0)` that `Base IR field` supplies, so each stops cooling at its own radiative equilibrium temperature instead of running the layer down; the two keys belong together and `input_read` warns if only one is set. Writes the `H2_IR`, `H2O_IR` and `CO_IR` net columns and the Planck-mean optical depths to `output/Cooling_breakdown.txt`, and a column-integrated infrared block to `output/FUV_bands.txt`. TO_BE_DONE.md item (G); `docs/lower_atmosphere_coupling.pdf` section 10 | 18 `+ Base IR field: True`, `+ Molecular IR bands: True` |
| `17_lower_profile` | **Lower-atmosphere profile handoff** (HD209458b): the lower atmosphere handed over as a table over an interval of pressure instead of the scalars of `base.inp`. The base state (`T0`, `R0`, `p_base`, `q_H2`), the elemental reservoirs (He/H and solar C/N/O) and `K_zz(r)` all come from the one file; `make_example_profile.py` regenerates it. See `docs/input_schema.md` section 2d | 14 `+ Lower atmosphere profile: lower_atmosphere_profile.dat`, `- He_Kzz` (the profile carries K_zz), `- Log10 lower boundary number density` (the profile's matching level fixes the base level, and a density key beside it is refused unless it agrees within 1%), no `metals.inp` (the profile carries C/N/O) |

Notes
- Wind-AE IC (`11`/`12`, `docs/wind_ae_solver.pdf`): `IC mode: windae`
  builds the IC in-process from a shipped seed (`inputdata/windae_seed.csv`;
  folders `11` and `12` carry an `inputdata` symlink to the repo-root
  directory), via the self-consistent-BC continuation.
  This converges seed-adjacent hot Jupiters (`12`, HD209458b) and
  far-from-seed planets (`11`, HD189733b) alike. Whether EXHALE then
  time-integrates the result cleanly is a separate matter: HD189733b hits an
  EXHALE-side base-breathing instability regardless of the IC source.
- `04_newton_from_state` needs a state to start from: copy a converged
  `Hydro_ioniz.txt` / `Ion_species.txt` (e.g. from `03_newton/output/`) to
  `output/Hydro_ioniz_IC.txt` / `output/Ion_species_IC.txt` first.
- Metals on/off is a runtime switch: any folder becomes metals-on by
  copying a `metals.inp` into it (and metals-off by removing it).
- HD189733b has a "breathing" (slightly inflowing) base. The
  `||R|| ~ 3e-4` Newton floor it used to show (`docs/steady_solver_memo.pdf`)
  is gone since the 2026-08-11 line-search and base-energy fixes: a warm
  re-convergence now reaches `||R|| = 1.2e-5` in 11 iterations (re-measured
  2026-08-15, `backup/regression/jfnk_hd189_tight`).
- A molecular run (`15`, `16`) needs three keys the atomic examples do not.
  `Solver: Newton 5.0e-2` raises the hand-off threshold, because the `du`
  descent of a molecular run is not monotonic and the run can spend its whole
  step budget above the `1e-2` default. `Resid tol: 1.0e-5` is what actually
  converges the molecular layer: `||R||` is set by the two or three cells just
  above the base, so a run stopped at the `1e-3` default leaves the layer still
  cooling. `Max steps: 150000` covers the marching warm-up. Marching alone
  never reaches the layer (the H3+ cooling time there is ~1e7 CFL steps), so
  the Newton finish is not optional here. Details and the converged numbers:
  `docs/lower_atmosphere_coupling.md`, section "Converged Tier-2 solution".
  The three keys converge `15` (metals off) from a cold start. The
  `15_molecular/output/` stored here is the current-binary converged state
  (replaced 2026-08-15; the earlier binary's state is kept as
  `15_molecular/output_pre_item3_20260813/`). After any further code change,
  re-converge `15` before using it as a seed, see the next note.
- **Converging `16` (molecular + metals) needs a warm restart, and the seed
  matters.** The three keys alone do not converge it: with metals on, the cold
  run reaches the hand-off and the JFNK line search bottoms around
  `||R|| ~ 1e-3` (worst cell r = 1.021, momentum) and falls back to marching,
  at every `Low-Mach damping` coefficient tried, `off`, `5e-3`, `1e-2`,
  `2e-2` (measured 2026-08-15). The cold path for this case is closed. The
  recipe that does work:

  1. Converge `15` **with the binary you are about to use**. A state converged
     by an older binary is a different seed and it fails (`info = 2`, measured).
  2. Copy that run's `output/Hydro_ioniz.txt` and `output/Ion_species.txt` into
     the `16` run directory as `output/Hydro_ioniz_IC.txt` and
     `output/Ion_species_IC.txt`.
  3. Add `Load IC? True` and `Low-Mach damping: 5.0e-3` to `16`'s `input.inp`,
     alongside `15`'s three keys. Run the pair `5.0e-3` **and** `1.0e-2` and
     compare: a single coefficient is not trustworthy here (caveat below).
  4. Expect `info = 0` at `Resid tol: 1.0e-5` and `log10 Mdot` 9.65-9.67.

  Caveat: the two coefficients give two solutions that agree in the wind
  (`log10 Mdot` within 0.02, `rho` within 5.1% above 1.3 R_p) and differ inside
  the stagnant layer below the escape radius by tens of percent in `rho` and
  `T`, because the convergence measure is taken over `r >= r_esc` and does not
  look there. That propagates to the Balmer lines at the 5-8% relative level
  (H-alpha peak 2.098% against 2.199%, H-beta 1.046% against 1.134%); He I
  10830 is insensitive (40.19% against 40.26%). Full account:
  `docs/hd209_metal_stagnation.md` section 9, `docs/Update_EXHALE_stage1.md`
  section 61. The hot-Uranus molecular+metals case converges cold, so this is
  specific to HD 209458 b.

## 13_lower_atmosphere/
Lower-atmosphere connection examples for HD 209458 b, HD 189733 b,
WASP-121 b and WASP-52 b: one `input.inp` per planet + driver-generated
`base_iso.inp` / `base_guillot.inp` handoff files, with a results table and
regeneration commands in its own README.  Full description:
`docs/lower_atmosphere_coupling.pdf`.

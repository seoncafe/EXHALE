# Ready-made example configurations (HD189733b)

One folder per solver/physics combination. Folders `01`–`11` are the
same planet (HD189733b) so the effect of each option can be isolated;
`12_windae_ic_hd209`, `14_diffusion` and `15_molecular` use HD209458b, where
those lower-atmosphere / diffusion options are validated (`11` uses HD189733b
far from the seed to exercise the self-consistent-BC continuation). Each
folder is self-contained: `cd` into it and run the repo-root binary,

```sh
cd examples/03_newton
OMP_NUM_THREADS=8 ../../EXHALE.x
```

Outputs land in the folder's own `output/`. The folders differ from the
baseline `01_legacy_marching/input.inp` only by the lines quoted below
(`diff 01_legacy_marching/input.inp <other>/input.inp` shows exactly what
each option adds). See the user manual (`docs/EXHALE_user_manual.pdf`)
Sect. "Worked examples" for the corresponding table, and Sect. 2 for what
each option does.

| Folder | Demonstrates | Added / changed lines |
|---|---|---|
| `01_legacy_marching` | Classic ATES v2 baseline: PLM marching, H/He only, Roche domain, `du`-based stop | (none — baseline) |
| `02_two_stage` | Two-stage PLM -> WENO3 marching | `Reconstruction scheme: PLM+WENO3`, `du_th [PLM,WENO3]: 0.5 1.0e-3` |
| `03_newton` | **Recommended default**: two-stage warm-up + JFNK Newton finish to the true steady state | + `Solver: Newton` |
| `04_newton_from_state` | Newton re-convergence of an existing state (no marching) | `Load IC? True`, `Valve eps: 1.0e-4`, `Resid tol: 1.0e-3`; run with `EXHALE_PTC=1 EXHALE_PTC_JFNK=1 EXHALE_PTC_DTAU0=1.0 ../../EXHALE.x` |
| `05_metals` | Trace-metal cooling (solar C/N/O/Mg/Ca/Na/Fe) | `metals.inp` present |
| `06_he23s` | He 2^3S metastable level (He I 10830 line) | `Include He23S? True` |
| `07_balmer_lya` | Non-LTE H(n=2) + Ly-alpha pumping (H-alpha/H-beta) | + `Stellar Teff/radius`, `Deexc heat`, `Jlya escape-prob`, `Stellar Lya flux/halfwidth/boost` |
| `08_full` | Everything on (= the `HD189733b/` planet folder + Newton); feeds `EXHALE_transit.py` | 07 + `metals.inp` |
| `09_spherical` | Spherical domain instead of the default Roche/L1 truncation | `Domain mode: Spherical`, `Outer radius [R_p]: 10.0` |
| `10_warm_seed_ic` | Warm-seed (hot Parker overlay) initial condition | `Hot Parker IC: 10000` |
| `11_windae_ic` | In-process Wind-AE IC for HD189733b — far from the seed, so the **self-consistent-BC continuation** (base-BC re-convergence + molecular-layer turn-off) is exercised; the Wind-AE ramp converges and writes the IC. EXHALE's *own* HD189733b base-breathing instability (separate from the IC) then limits the warm start | `+ IC mode: windae`, `+ Solver: Newton` |
| `12_windae_ic_hd209` | In-process Wind-AE warm-start IC that **works** — HD209458b (not HD189733b), close to the shipped seed, so the ramp converges and EXHALE warm-starts cleanly (spherical 10 Rp) | HD209458b params `+ Domain mode: Spherical`, `IC mode: windae`, `Solver: Newton` |
| `13_lower_atmosphere` | **Lower-atmosphere connection** for four planets: an analytic 1-ubar base column plus a `base.inp` handoff (isothermal or Guillot T(p) generator) — see the folder's own README (multi-planet, not the HD189733b baseline) | `Lower column: <R_1bar>`; driver-generated `base_iso.inp` / `base_guillot.inp` |
| `14_diffusion` | **Diffusive separation of He and metals** (HD209458b): the He/H ratio declines with altitude and each trace metal settles independently, reshaping the He 10830 line | HD209458b params `+ Include He23S? True`, `+ He_diffusion: True`, `+ He_metal_diffusion: True`, `+ He_Kzz: 1.0e9`, `+ He_alphaT: 0.0`, `metals.inp` present |
| `15_molecular` | **Full molecular chemistry** (HD209458b): H2/H2+/H3+/HeH+ in the coupled ionization equilibrium; a sharp H2->H front forms above a thin molecular base (front near r ~ 1.02 R_p), the wind above it essentially atomic | HD209458b params `+ Molecular chemistry: True` (no `metals.inp`; mutually exclusive with metals and `He_diffusion`) |

Notes
- Wind-AE IC (`11`/`12`, `docs/wind_ae_solver.pdf`): `IC mode: windae`
  builds the IC in-process from a shipped seed (`inputdata/windae_seed.csv`,
  symlinked into each folder), via the self-consistent-BC continuation.
  This converges seed-adjacent hot Jupiters (`12`, HD209458b) and
  far-from-seed planets (`11`, HD189733b) alike. Whether EXHALE then
  time-integrates the result cleanly is a separate matter: HD189733b hits an
  EXHALE-side base-breathing instability regardless of the IC source.
- `04_newton_from_state` needs a state to start from: copy a converged
  `Hydro_ioniz.txt` / `Ion_species.txt` (e.g. from `03_newton/output/`) to
  `output/Hydro_ioniz_IC.txt` / `output/Ion_species_IC.txt` first.
- Metals on/off is a runtime switch: any folder becomes metals-on by
  copying a `metals.inp` into it (and metals-off by removing it).
- HD189733b has a "breathing" (slightly inflowing) base; the practical
  Newton residual floor appears to be ||R|| ~ 3e-4 (see
  `docs/steady_solver_memo.pdf`), which is below the default
  `Resid tol: 1.0e-3` and therefore harmless here.

## 13_lower_atmosphere/
Lower-atmosphere connection examples for HD 209458 b, HD 189733 b,
WASP-121 b and WASP-52 b: one `input.inp` per planet + driver-generated
`base_iso.inp` / `base_guillot.inp` handoff files, with a results table and
regeneration commands in its own README.  Full description:
`docs/lower_atmosphere_coupling.pdf`.

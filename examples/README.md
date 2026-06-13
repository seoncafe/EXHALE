# Ready-made example configurations (HD189733b)

One folder per solver/physics combination. Folders `01`–`11` are the
same planet (HD189733b) so the effect of each option can be isolated;
`12_windae_ic_hd209` uses HD209458b (a seed-adjacent planet) to show a
*working* Wind-AE IC, contrasted with the deliberately non-converging
`11`. Each folder is self-contained: `cd` into it and run the repo-root
binary,

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
| `02_two_stage` | Automatic two-stage PLM -> WENO3 marching | `du_th [PLM,WENO3]: 0.5 1.0e-3` |
| `03_newton` | **Recommended default**: two-stage warm-up + JFNK Newton finish to the true steady state | + `Solver: Newton` |
| `04_newton_from_state` | Newton re-convergence of an existing state (no marching) | `Load IC? True`, `Valve eps: 1.0e-4`, `Resid tol: 1.0e-3`; run with `ATES_PTC=1 ATES_PTC_JFNK=1 ATES_PTC_DTAU0=1.0 ../../EXHALE.x` |
| `05_metals` | Trace-metal cooling (solar C/N/O/Mg/Ca/Na/Fe) | `metals.inp` present |
| `06_he23s` | He 2^3S metastable level (He I 10830 line) | `Include He23S? True` |
| `07_balmer_lya` | Non-LTE H(n=2) + Ly-alpha pumping (H-alpha/H-beta) | + `Stellar Teff/radius`, `Deexc heat`, `Jlya escape-prob`, `Stellar Lya flux/halfwidth/boost` |
| `08_full` | Everything on (= the `HD189733b/` planet folder + Newton); feeds `TPM.py` | 07 + `metals.inp` |
| `09_spherical` | Spherical domain instead of the default Roche/L1 truncation | `Domain mode: Spherical`, `Outer radius [R_p]: 10.0` |
| `10_warm_seed_ic` | Warm-seed (hot Parker overlay) initial condition | `Hot Parker IC: 10000` |
| `11_windae_ic` | In-process Wind-AE warm-start IC — **deliberate non-converging case** (HD189733b is strongly bound and far from the shipped seed, so the static-BC ramp stalls; use `IC mode: auto` for this planet) | `+ IC mode: windae`, `+ Solver: Newton` |
| `12_windae_ic_hd209` | In-process Wind-AE warm-start IC that **works** — HD209458b (not HD189733b), close to the shipped seed, so the ramp converges and EXHALE warm-starts cleanly (spherical 10 Rp) | HD209458b params `+ Domain mode: Spherical`, `IC mode: windae`, `Solver: Newton` |

Notes
- Wind-AE IC (`11`/`12`, `docs/wind_ae_solver.pdf`): `IC mode: windae`
  builds the IC in-process from a shipped seed (`inputdata/windae_seed.csv`,
  symlinked into each folder). It suits hot Jupiters near the seed
  (`12`, HD209458b); a strongly-bound far-from-seed planet (`11`,
  HD189733b) stalls the static-BC ramp — use `IC mode: auto` there.
- `04_newton_from_state` needs a state to start from: copy a converged
  `Hydro_ioniz.txt` / `Ion_species.txt` (e.g. from `03_newton/output/`) to
  `output/Hydro_ioniz_IC.txt` / `output/Ion_species_IC.txt` first.
- Metals on/off is a runtime switch: any folder becomes metals-on by
  copying a `metals.inp` into it (and metals-off by removing it).
- HD189733b has a "breathing" (slightly inflowing) base; the practical
  Newton residual floor appears to be ||R|| ~ 3e-4 (see
  `docs/steady_solver_memo.pdf`), which is below the default
  `Resid tol: 1.0e-3` and therefore harmless here.

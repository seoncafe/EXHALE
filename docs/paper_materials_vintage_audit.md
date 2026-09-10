# Vintage audit of the paper and poster materials

**Written 2026-08-30/31. Read-only: no run was repeated, no output regenerated,
no source touched.** Every date below is a file mtime measured in this tree;
every impact figure is quoted from the changelog section that measured it,
and the section is named so the measurement can be found.

The question this answers is not "are the numbers right" but **"which
physics is each stored result standing on, and what has to be re-solved
before it can be quoted again."** The product is a work list for a
re-convergence campaign, not the campaign.

---

## 1. Verdict

**Every production result the paper and the poster rest on predates all six
of the changes that moved a golden or an observable. Nothing in
`paper/`, `poster/`, `python/paper_data.py` or the four planet folders'
current `output/` can be quoted as current.**

Counted directly (`find ... -name Hydro_ioniz.txt` over `HD209458b/`,
`HD189733b/`, `WASP-52b/`, `WASP-121b/` and `benchmarks/`): **91** stored
solutions. Exactly **four** are newer than 2026-08-19, and all four are the
oxygen-chemistry runs `HD209458b/oi1302*` of 2026-08-30 09:03-09:39. The
other **87 are older than section 75** -- the earliest of the six gates. The
most recent of those 87 is the four-planet production set of
**2026-08-19 14:56-16:59**, which is what the paper quotes.

The campaign that follows is far smaller than 91. **36** of those
directories are named `output/`, i.e. they are the current state of their
case; the other **55** are `output_pre_*` / `output_<label>_<date>` archive
snapshots that exist to record a previous state and must **not** be
regenerated. Of the 36 live directories, **15 carry a `metals.inp`** and are
therefore the section-107 targets.

Two findings are separate from vintage and change what the campaign has to
decide first:

- **The LaRT chain and the paper describe different models, not different
  vintages.** Everything Lyman-alpha reads `benchmarks/`, and
  `benchmarks/hd209` runs with `He_diffusion: True`, `He_metal_diffusion:
  True`, `He_Kzz: 1.0e9` and no `Jlya escape-prob`, while `HD209458b/` runs
  with `Jlya escape-prob: True` and no diffusion. Their mass-loss rates
  differ by **0.93 dex** (10.40 against 9.47). Re-converging the paper's
  four models and re-running LaRT on `benchmarks/` would leave that
  disagreement exactly where it is.
- **The WASP-52 b transit curves are older than the winds they were
  synthesized from.** Every `WASP-52b/<case>/tpm_He10830.txt` is dated
  2026-08-09 07:52 while its own `output/Hydro_ioniz_adv.txt` is
  2026-08-15 20:55-23:40. The "companion analysis" factors of 3-10 (L1
  truncation) and about 2 (helium abundance) that `paper/ms.tex` quotes at
  lines 859-871 come from that chain, so they do not describe the winds
  stored beside them.

---

## 2. The six gates: dates, reach, and measured size

Dates are the section headings of `Update_EXHALE_stage1.md`. "Reach" is which runs
the change can touch at all; a run outside the reach needs no re-solve on
that account.

| gate | date | reach | measured size, and where measured |
|---|---|---|---|
| **75** He 2^3S double count | 2026-08-27 | every `Include He23S? True` run | **Negligible for these planets.** `wasp_full` (WASP-121 b): max rel. `drho` 1.24e-6, `dT` 1.32e-7, `dMdot/Mdot` -6.5e-8. The move is the triplet mass fraction, which is 1.3e-6 on a hot Jupiter. Sec. 75, "What moved" |
| **82** three-body M was a mass density | 2026-08-27 | molecular runs only | **Does not reach any of these runs.** No `input.inp` in the four planet folders or `benchmarks/` sets `Molecular chemistry: True` (checked, 0 hits). Sec. 82 |
| **87** Penning rate -> Garcia Munoz | 2026-08-28 | every `Include He23S? True` run | **The He 10830 gate.** `wasp_full`: scale-relative `rho` 8.6e-5, `T` 6.4e-4, heating 2.6e-3, `log10 Mdot` 13.23 unchanged -- but **peak n(He 2^3S) rises 3.1 %** there, 6.9 % on `lower_profile`, and the He 10830 equivalent width of the LHS 1140 b reference winds rose 7.6-18.5 %. The wind barely moves; the metastable line does. Sec. 87 |
| **88** `_adv` population as a difference; `info` never read | 2026-08-29 | every `_adv` file, hence every transit product | **No golden can see it** -- no regression case compares an `_adv` file. On healthy solutions the `_adv` species move 1e-4 to 1e-3; on a run whose advection solve was failing (310 of 4660 cells returning `info /= 1`) the He I and 2^3S `_adv` profiles moved by **0.25**. `EXHALE_transit.py` reads `Hydro_ioniz_adv.txt`, `Ion_species_adv.txt` and `OI_levels_adv.txt`, so this gate sits directly under every `tpm_*` file. Sec. 88 |
| **96** SvS85 secondary branching on the cell's own composition | 2026-08-30 | every run | `wasp_full`, as the whole series 94-96 against the golden: `rho` 2.3e-3, `p` 1.0e-3, `T` 6.8e-4, `H II` 1.0e-2, **`He II` 7.6e-2**, `log10 Mdot` 13.23 unchanged. The 7.6 % is the base cell, where `n(He I)/n(H I) = 0.0851` departs 15 % from SvS85's 0.1; it falls to 0.4 % at 1.2 R_p and 6e-4 at 1.4 R_p. Sec. 96 |
| **107** `cx_O2p_H` on by default | 2026-08-30 | metals-bearing runs only | **The largest mover for these runs.** `wasp_full`: max abs rel. **`rho` 3.9e-2 / `T` 5.0e-2**; `wasp_he23off` 4.0e-2 / 5.1e-2; `lower_profile` (HD 209458 b) 2.5e-2 / 2.6e-2; `mol_metals` 4e-4 (shielded molecular base). O III is suppressed by 4.1e-5 at 1.05 R_p and 5.2e-1 at 1.5 R_p on `wasp_full`; the oxygen lands in O II, which is a coolant and O III is not, so n(O I) rises 3.4x and n(O II) 3.2x near 1.5 R_p, cooling rises 16 % and `T` falls 5 % there. `log10 Mdot` unchanged to two decimals on all four cases. The three metals-free molecular cases stayed byte-identical. Sec. 107 |

**Do not treat these as six patches to apply.** Sections 94, 95, 97, 98 and
100-105 landed inside the same interval, and the LHS 1140 b refresh campaign
states the point plainly: what a re-solve applies is *the current binary*,
not a selected subset (`LHS1140b/exhale/refresh_j96/basemetals/results.txt`,
"sections 94, 95, 97 and 100-105 are also inside the interval and are not
separated out here"). The gate list is for deciding **which** runs must be
re-solved and **why**; the re-solve itself carries everything.

### A note on what the setup report can and cannot tell you

Vintage cannot be read off `EXHALE_setup.out` or `EXHALE_resolved.out`.
Checked in `src/modules/files_IO/write_setup_report.f90`: the first is a
human-readable summary (planet and star parameters, He 2^3S on/off, the
metal element count, `N_eq`, the grid, `Mdot`) and the second is the resolved
geometry and composition the transit tool consumes. The keys that would
identify the physics generation -- `use_sec_ion`, `he_h_charge_exchange`,
`cx_o2p_h_scale`, `ates_photoion_rate` -- appear only in `parse_dump.txt`,
written only under `EXHALE_PARSE_DUMP=1`, and **there is no `parse_dump.txt`
anywhere in the five trees surveyed**. Vintage here is therefore established
from output mtime against the section dates, from the presence of the
`cx_O2p_H` key in `metals.inp`, and from the run logs.

`EXHALE_setup.out` also sits in the run directory root and is overwritten by
every run, so it describes the newest `output/` only -- never the
`output_pre_*` snapshots beside it.

---

## 3. The vintage table

### 3a. The paper's production set -- 2026-08-19, pre-section-75

This is the set `python/paper_data.py` hard-codes and `paper/ms.tex` quotes.

| run directory | `output/` mtime | generation | metals | products standing on it |
|---|---|---|---|---|
| `HD209458b/output/` | 2026-08-19 14:56:07 | pre-75 | on, 7 elements | `paper_data.py` `log10 Mdot = 9.47`, `\|\|R\|\| = 1.77e-4`; `paper/fig_struct_hd209.pdf`; `HD209458b/tpm_*.txt` (7 files, 17:01:42) -> `fig_transit_{He10830,Halpha,Lya}.pdf`; `ms.tex` ll. 721-750 |
| `HD189733b/output/` | 2026-08-19 14:56:11 | pre-75 | on, 7 | `9.14`, `8.35e-4`; `fig_struct_hd189.pdf`; `tpm_*` (17:01:36); `ms.tex` ll. 779-818 |
| `WASP-121b/output/` | 2026-08-19 14:58:36 | pre-75 | on, 7 | `13.20`, `7.75e-4`; `fig_struct_wasp121.pdf`; `tpm_*` (17:01:49); `ms.tex` ll. 885-920 |
| `WASP-52b/output/` | 2026-08-19 16:59:02 | pre-75 | on, 7 | `11.70`, `8.66e-4`; `fig_struct_wasp52.pdf`; `tpm_*` (17:01:47); `ms.tex` ll. 838-849 |

Cost of that set, measured from the `Newton finish at step` and
`Execution Time` lines each `run_20260819_const.log` prints:

| planet | marching steps | wall clock |
|---|---|---|
| HD 209458 b | 2002 | 0 h 00 m 55.6 s |
| HD 189733 b | 2002 | 0 h 00 m 59.1 s |
| WASP-121 b | 8241 | 0 h 03 m 24.6 s |
| WASP-52 b | **178200** | **2 h 03 m 51.5 s** |

All four finished at `info = 0`. **A full cold-start redo of the paper's
four models is about 2 h 10 m of wall clock, and WASP-52 b is 95 % of it.**

All four are metals-on, so all six reaching gates apply: 75 (negligible), 87,
88, 96, 107. Section 82 does not reach them.

Derived products and their own dates: `paper/ms.tex` 2026-08-28 12:47,
`ms.pdf` 12:52 -- i.e. **the manuscript was last compiled nine days after the
runs and two days before the last gate**, and every structure and transit
number in it was verified against the stored outputs at that time. The
numbers are internally consistent; they are consistent with a superseded
binary.

### 3b. The poster and the "companion analysis" -- three vintages in one chain

| layer | date | note |
|---|---|---|
| `WASP-52b/<case>/output/` (22 cases) | 2026-08-15 20:55 - 08-16 00:25 | pre-75; 18 metals-off, 4 metals-on (`*_solar_met*`) |
| `WASP-52b/<case>/tpm_He10830.txt` | **2026-08-09 07:52** | **older than the winds above** |
| `poster/figs/*.pdf` (7 figures) | 2026-08-09 07:34 | older still |
| `poster/WASP-52b_poster.{tex,pdf}` | 2026-08-11 18:44 | |
| `poster/make_poster_figs.py` | 2026-08-29 22:23 | the script has moved since the figures were drawn |

`make_poster_figs.py` reads `WASP-52b/fxuv0p25_{he98,solar}_L1/output/Hydro_ioniz_adv.txt`
for `profiles.pdf`, the `tpm_He10830.txt` of five cases for `he_compare.pdf`
and `hescan.pdf`, and three LaRT `.h5` fields for `halpha.pdf` /
`halpha_matched.pdf`. It also defaults its observation directory to
`WASP-52b/observations/`, **which does not exist** -- the script cannot run
without `EXHALE_OBSDATA` set.

`paper/ms.tex` ll. 859-871 quotes this chain: the L1-truncation factor 3-10
and the helium-abundance factor of about 2. Those come from `tpm_*` curves
that predate their own winds by six days, so they are stale twice over.

### 3c. The LaRT coupling -- a different model, then a vintage problem

`benchmarks/` is what everything Lyman-alpha reads.

| run | `output/` mtime | metals | `log10 Mdot` | LaRT case built from it |
|---|---|---|---|---|
| `benchmarks/wasp121/` | 2026-08-15 18:43:04 | off | 13.07 | `lart_runs/bm_wasp121/`, `insitu_wasp121/` |
| `benchmarks/hd209/` | 2026-08-15 19:04:17 | on, 3 elements | 10.40 | `lart_runs/bm_hd209/`, `insitu_hd209/` |
| `benchmarks/hd189/` | 2026-08-16 00:12:53 | off | 9.58 | `lart_runs/bm_hd189/`, `insitu_hd189/` |
| `benchmarks/wasp52/` | 2026-08-16 00:16:19 | off | 11.87 | `lart_runs/bm_wasp52/`, `insitu_wasp52/` |

LaRT executions: `bm_*` 2026-08-16 01:25 - 11:31 (2e6 photons each);
`insitu_*` 2026-08-16 11:31 - 08-17 08:36 (1-2e5 photons). The eight
`*_pre_kb_20260815/` directories hold the complete July runs (executed
2026-07-06/07) kept aside when the CODATA `k_B` rebaseline landed -- the
precedent that a changed wind means a redone LaRT run, not a reused one.
The `WASP-52b/LaRT_lya*` set (three cases plus three predecessors) ran
2026-08-17 08:40 - 08-19 11:03 off `WASP-52b/fxuv0p25_he98{,_L1}`.

**A re-converged wind forces a LaRT re-run.** The LaRT inputs *are* the
EXHALE profiles: `exhale_to_lart.py` `write_profiles()` slices
`dens_profile.txt` = `n_HI(r)`, `temp_profile.txt` = `T(r)` and
`velo_profile.txt` = `v_r(r)` straight out of `output/*_adv.txt`, and
`insitu_emiss.txt` is built from `n_HI`, `n_HII`, `n_e`, `T` of the same
files; the namelist geometry (`distance2cm`, `a/R_p`, `R_star/R_p`, `rmax`)
comes from the same run's `input.inp`. Nothing in the `.h5` is independent
of them -- `Pa_2D` / `Pa_1D` is the Lyman-alpha scattering rate *in that*
field, at a pole optical depth of order 2e9. The `k_B` rebaseline alone
already moved `N(H I)_pole` from 1.658e22 to 1.431e22, a 16 % change, and
that is smaller than what sections 87-107 will do.

Note that the `.in` files carry a 2026-08-19 22:55 mtime from a
**comment-only** edit (source path and executable name); no parameter
changed, so the stored fields remain valid for their own namelists.

**Settle the model question before the vintage question.** `paper_data.py`
and both figure scripts point at the planet folders; `exhale_to_lart.py`,
`lya_insitu_emissivity.py` and `tpm_halpha_lart2d.py` point at
`benchmarks/`. Re-pointing `benchmarks/` at the planet folders, or keeping
it as a deliberately separate diffusion-on comparison set, is a decision the
code does not contain.

### 3d. The only current-generation runs in the tree

| run | `output/` mtime | generation |
|---|---|---|
| `HD209458b/oi1302/` | 2026-08-30 09:03:53 | past 96; `cx_O2p_H` **off** (the section-104 default-off pair) |
| `HD209458b/oi1302_diffusion/` | 2026-08-30 09:07:27 | past 96; `cx_O2p_H` off; the other half of that pair |
| `HD209458b/oi1302_cxO2p_1/` | 2026-08-30 09:38:59 | **the only run whose physics equals the current default** -- its `metals.inp` carries `cx_O2p_H 1.0`, and section 107 names it as the reference the new default reproduces |
| `HD209458b/oi1302_cxO2p_3/` | 2026-08-30 09:39:01 | scale 3.0 sensitivity arm; no `tpm_*` |

The source of section 107 (`charge_exchange.f90`, `metals_input_read.f90`)
is dated 2026-08-30 12:02, after all four, which is consistent: they are the
runs of sections 101-105 that the decision was taken on.

### 3e. Archive snapshots -- do not regenerate

55 of the 91 directories are archive snapshots kept to record a previous
state, counted by directory name:

```
27  output_pre_kb_20260815
 4  output_pre_rootfix_20260811     4  output_pre_ppadv_20260812
 4  output_pre_newtonfix_20260810   4  output_pre_lyafix_20260812
12  one-off labels, one directory each: output_pre_metals_20260811,
    output_pre_lsfix_20260811, output_pre_gridfix_20260811,
    output_pre_ghostT_20260811, output_nan_20260810,
    output_march1e6_20260811, output_metals_on, output_metals_off,
    output_jun7, output_caseA2058, output_caseA2058_jun7, baseline
```

They are pre-section-75 by construction and that is their purpose. A further
five snapshots hold only `*_adv.txt` with no `Hydro_ioniz.txt` and so do not
appear in the 91 at all; five of the 2026-08-09 recovery set carry the
pre-`# columns` output format.

---

## 4. The re-convergence work list

Counted by what something reads, not by what exists.

| # | what | runs | why it must be re-solved | cost anchor |
|---|---|---|---|---|
| 1 | the four paper models: `HD209458b/`, `HD189733b/`, `WASP-52b/`, `WASP-121b/` top-level `output/` | **4** | metals-on, He 2^3S on: gates 87, 88, 96, 107 all reach them, and 107 is measured at 4 % in `rho` / 5 % in `T` on the WASP-121 b analogue | **2 h 10 m** cold start, measured (section 3a); far less seeded (see below) |
| 2 | transit synthesis on the four | **4** | `EXHALE_transit.py` reads the `_adv` files that gate 88 changed; the metal doublets the paper quotes (Mg II, Ca II, Na I) sit directly on the oxygen redistribution of gate 107 | minutes; `Do only PP` then `EXHALE_transit.py` |
| 3 | the four `benchmarks/` runs | **4** | the LaRT chain's only inputs | as 1 |
| 4 | LaRT re-run on the new benchmark profiles | **8** cases (`bm_*` + `insitu_*`) | inputs are the profiles themselves | 2026-08-16/17 wall clock: `bm_*` ~10 h total at 2e6 photons; `insitu_*` ~21 h at 1-2e5 photons |
| 5 | the WASP-52 b case scan behind the poster and the "companion analysis" | **22** sub-cases (verified by listing; 4 carry a `metals.inp`) | gates 87, 88, 96 reach all 22; 107 reaches the 4 | the expensive block: WASP-52 b is the slow marcher of the four, so cold-starting 22 of them is tens of hours. Seed them |
| 6 | WASP-52 b LaRT (`LaRT_lya`, `LaRT_lya_L1`, `LaRT_lya_L1_fwhm`) | **3** | as 4 | 2026-08-17/19 wall clock: roughly 1 h, 19 h, 22 h |
| 7 | poster figures + `paper/figs/fig_lya_insitu.pdf` | 8 figures | drawn from 1-6 | minutes, once 1-6 land |

**Minimum for the paper alone: items 1 and 2 -- four winds and four transit
passes.** Items 3, 4 and 6 are the Lyman-alpha side and are blocked on the
model decision of section 3c. Item 5 is what the "companion analysis"
paragraph of `ms.tex` needs, and it is the largest single block.

### The recipe already exists

`LHS1140b/exhale/refresh_j96/resolve_wind.sh` is the tested pattern: copy
`input.inp` **and every runtime-configuration file the source carries**
(`lower_atmosphere_profile.dat`, `base.inp`, `metals.inp`, `opacity.inp` --
EXHALE is configured by file presence, so a dropped file silently changes
the run), seed `output/` from the stored solution, run N PTC/JFNK passes
with `EXHALE_PTC=1 EXHALE_PTC_JFNK=1`, then one `Do only PP` pass, then
`EXHALE_transit.py`. Section 106 ran **82 JFNK re-solves in seven series
plus three end-to-end flux closures** through it in a day, on winds seeded
from their own converged state. Seeding is what keeps item 5 affordable: a
cold WASP-52 b start is 178200 marching steps and two hours, and item 5 is
22 of them.

Two traps that campaign already paid for and this one should not:

- **A `Do only PP` re-run is not free.** It advances the stored solution by
  one RK step before post-processing, which by itself moved a stored
  equivalent width by 4.4e-4 (section 107). Compare new against old on two
  runs that took the *same* path, never against the file sitting in the
  directory.
- **A converged-looking `info = 0` at a loose `Resid tol` may be a wind
  that took zero Newton steps.** Section 106 found one arm at
  `Resid tol: 5.0e-3` handing back `||Fs||2 = 1.27` against 0.018-0.074 on
  its siblings. Check the achieved residual, not only `info`.

And one property of these solutions that no re-solve removes: section 106
measured a **seed hysteresis of 4.6 %** in equivalent width -- one identical
configuration reaching three different fixed points depending on which
converged wind it started from, each at `info = 0` and `||R|| ~ 1e-4`. That
is larger than what gate 96 is worth on the same arms. It was measured on
LHS 1140 b and has not been looked for on these four planets.

---

## 5. How much will sections 96 and 107 move each planet?

Neither has been measured on these four runs -- that is what the campaign is
for. What **has** been measured, and is the best available prior:

- **Section 96, on the WASP-121 b regression case** (`wasp_full`, the same
  planet and the same metals-on, He 2^3S-on configuration as
  `WASP-121b/output/`): `rho` 2.3e-3, `p` 1.0e-3, `T` 6.8e-4, `H II` 1.0e-2,
  `He II` **7.6e-2**, `log10 Mdot` 13.23 -> 13.23. The 7.6 % is the base
  cell and decays outward to 6e-4 by 1.4 R_p, tracking the departure of the
  local `n(He I)/n(H I)` from SvS85's 0.1. The two WASP-121 b golden cases
  were re-snapshotted; five other cases stayed byte-identical because their
  composition sits at the SvS85 value.
- **Section 107, on the same case**: `rho` **3.9e-2** / `T` **5.0e-2**
  (`wasp_full`), 4.0e-2 / 5.1e-2 (`wasp_he23off`), 2.5e-2 / 2.6e-2 on the
  HD 209458 b `lower_profile` case, 4e-4 on the shielded molecular base.
  `log10 Mdot` unchanged to two decimals throughout -- the escape rate
  absorbs it, the temperature does not. The mechanism is stated rather than
  absorbed: oxygen removed from O III lands in O II, which is a coolant and
  O III is not, so near 1.5 R_p n(O I) rises 3.4x, n(O II) 3.2x, cooling
  16 %, and `T` falls 5 %.

**Expected reach on the four planets.** WASP-121 b has a direct analogue in
`wasp_full` and should move by about those amounts. HD 209458 b has one in
`lower_profile` (2.5 % / 2.6 %) though that case carries its oxygen through
a photochemical column rather than `metals.inp`. HD 189733 b and WASP-52 b
have no measured analogue at all; both are metals-on with solar oxygen, so
the `wasp_full` figure is the prior to carry, not a result.

**A 4-5 % move in `T` is not a 4-5 % move in a line depth.** The metal
doublets `ms.tex` quotes (Mg II 11.8/6.9 % on HD 209458 b, 71.2/62.8 % on
WASP-52 b; Ca II 3.6 % and 16.9 %) sit on the same metal ionization balance
that gate 107 redistributes, and the He 10830 depths sit on a metastable
population that gate 87 raised by 3.1 % at peak on the WASP-121 b analogue.
Those are the numbers to watch, and none of them can be predicted from the
hydro move.

**LHS 1140 b is the counter-example worth keeping.** Section 107 was
measured there and the red-pair equivalent width moved by **-4.2e-7** in
relative terms, because the photochemical column's cold trap has already cut
`X_O` by a factor 685 at the matching level. The size of a gate is a
property of the atmosphere, not of the gate.

---

## 6. H2 collision-induced absorption: what activation would need

**No implementation is proposed here.** This section records where the
exclusion was decided, what the physics would be, where in the code it would
go, what data exists, and over what range it is valid.

### The decision on record

`docs/oxygen_chemistry_new_plan.md` section 4.1: *"(G)'s H2 CIA component is
already decided out of the paper's scope by the user; nothing here reopens
that."* CIA is the continuum half of item (G) of `TO_BE_DONE.md`, and (G)
states the gap precisely: the metals-on molecular layer, and the metals-off
one beyond the H2 -> H front (229.7 K at r = 1.095 on the hot Uranus case),
sit far below any radiative equilibrium temperature, and the `Base IR field`
closure of section 55 cannot reach them, because it gives the **existing**
line coolants their incident field and past the front there are no molecular
line coolants left. The whole radiative budget between the base and the
front is line channels with a radiative time of 3e8 s against a flow time of
5e6 s.

### What it would be, and where it would go

CIA is a **continuum** thermal-infrared opacity, so it belongs on the
cooling and energy side, not on the stellar-photon side:

- **Not `opacity.inp`.** That schema
  (`docs/input_schema.md` section 5) carries exactly four species scale
  factors and four table paths -- H I, He I, He II, He 2^3S -- plus the
  Robinson & Catling pressure-broadening parameters. It is the photo-opacity
  of the ionizing continuum and has no thermal-infrared band.
- **The natural home is the `Base IR field` path**:
  `parameters.f90:313-323` (`base_ir_field`, default `.false.`),
  `Cool_coeff.f90:2074, 2279, 2342` (`fine_structure_line_transfer`),
  `util_ion_eq.f90:959-1065`, `T_equation.f90:33`, with the key parsed at
  `input_read.f90:657-668` and echoed at `write_setup_report.f90:208-212`.
  That path already computes a diluted `B_nu(T0)` incident from below and
  already decides line by line whether the gas sees it. A continuum term is
  the same construction with a wavelength-resolved opacity instead of a line
  profile, and it would supply exactly what (G) says is missing: a coolant
  and an absorber that survive past the H2 -> H front.
- The same key would then control both, which is the right shape -- the
  question "does the gas see the atmosphere below it" has one answer, not
  two.

### Data, already on disk

`photochem_clima_data/photochem_clima_data/data/CIA/` carries 14 HDF5
tables, including the two that matter: **`H2-H2.h5` and `H2-He.h5`**.
Measured directly from the files rather than quoted:

| | value |
|---|---|
| shape | `log10xs(500, 30)`, `wavelengths(500)`, `T(30)`, plus a `notes` string |
| wavelength range | **0.607 - 247.30 um** |
| temperature range | **100 - 3000 K** |
| units | `cm^2/molecule` (stated in the file's own `notes`) |
| provenance | Molliere et al. (2019), the petitRADTRANS compilation, itself drawn from HITRAN |

The temperature range covers the whole problem: the collapsed molecular
layer sits at 190-300 K and `T_eq` is 1140 K on the hot Uranus case, 397 K
on the `examples/15_molecular` HD 209458 b case. `H2-CH4`, `CO2-H2`,
`N2-H2`, `H2O-H2O` and `H2O-N2` are in the same directory if the H2O/CH4/CO
bands (G) also names are wanted later.

The alternative source is HITRAN CIA directly (`hitran.org/cia/`), which is
what most of these tables came from. The local files' one stated limitation
is that H2-H2, H2-He and N2-N2 stop at 250 um while HITRAN has absorption
beyond it -- a small error for emitters below about 250 K, which is the
regime of the collapsed layer, so it is worth checking rather than assuming.

### Validity, and the reason it was excluded

`references/` already holds the anchor: **Lavvas & Arfaux (2021, MNRAS 502,
5643)** treat middle-atmosphere thermal structure with CIA in the radiative
transfer and note that CIA becomes significant at **p > 1 bar**. The EXHALE
domain begins at the 1 microbar match (a 9.0 microbar base on the hot Uranus
case) -- six orders of magnitude lower in pressure than where CIA starts to
matter. **CIA's own opacity inside the EXHALE domain is therefore
negligible, and that is the argument that excluded it.**

What (G) asks for is a different quantity: not CIA as a local coolant but
CIA as the source of the **incident continuum field** the deep atmosphere
radiates upward into the domain, which is exactly where the column is
optically thick. So an implementation would be a downward boundary condition
built from the deep column's CIA optical depth, not a term added to the
local cooling sum. Getting that distinction wrong is the way to spend the
effort and change nothing.

Two further code-level anchors already in `references/`: **Robeling et al.
(2026, Kompot)**, a 1-D self-consistent thermochemical upper-atmosphere model
carrying H2-H2 and H2-He CIA (Abel et al. 2011) plus H2O/CO/CO2/CH4 in the
radiative budget -- the closest published analogue of what this would build;
and **Koskinen et al. (2013a)**, the paper-level original, which includes
H3+, CO, H2O and CH4 as strong infrared coolants on HD 209458 b.

One piece is already in the tree and is not the same thing:
`src/utils/radiative_convective_column.py` runs Photochem's `clima`, whose
correlated-k radiative transfer **does** include CIA from these tables. That
is the lower-atmosphere climate solve that produces `T(p)` for the handoff
profile -- below the EXHALE domain, not inside it. It does not close (G),
and it does not make (G) cheaper except that the tables and their reader are
already here.

---

## 7. What this audit did not check

- **Whether the Python tools changed their output.** `EXHALE_transit.py`
  (2026-08-30 09:46), `exhale_transit_lib.py` (09:43) and
  `examples/exhale_io.py` (15:14) were all edited after the
  2026-08-19 17:01 `tpm_*` curves were written, as were
  `examples/exhale_to_lart.py` (2026-08-29 22:23) and
  `lya_insitu_emissivity.py` (2026-08-27 21:03) relative to the LaRT inputs
  they produced. Whether any of those edits move a number was not
  determined; it cannot be settled by reading the current files alone.
- **The size of gates 87, 88, 96 and 107 on HD 189733 b and WASP-52 b.**
  No measured analogue exists for either; section 5 carries the
  `wasp_full` figure as a prior and says so.
- **Anything in `examples/`, `LHS1140b/` or `backup/regression/`.** The
  LHS 1140 b tree was refreshed onto the current binary by section 106 and
  the regression goldens by section 107; neither is in scope here.

## 8. Loose ends found in passing

Nothing here was touched -- this audit is read-only -- and all three are
recorded so they are not discovered mid-campaign.

- **`benchmarks/hd209/input.inp` is newer than its own output.** Its mtime
  is 2026-08-25 12:21, ten days after the `output/` beside it
  (2026-08-15 19:04). The stored output was not produced by the `input.inp`
  now sitting next to it, so re-running that directory as it stands does not
  reproduce the stored result.
- **`poster/make_poster_figs.py` points at a directory that does not
  exist.** It defaults its observation path to `WASP-52b/observations/`
  (checked: absent), so the script cannot run without `EXHALE_OBSDATA` set.
- **`paper/figs/fig_lya_insitu.pdf` is an orphan.** Dated 2026-08-19 17:02
  and produced by `examples/lya_insitu_emissivity.py`, it is referenced by
  no `\includegraphics` in `ms.tex` (checked: 0 matches).

---

## 9. Addendum 2026-09-03: a seventh reason, and it is not a vintage

**stale: computed with ghost rows included (P48, Update_EXHALE_stage1 section 137);
depths move 0.2-3.8%; regeneration awaits instruction.**

The six gates above are physics changes that a re-convergence would answer.
This one is not: until 2026-09-03 every Python reader of an EXHALE profile
returned the file's **ghost rows** -- the two boundary rows at each end of
every `Hydro_ioniz*.txt` / `Ion_species*.txt`, a fixed base state below and a
zero-gradient / WENO3 extrapolation above -- as if the solver had converged
them, and `EXHALE_transit.py` and `EXHALE_plots.py` had the same trap
independently. So every stored transit curve was integrated over a domain two
cells too tall, and every stored profile figure carries one boundary point at
each end of every curve. Re-solving does not fix it; re-running the readers
does, and they are fixed now.

**Measured** on three converged states (section 137.2 of the changelog):
line-center depths move by -3.81 per cent (hot Uranus Ly-alpha, the largest
single move), -2.55 per cent (HD 189733 b Ly-alpha), -1.95 per cent (hot
Uranus He 10830), -1.69 to -1.62 per cent (WASP-121b He 10830 / Ly-alpha /
Mg II / Ca II), and by less than 0.2 per cent on the weak lines; 4 A band
depths by -2.5 to +0.5 per cent.

**Marked in place, values untouched, nothing regenerated:**

- `paper/ms.tex` -- eight `%% STALE` comments: the abstract's depth range, the
  four per-planet results paragraphs that quote depths, and the three transit
  figures (`fig_transit_He10830/Halpha/Lya.pdf`).
- `poster/WASP-52b_poster.tex` -- appended to the existing `%% Model vintage:`
  block as a second, independent reason.
- `poster/make_poster_figs.py`, `examples/make_figures.py`,
  `python/make_transit_figures.py` -- docstring notes; all three are producers
  of stale figures.
- `docs/transmission_spectrum.tex` -- three sites (the H-alpha/H-beta depth
  example, the validation section, the metal-doublet figure).
- `docs/EXHALE_user_manual.tex` -- the He 10830 stellar-radius example and the
  `tpmspec` transit figure.
- `docs/figures/README_STALE_P48.md` -- the affected figure files, split into
  those whose depths move and those where it is the plotted range.
- `HD209458b/`, `HD189733b/`, `WASP-52b/`, `WASP-121b/`, `LHS1140b/` --
  `STALE_P48.md` in each, naming the count of saved `tpm_*` files and the
  notebooks whose stored outputs are affected.

**Not in scope of the marking:** `benchmarks/`, `backup/`, and the
`output_pre_*` / `tpm_*_<date>` archive snapshots, which exist to record a
previous state and must not be regenerated.

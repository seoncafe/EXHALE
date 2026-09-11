# Session handoff, 2026-09-07

State of the tree at the end of the PLAN rev2 series (Steps A, B, C and the
items that grew out of them; `docs/Update_EXHALE_stage2.md` section 7 is the record,
`docs/Update_EXHALE_stage2.{tex,pdf}` its generated twin).

## Verified state

- Shared build: `EXHALE.x` md5 `ba3d34e6eaee`, gfortran 16.2 (conda) +
  OpenBLAS, `make -q` up to date.
- `make test`: every suite green (the two steady suites need
  `EXHALE_STEADY_RUNLOGS=backup/regression/wasp_full_newton/run.log`; the
  `attempted_step` whole-binary rows need `EXHALE_ATTEMPTED_STEP_EXE`).
- Goldens: `backup/regression/golden/` refreshed 2026-09-07 13:19 from this
  gate's outputs; the previous set is `golden_pre_series_20260907/`.
  `make check`: REGRESSION PASS, 16 cases, 64 files byte-identical
  (13:21-16:25 KST). `run_fcheck.sh`: CLEAN.
- Movement of every case against the previous goldens: the table in the
  "Series gate" entry of `Update_EXHALE_stage2.md` (log10 Mdot within +-0.05 dex;
  `wasp_full_newton` 13.30 certified; `wasp_he23off` exits 2 by design, its
  stationary claim refused at the default tolerance).

## Cost

Since B3c (the coupled source step) a marching step costs 5-10x: `wasp_full`
to its du stop 3 h at one thread; the 12000-step molecular cases 40-80 min
each. `make check` runs the cases eight at a time. COST3 measured the lever
(the cell's own optical depth in its residual, `System_HeH*`); it is B5 work.

## Open items, in order

1. B5 (species rows, PTC certification) with the COST3 lever.
2. THREAD-DET: which routine turns a thread-count-dependent block boundary
   into a moved last bit (6.6e-14 between 1 and 16 threads).
3. DOCS-LINES: line-number citations in the 2026-09-06 design documents
   drifted by up to a thousand lines; re-anchor on routine names.
4. HYG-PY: `python/paper_data.py` `HEAT_FIXED` and `make_struct_figures.py`
   channel selection lack the three new heating channels (21 columns now).
5. B4 design's open measurement: the converged 12000-step element-separation
   pair for B4-1; WENO3-REF's recorded `eps` prescription experiment.
6. `wasp_he23off` at the default `Resid tol`: leave as the refused-claim
   (exit 2) case in the matrix, or state a tolerance (user's call).

## User-gated

- D-physical reruns of the planet folders and `benchmarks/` on the loaded
  SEDs with the new physics, Newton-finished: every H-alpha/H-beta transit
  number in the repository from a `Jlya escape-prob: True` run is superseded
  by LYA-BETA (n_2p up to 84x at the base), and the metal changes move Mg II
  (-84 percent at 1.18 Rp on `wasp_full` 300). He 10830 is not affected
  through those routes.
- `Load IC` cases (`jfnk_hd189`, `jfnk_hd189_tight`, `ptc_warm`,
  `armD_D2_newton`, the planet folders, `benchmarks/koskinen2022_model_a/*`)
  need their IC regenerated on the current grid (grid guard).
- A converged run reaching 3000 K to check the CO domain nesting.
- Paper/poster re-convergence, LHS 1140 b write-up, GitHub public switch:
  only on explicit instruction.

## Not in git

Every `src/tests/*/` suite of the series (19 directories), the 2026-09-05/06
design and plan documents, `src/utils/h2_shielding_lbl/`,
`src/utils/update_log_to_tex.py`, `references/` additions (11 papers, the
Verner, Mao and Badnell data directories, Krause 1979, Locci, Cecchi-Pestellini,
Fox, Harrington, de Jong, Osterbrock), `docs/session_handoff_20260907.md`.
`git status --short | grep '^??'` lists them.

## Decision document

`docs/To_be_determined_by_user_20260906.md`: all seven items decided
(1 b, 2 a, 3 accepted, 4 C, 5 executed, 6 B, 7 b).

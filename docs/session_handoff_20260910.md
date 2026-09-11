# Session handoff, 2026-09-10

The state to restart from. Read `docs/code_status_20260910.md` first, then
`docs/ISSUES_20260909.md` (sections 3.1 and 5 carry the stage-2 solver
account to its close), `docs/PLAN_20260909_rev1.md` (the plan of record),
`docs/To_be_determined_by_user_20260906.md` (decisions 1 to 23) and
`docs/Update_EXHALE_stage2.md` section 7 (the dated record, items B5e to N38). The
previous handoffs (`session_handoff_20260908.md`, `_20260907.md`,
`_20260905.md`) are records and are not updated.

## Verified state

- Tree: branch `v1.00`, commit `3414478` (2026-09-10 12:06 KST, "Stage 2
  species-row solver, restart contract, certification by regime", 664
  files) pushed to `origin/v1.00` (the remote's default branch; a `main`
  created by mistake the same minute was deleted). On top of it,
  UNCOMMITTED at the time of writing: items N31 to N38 (see below) and the
  documents of this closing.
- Goldens: refreshed 2026-09-10 00:40 (Stage D gate; previous set
  `backup/regression/golden_pre_stageD_20260910/`), `lower_profile` again at
  11:10 for N29 (`golden_pre_n29_20260910/`); `make check` PASS on the
  shared build at that gate. Nothing since N30 moves a default path: every
  N31 to N38 change is a default-off hook, option or key, verified bitwise on
  `backup/regression/wasp_full_newton/IC` and on the first 40 iterations of
  the atomic reload after each item.
- Suites: 25 `src/tests/*/run.sh`; `krylov_and_dogleg` grew 174 -> 264 rows
  (N31 to N36) plus N37/N38's rows; `residual_determinism` has one row that
  FAILs by design (`closure_spread_within_the_row_tolerance_atomic_elem_newton`).
- Fixtures: `backup/regression/atomic_elem_newton` (README with a 2026-09-10
  state section), `backup/regression/carrier_elem_newton` (NEW, N31: the
  carrier reload pinned with its control), `backup/regression/wasp_full_newton/IC`
  (the certified atomic state, the bitwise guard).
- Documentation brought to the stage-2 state 2026-09-10 (DOC1 to DOC3, in
  the commit): the seven physics tex documents and their PDFs, `README.md`,
  `README_HOWTO.md`, `docs/code_status_20260910.md`, `TO_BE_DONE.md`,
  `steady_solver_design.md` section 22 (22.4 corrected after N34).

## What this session settled (N31 to N38, all measured, nothing adopted)

The species-row stationary solves (atomic element reload, molecular carrier
reload) do not converge, and the reason is now stated by measurement rather
than by reading:

1. The finite-difference Jacobian action is not additive along the
   preconditioned Krylov directions because the residual carries a
   ROUNDING FLOOR (N31, N33): the flux difference over a base cell of width
   1.955e-4, amplified by the near-hydrostatic cancellation at the face
   (`epsilon x face state x r^2/dV`, followed over 4.6 decades along the
   column). It follows no inner tolerance (N32). Quadruple precision inside
   the assembly lowers it by 2.1 only, onto the double representation of
   what the assembly is handed (N34).
2. With a faithful image (N34's quadruple-precision assembly) the Krylov cycle STILL stalls: the
   banded preconditioner misses no coupling, the preconditioned operator is
   near-singular on the species rows themselves (Ritz ratio 6.6e5 carrier,
   1.1e4 atomic; the binding carrier row's diagonal 3.3e3 below its coupling
   to the hydrodynamic unknowns of its own stencil), and a longer cycle
   returns a WORSE true step (optimum 60 to 80 products) (N35).
3. The linear model's own row scaling and a true-residual cycle both fail
   acceptance (N36). The linear algebra is exhausted as a lever.
4. Decision 23 (user, "(1)" = the recommendation): a well-balanced flux
   difference for the near-hydrostatic layer, as the default-off key
   `Well balanced:` (N37, Kappeli and Mishra 2014/2016 local hydrostatic
   reconstruction; PDFs in `../references/`). N38 measures why the binding
   species row is nearly independent of its own unknown. RESULTS: N37 is exact
   on the discrete equilibrium (4e-14 of rho g) but the rounding floor
   does not move and the key stays off (a `store_row_terms` prerequisite
   for any default-on, ISSUES 3.7); N38 shows the small diagonal was a
   column-scale reading (sound speed against Mach 1e-3 to 1e-5) and
   column scalings are invisible to the cycle. Both logged.
5. The user closed the program 2026-09-10 evening: "어떤 것을 해도, 끝이
   없어보이니, 엄청난 차이가 나오지 않는 한 이번 건만 마무리 하고 끝내기
   바래." No further solver item is to be started without a new instruction.

## A note for the next session: the different approach, if one is taken

The full analysis is `docs/solver_approach_analysis_20260910.md` (what
was measured lever by lever, the diagnosis, four alternatives compared,
the one measurement that would discriminate between them). Short form:

The user's own reading at the close ("아마도 완전히 다른 접근이 필요할 지도
모르겠음") matches the measurements. What failed is not the hydrodynamic
solve (the three-unknown route converges and certifies; its rows are well
conditioned) and not the species transport operators, but their COUPLING
in one Newton-Krylov space: the species block is near-singular there, and
N21 showed the marching path relaxes what the coupled Newton cannot. Two
routes, neither started:

- **A segregated (partitioned) solve**: the hydrodynamics by the certified
  JFNK, each element and carrier by its own implicit tridiagonal transport
  solve (the "tridiagonal Newton of the implicit step" already exists in
  `binary_element_diffusion.f90`), an outer Picard iteration with
  under-relaxation to the joint fixed point, certification unchanged. Linear
  convergence instead of quadratic, which the present solve does not enjoy
  anyway (chaotic at the ulp level, N26c). Two variants: a Picard fixed
  point between the two solves, or the CETIMB way (READ from the session
  memory of Koskinen et al. 2013 and Huang et al. 2023, to be re-checked in
  the published method sections before it is cited): no stationary Newton
  at all, the time-dependent equations marched to steady state with the
  species continuity equations and the hydrodynamics advanced by operator
  splitting in each step (van Leer advection, Crank-Nicolson for the stiff
  diffusion and conduction terms, a periodic Shapiro filter against the
  base sound wave). In EXHALE that would mean taking the species rows out of
  the Newton unknowns, advancing each species by an implicit tridiagonal
  step inside the present marching, and finishing the hydrodynamics alone
  with the certified JFNK; N21's observation that the marching relaxes what
  the coupled Newton cannot points the same way. Most of the present code
  is reused.
- **A different treatment of the layer**: the near-hydrostatic layer as a
  boundary condition rather than a resolved region (Koskinen's models start
  at 1 microbar), or a mass or log-pressure coordinate that removes the
  `r^2/dV` amplifier. Physically the deeper change; grid, boundary, initial
  conditions and goldens all move; a rewrite.

## Open items, in order (short form; evidence in ISSUES 3.1 to 3.8)

1. The species-row stationary solve (above; closed for now).
2. The layer's element flux conservation (9.8 in the layer, 2.7e-2 in the
   wind; certification gates the wind only, decision 22).
3. The element operator's discretization order (1.59; N = 1000 halves the
   row).
4. Hygiene for a gate: the global `count` shadowing the intrinsic `COUNT`
   (`parameters.f90` line 38); `newton_dense`'s hardcoded step test; the
   dead `dp_bc`/`ntot_bc` refreshes after the sweep
   (`ionization_equilibrium.f90` 2826 to 2858); `golden/mol_ir_bands/`
   without `Cooling_breakdown.txt`; `examples/README.md`; the `.md` twins of
   `lower_atmosphere_coupling` and `transmission_spectrum` parted from
   their `.tex`.
5. The `_adv` product on non-stationary states; the transit numbers in the
   repository predate schema 2 (marked stale where quoted).

## User-gated (unchanged)

Physical reruns of the planet folders and `benchmarks/` (every H-alpha and
H-beta transit number from a `Jlya escape-prob: True` run is superseded by
LYA-BETA, and Mg II by the metals work); `Load IC` case regeneration on the
current grid; a converged run reaching 3000 K for the CO domain; paper and
poster re-convergence, LHS 1140 b, the GitHub public switch.

## Not in git

`backup/` (goldens, fixtures, the archived golden sets), `build*/`, `*.x`,
the planet folders, `TO_BE_DONE.md` (ignored by `.gitignore`, which is itself
untracked and is left alone), the scratch directory of this session
(`/tmp/claude-1000/.../scratchpad`, briefs `planrev1_brief_N*.md`, reports
`planrev1_N*_report.md`, run directories).

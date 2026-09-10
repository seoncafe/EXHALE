# Common brief for the independent source audit of EXHALE v1.00 (2026-09-05)

READ FIRST: `~/.claude/CLAUDE.md` and `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/CLAUDE.md`, and obey every rule there. Key rules repeated:
- Answers come from the actual code, not from comments, docs or memos. A comment is a claim; check it against the code it describes.
- Do NOT modify any file under `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/` (source, inputs, outputs, goldens, docs). This is a READ-ONLY audit. Do not run `make` in the repository, do not run EXHALE.x, do not run the regression matrix.
- You MAY compile small probe programs in your own scratch directory (given in your task prompt) with `gfortran -O0 -g -fcheck=all -fbacktrace -fopenmp -J. -I.` linking the production `.f90` files you need (see `docs/audit_20260905/run_checks.sh` for the pattern). Probes are optional; use them when a claim needs a number.
- Forbidden wording anywhere in your report: the words the global CLAUDE.md bans, and the hyphenated "per-X" noun pattern (write "rate of each cell" instead of the hyphenated form). Do not use em-dashes. Write in English, American spelling.
- Do not speculate. Label every item CONFIRMED (you read the code path end to end and, where numbers are involved, evaluated it) or SUSPECTED (you see a likely problem but could not close the argument). Say what you did not check.
- Report only NEW findings. The list below is already known; do not re-report it (but DO report if you find that a known item is worse, or reaches a caller nobody noticed).

## Repository
Root: `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` (HEAD 35d9dd5, clean tree apart from docs). Source under `src/`. Published papers under `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/references/` (use `pdftotext -layout` when you need to check a coefficient or an equation against its source). Documentation you may consult for context (not as evidence): `docs/code_status_20260905.md`, `docs/physics_numerics_audit_20260905.md`, `docs/development_plan_20260905_rev3.md`, `docs/input_schema.md`, `README.md`.

Code conventions worth knowing: Fortran, free form, `implicit none`; global parameter block in `src/modules/init/parameters.f90` (note the global temperature scale `T0`, case-insensitive: a local `t0` shadows it; the dummy names `G2s`/`G2p` once shadowed module parameters `g2s`/`g2p`, so check every dummy/local name against module scope); adimensional units in the hydro (`T0`, `q0`, code time) with conversions in `UW_conversions.f90`, `composition.f90`; cells `1..N` with `Ng` ghost cells each side (`1-Ng:N+Ng`); OpenMP over cells in the ionization sweep, cooling, hydro stages, carrier Jacobian (check threadprivate/shared correctness and reduction order where you see `!$omp`).

## Already known (do not re-report)
Hydro/numerics: LLF speed `max(|v+a|)` omits `v-a` (Num_Fluxes.f90 287, also reached by positivity_limited_fluxes); ROE two-rarefaction uses sqrt instead of exponent 1/z, shock densities not multiplied by upstream rho, PVRS branch entered whenever p_max/p_min <= 2 giving NaN for equal-pressure expansion, star values feed only the entropy fix in Num_flux, production ROE flux not positivity preserving for vacuum-producing expansion; marching loop rebuilds p at fixed T after the ionization sweep (composition projection changes thermal energy without a source, EXHALE_main.f90 ~1168-1212); energy_semi_implicit does 2 fixed iterations, recomputes cooling only after iter 1, returns cooling of the previous iterate, no residual test, |dC/dT| damping; retry_step covers hydro RK stages only, carrier step and write_back have no status/retry; positivity counter n_faces_flux_positivity_limited counts repaired interfaces only.
Carriers (diffusive_photochemistry.f90): H+ carrier (ic_Hp=5) uses oxygen reference scale nO_free in Jacobian step and residual floor; pct_newton_resid set before limit_to_element_budget; CO has src=0 and an equilibrium ceiling without kinetics, the limiter rescales OH/H2O/H2/H+ too, pct_co_ceiling reset each call and counts only over>1e-10; proton molecular diffusion set to zero (no ambipolar term); carrier advection uses cell-centered derivative and r^2 dr volume, base inflow direction from wind average of rho v r^2; hard-sphere Blanc mixture diffusion; transported-H+ constraint applied only when bg_ready or restart.
Steady solver: self-consistency loop stops on 1e-3 change of norm (steady_newton.f90 ~2636-2655); JFNK action is finite difference; unknown packing nvar_jac=3 or 4 (one H2 row), no H+ row; input_read refuses Ionization transport with Newton/coupled carrier/atomic gas; eval_residual keeps interior U, rewrites ghosts, assemble_residual recomputes pressure via U_to_W under caloric EOS.
Radiation/chemistry energy: H2 double photoionization charges 51.4 eV threshold, chemical product energy is 31.675 eV, 19.7 eV fragment energy has no recipient; double channel reuses single-dissociative degradation row; neutral dissociation channel gives whole remainder to heat (excited fragments H(2p)+H(2s) not represented); 51.4 eV is a vertical threshold, appearance ~47-48 eV; h2_photo_channels reads Chung 1993 Table II as (N_S+2N_D)/N_M but detector gives one pulse per double event so it is (N_S+N_D)/N_M; oxygen network has no collisional reaction-heat ledger, O(1D) energy absent; He(2^3S)/HeH+ associative heat omitted; recombination-cooling fits vs ionization-potential cooling mapping to a full ledger not derived.
H3+ cooling (h3p_cooling.f90): density clamped to 1e6 cm^-3 (nonzero cooling at n_H2=0), polynomial joins discontinuous (43% at 300 K), non-LTE factor jumps 0.9985->1 at 1e14, LTE fit follows Miller unscaled Table 4 while Table 6 matches scaled, evaluated down to 30 K while Table 4 starts at 500 K.
Thermo: q_rovib_H2 (mol_rates) and caloric EOS use different H2 level models (identity u = kT^2 dlnZ/dT violated 0.6-8%); Shomate H2 in oxygen reverse rates is a third representation; H2O/CO/OH/molecular ions have no rovibrational heat capacity; IR emission/absorption tables not exactly balanced at Tgas=Trad, W=1.
Other: molecular post-processing (_adv) excludes molecules; lower-profile schema molecules not evolved; Knudsen/coupling diagnostics not acted on; He+ and molecular ions local equilibrium in hydro; ATES-style post-processed advection for the atomic line; base limit cycle at default CFL; wasp_full_newton stalls at ||R||~5e-4; spectrum file lacks <5 A; literal H I/He I ionization potentials in eval_cool; callerless Cool_coeff wrappers.

## What to look for (new)
1. Wrong physics or mathematics: sign errors, wrong exponents, wrong constants or units (cgs vs code units, eV vs erg, per steradian vs total), factors of 2 or 4 pi, wrong statistical weights, rate coefficients whose formula or coefficients do not match the cited source (check the PDF), detailed-balance relations that do not close, wrong temperature (K vs adimensional) passed to a routine.
2. Conservation: element, charge, mass, energy, momentum not closed in a path (including ghost cells, boundaries, operator splits, restart round trips).
3. Index/array errors: off-by-one in `1-Ng:N+Ng` loops, ghost cells treated as interior, stale arrays read (a value computed for the previous state used for the current), uninitialized variables, save/threadprivate misuse under OpenMP, race conditions, reduction order dependence where determinism is claimed.
4. Name shadowing (case-insensitive) of module/global names by dummies or locals.
5. Dead, unreachable, or inconsistent paths: options that silently do something else than their key says, defaults in input_read that disagree with the setup report or docs, guards that can never trigger, floors/clamps that change physics without being reported.
6. Numerical hygiene: 0/0, log(0), sqrt of negative, catastrophic cancellation, tolerances in absolute units that break with unit changes, NaN propagation paths.
7. Comments and headers that contradict the code they describe (report as documentation defects, low severity, but report them).
8. Interfaces: argument order/kind mismatches across module boundaries, optional arguments assumed present, intent(in) arrays modified.

## Report format (write to the file named in your task prompt, then also return the same text)
For each finding:
```
### <short physical title>
Status: CONFIRMED | SUSPECTED
Severity: P0 (wrong answer in a default path) | P1 (wrong answer in an opt-in path, or conservation) | P2 (hygiene, docs, dead code)
Location: file:line-line (quote the exact lines)
What the code does:
What it should do (with the source: paper, equation, or invariant):
How it was verified (read the callers; probe output; PDF check):
Reaches: (which callers/options/default? which regression cases would move?)
```
End with: "Files read completely: ..." and "Not checked: ...". Be exhaustive within your group; read every file in the list fully, not by grep. The goal is completeness; a report with zero new findings is acceptable only if you state you read every line.

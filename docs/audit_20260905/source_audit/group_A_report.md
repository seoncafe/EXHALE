# Audit A: hydrodynamics, time stepping, grid, gravity, boundaries, steady solver

Root: `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/`

---

### The cell width `dr_j(j)` is the width of cell j+1, not of cell j

Status: CONFIRMED
Severity: P0 (wrong answer in a default path)
Location: `src/modules/init/define_grid.f90:147-154` and the identical repeat at `:180-187`

```
      ! Cell edges (N+2*Ng-1 points) - r_edg(j) = r_{j+1/2}
      r_edg(1-Ng:N+Ng-1) = 0.5*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0*r_edg(N+Ng-1) - r_edg(N+Ng-2)

      !--- Cell dimensions r_{j+1/2} - r_{j-1/2} --- !
      ! Cell size (N+2*Ng points) - dr(j) = dimension of cell j
      dr_j(2-Ng:N+Ng-1) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)
      dr_j(1-Ng) = dr_j(2-Ng)
      dr_j(N+Ng) = dr_j(N+Ng-1)
```

**What the code does.** The two array sections are element-aligned so that the
element written into `dr_j(j)` is `r_edg(j+1) - r_edg(j)`. With the stated
convention `r_edg(j) = r_{j+1/2}` (the comment on the line above, and the
convention every consumer uses), cell j spans `[r_edg(j-1), r_edg(j)]`, so
`r_edg(j+1) - r_edg(j)` is the width of cell **j+1**. `dr_j(j)` is therefore
shifted outward by one cell everywhere.

**What it should do.** `dr_j(j) = r_edg(j) - r_edg(j-1)`, i.e. the LHS section
should start one index higher (`dr_j(3-Ng:N+Ng) = r_edg(3-Ng:N+Ng) -
r_edg(2-Ng:N+Ng-1)`), leaving `dr_j(2-Ng)` and `dr_j(1-Ng)` to the fallbacks.
The invariant is the one `RK_rhs` itself uses two lines apart: it forms
`dV = (r_edg(j)^3 - r_edg(j-1)^3)/3` for the SAME cell j whose `dr = dr_j(j)`
it takes.

**How it was verified.**
1. Fortran section semantics, checked by compiling the two statements verbatim
   on a stretched grid (`auditA/grid_probe.f90`): `dr_j(j)/(r_edg(j)-r_edg(j-1))`
   equals the stretch ratio q exactly at every cell.
2. The whole `Mixed` branch replicated verbatim with the production defaults
   (N=500, N_low_cells=50, dr_base=2e-4, r_max=10, Ng=2;
   `auditA/mixed_probe.f90`): `dr_j` is exact in the 50 uniform base cells and
   **1.4508% too large from cell 51 to cell 495**, the ratio being exactly the
   Newton-solved stretch parameter x0 = 1.014508. Sum over the physical cells:
   8.93798 against a true domain thickness of 8.81471.
3. On the ACTUAL grids of three regression cases, rebuilt from the `r` column of
   `backup/regression/golden/*/Hydro_ioniz.txt` (504 rows = `1-Ng..N+Ng`):

   | case | max \|dr_j/w - 1\| | median ratio (cells 100..N-5) | sum dr_j / sum true widths |
   |---|---|---|---|
   | `wasp_full` | 6.84e-3 | 1.006844 | 1.006663 |
   | `mol_base_handoff` | 1.22e-2 | 1.012185 | 1.011935 |
   | `lower_profile` | 1.17e-2 | 1.011718 | 1.011480 |

**Reaches.** `dr_j` has 12 consumers; the ones that change an answer in a
default path are:

* `src/modules/functions/utilities.f90:461-500` (`calc_column_dens`,
  `calc_column_dens_one`) and `:523-529` (`calc_column_dens_metals`): the
  column integral is the rectangle rule `Ncol(j) = Ncol(j+1) + n(j)*dr_j(j)*R0`,
  so **every photoionization optical depth is too large by the local stretch
  factor** (0.7% on `wasp_full`, 1.2% on `mol_base_handoff`, 1.5% on a
  r_max = 10 HD 209458 b grid). This is on by default for H I, He I, He II,
  He 2^3S, H2 and every photo-ionizable metal ion.
* `src/modules/states/Source.f90:34` under PLM: the gravity source
  `-rho_bar (phi_{j+1/2} - phi_{j-1/2})/dr` is divided by `dr_j(j)` while the
  pressure lives inside the flux and is divided by `dV`, so **the discrete
  hydrostatic balance itself is off**: the effective gravity is 0.7-1.5% too
  weak over the stretched region. PLM is stage 1 of the recommended two-stage
  marching recipe and the endpoint the PLM->WENO3 continuation starts from.
* `src/modules/time_step/RK_rhs.f90:96,119` under WENO3: `(pR - pL)/dr` and the
  gravity source share the same wrong `dr`, so the hydrostatic ZERO is
  unaffected; what is wrong is the weight of the pressure-plus-gravity pair
  relative to the advective momentum-flux divergence `(A_p F_p - A_m F_m)/dV`,
  which carries no `dr`. The steady momentum equation the WENO3 residual drives
  to zero is therefore `(1/r^2)d(r^2 rho v^2)/dr + (1/q)(dp/dr + rho g) = 0`
  with q = 1.007-1.015.
* `src/modules/time_step/eval_dt.f90:42,45`: the CFL step is too large by q.
* `src/modules/radiation/Cool_coeff.f90:2514` and
  `src/modules/radiation/util_ion_eq.f90:2534,2642`: the path length `dl` of
  the escape-probability / self-shielding terms.
* Volume weights `r^2 dr_j` in `steady_residual.f90:614` (the integrated
  residual norm) and `diffusive_photochemistry.f90:928`; the carrier row scale
  `steady_residual.f90:494`; `write_output.f90:434,664,730` (`dr_cm`); the
  ESWENO3 smoothness regularizer `Reconstruction.f90:104-105` (negligible,
  it is an epsilon).
* Not affected: `base_boundary.f90:451` (uses `dr_j(1-Ng)`, inside the uniform
  base region where `dr_j` is exact) and `write_setup_report.f90:314` (uses
  `dr_j(1)`, likewise exact).

Every regression case would move. The two definitions of the cell width now in
the code (`dr_j(j)` and `r_edg(j)-r_edg(j-1)`, both used inside `RK_rhs`)
disagree by q, which is exactly the duplicated-definition failure the project
convention warns about.

Provenance note: the statement is inherited verbatim from upstream ATES v2.0
(`ATES/ATES-Code-main/src/modules/init/define_grid.f90:122`). EXHALE is a
heavily modified fork, not a reference copy kept for diffing, so this is
reported rather than left as an external-code exception.

---

### A momentum source of 1e-16 per cell per step, added to avoid a division by zero in a diagnostic

Status: CONFIRMED
Severity: P1 (an unreported source term in the default path)
Location: `src/EXHALE_main.f90:1333`

```
      	u(2,:) = u(2,:) + 1.0e-16 ! To avoid division by zero
```

**What the code does.** Immediately before the `dtu` infinity-norm is formed as
`maxval(abs(1.0 - u(k,j)/u_old(k,j)))`, the momentum density of every cell
AND every ghost is raised by the absolute constant 1e-16 in code units. The
state is not restored afterwards: the next loop iteration takes `u_old = u` and
starts the RK stages from this `u`, so the kick is a persistent source of
1e-16 of momentum density per cell per step, and it is not paired with any
energy change (`u(3,:)` is untouched, so the internal energy implied by
`u(3) - u(2)^2/(2 u(1))` moves too).

**What it should do.** A guard for a diagnostic belongs in the diagnostic:
`abs(1 - u(2,j)/u_old(2,j))` should be formed with a denominator floor, or the
kick applied to a copy. Nothing in the conservation laws licenses adding
momentum to the gas.

**How it was verified.** Read the loop body end to end: the only writes to `u`
between line 1333 and the next `u_old = u` (line 999 of the following pass) are
`u = u_old` inside `retry_step`, which restores the kicked state, not a clean
one. `W` is formed at line 1269, before the kick, so the kick is invisible to
`eval_dt` and to every reported profile, and only enters through the state
itself. The constant is absolute in code units (`rho0 v0`), so its significance
scales inversely with `n0`: harmless where the wind carries `rho v ~ 2e-3`
(measured from the `wasp_full` golden: n0 = 3.09e12 cm^-3, rho/n0 ~ 6e-4,
v/v0 ~ 3.2, so the 14408-step total kick is 7e-10 of it), but on a tenuous
outer wind at `rho v ~ 1e-8` a 1e5-step run accumulates 1e-11, i.e. 1e-3 of the
local momentum. I did not run the code, so I have not measured the actual
outer-cell momentum of an r_max = 10 model.

**Reaches.** Every marching run. Also inherited verbatim from upstream ATES
(`ATES/ATES-Code-main/ATES_main.f90:293`).

---

### The Crank-Nicolson conduction stage clamps the temperature silently

Status: CONFIRMED
Severity: P2 (opt-in path, uncounted clamp)
Location: `src/modules/time_step/viscous_conduction.f90:527-532`

```
         do j = 1, N
            Tnew(j) = max(sol(j), 1.0d-2)
            W(3,j)  = n_part(j)*Tnew(j)
```

**What the code does.** The implicit temperature solve floors every cell at
0.01 T0 and reports nothing. The identical floor in the explicit-source energy
update (`energy_semi_implicit.f90:256-267`) is counted three ways
(`n_energy_floor_hits`, first/last step, cell histogram) precisely because,
as its own header says, "a cell that lands on the floor has NOT converged to a
physical temperature". The same statement is true here and the same run-end
report does not see it.

**What it should do.** Either count and report the activations as
`energy_semi_implicit` does, or let the state go non-positive and be caught by
the run's positivity detector.

**How it was verified.** Read both routines; `viscous_conduction` has no
counter module variable and no write of one. Also, unlike the explicit update,
the clamp here breaks the claimed fixed-point identity of the module header
("at the fixed point Tnew = Tcell the cap terms cancel"): a clamped cell is not
a zero of `viscous_conduction_sources`, so the marching stage and the steady
residual have different fixed points in exactly the cells that hit the floor.

**Reaches.** Only `Viscosity:` / `Conduction:` runs; both default off, so no
regression case moves.

---

### `minloc(..., dim=1)` result used as a declared subscript in the IC density probe

Status: CONFIRMED
Severity: P2
Location: `src/modules/init/set_IC.f90:116-118`

```
			r_half = 0.5e0*(r_max + 1.0e0)
			i_rhalf = minloc(abs(r-r_half), dim = 1)
			minrho = W(1,i_rhalf)
```

**What the code does.** `r` is declared `r(1-Ng:N+Ng)`. `MINLOC` returns the
position counted from 1, not the declared subscript, so `i_rhalf` is the
intended cell index plus `Ng`, and `W(1,i_rhalf)` samples a cell two positions
farther out than the one nearest `r_half`. When the nearest cell to `r_half`
lies within `Ng` of the top of the grid the read is out of bounds
(`i_rhalf` can reach `N+2*Ng` while the upper bound is `N+Ng`).

**What it should do.** `i_rhalf = minloc(abs(r-r_half), dim=1) + lbound(r,1) - 1`.

**How it was verified.** Compiled probe `auditA/minloc_probe.f90`: for
`r(-1:12)` and a target at `r(6)`, `minloc(..., dim=1)` returns 8 and `r(8)` is
two cells away from the target.

**Reaches.** Only the `b0_eff` ramp of the cold hydrostatic IC (it decides
whether the mid-domain density is above 1e-8), so the effect is a slightly
different `b0_eff` on some grids. No regression case is expected to move,
because the ramp exits on a coarse 1e-8 test.

---

### `RK_rhs` never writes `dF(:,1-Ng)` and `S(:,1-Ng)`, and both are read

Status: CONFIRMED
Severity: P2 (undefined-value read; overwritten before use)
Location: `src/modules/time_step/RK_rhs.f90:92-93` (`do j = 2-Ng,N+Ng`) against
`src/EXHALE_main.f90:1049-1051,1084-1086,1117-1119` and
`src/modules/time_step/steady_residual.f90:169-171`

**What the code does.** The face loop covers `1-Ng..N+Ng` but the cell loop
starts at `2-Ng`, so the `intent(out)` arrays `dF` and `S` are left at whatever
the caller's arrays held at index `1-Ng`. Three sites then read that element:
the RK stage updates `u1(k,:) = u(k,:) - dt_loc*(dF(k,:) - S(k,:))` (whole-array
expressions), `assemble_residual`'s `R(k,:) = dF(k,:) - S(k,:)`, and the
homotopy blend `dF = om*dF + lam*dFw` in
`reconstruction_continuation_rhs`. In every case the value is overwritten
(`Apply_BC`) or never read again (`R(:,1-Ng)` is outside every norm window), so
no answer moves, but a build with `-finit-real=snan` would trap on the first
step.

The same pattern is at `src/modules/states/Apply_BC.f90:231-232`
(`WL_out = WL_in`): neither `PLM_rec` (`do j = 2-Ng,N+Ng-1`) nor the ESWENO3
branch (`do j = 0,N+1`) ever writes `WL(:,1-Ng)`, and the whole-array copy reads
it before `Rec_BC` overwrites it three lines later.

**How it was verified.** Read the loop bounds of `PLM_rec`, the WENO3 branch,
`Rec_BC` and `RK_rhs`, and traced each of the four writers/readers of index
`1-Ng`. `positivity_limited_faces` runs after `Rec_BC`, so it never sees an
undefined element.

---

### `refresh_row_terms` rebuilds the row scales without the viscous source and with a monotone `max` on the energy row

Status: CONFIRMED
Severity: P2
Location: `src/modules/time_step/steady_residual.f90:298-303` against
`:250-257`

```
      do j = 1-Ng, N+Ng
         face_mass_flux_r2(j)     = face_flux(1,j)*r_edg(j)*r_edg(j)
         momentum_largest_term(j) = max(abs(dF(2,j)), abs(S(2,j)))
         energy_largest_term(j)   = max(abs(dF(3,j)), abs(S(3,j)),       &
                                        energy_largest_term(j))
      enddo
```

**What the code does.** Two differences from `store_row_terms`, which is the
routine this one is supposed to reproduce for a state nobody has evaluated:
(i) the operator-split viscous/conduction terms `Smom`, `Sene` are dropped from
the momentum and energy scales, so with `Viscosity:`/`Conduction:` on the two
paths give different scales for the same state; (ii) the energy scale is
`max(new hydro terms, PREVIOUS value)` rather than an overwrite, so across a
sequence of probe states it can only grow, and a residual measured after an
inflating probe is divided by a scale that belongs partly to a state the solve
rejected. `store_row_terms` overwrites, so an accepted state resets it.

**What it should do.** Either recompute the radiative terms, or state the scale
as belonging to the last stored state (as the header already does) without the
`max` accumulation, and carry `Smom`/`Sene` through the same way.

**How it was verified.** Read both routines and every caller of
`residual_row_scale` / `mass_flux_row_scale` / `momentum_row_scale` /
`energy_row_scale`. The header of `refresh_row_terms` documents the heat/cool
staleness but not the `max` accumulation and not the missing `Smom`.

**Reaches.** Only states that reach `refresh_row_terms` with
`state_of_row_terms /= u`, which the header measures as ~15 in 1e5 scale
requests. The one place it is systematic is the self-consistency loop at
`steady_newton.f90:2640-2647`: there `eval_residual` moves the composition (and
therefore `n_part_cell1`, and therefore the ghosts `Apply_BC` writes) on every
sweep while the caller's `u` keeps the ghosts written once at line 2587, so
`all(state_of_row_terms .eq. u)` fails and the row scales of the number the
gate finally reads are rebuilt on a ghost state one sweep old.

---

### The self-consistency sweep count is reported one too high when it does not converge

Status: CONFIRMED
Severity: P2 (report defect)
Location: `src/modules/time_step/steady_newton.f90:2640-2649`

`it_sc` is the DO variable of `do it_sc = 1, n_selfconsistent_max`; on normal
loop completion Fortran leaves it at `n_selfconsistent_max + 1`, and that is the
number printed as "`, it_sc, ' sweeps: ||R||='`". Only the non-converged case is
affected; a loop that exits early prints the right count.

---

### `j_min` has no lower clamp while `j_flux` does

Status: CONFIRMED
Severity: P2
Location: `src/modules/init/define_grid.f90:194-244`

The `j_flux` search ends with `j_flux = max(j_flux, 1)` and the `j_min` search
does not. With `Escape radius [R_p]:` at or below `r(1-Ng)` the DO loop exits on
its first trip and `j_min` is `1-Ng`; the wind window `[j_min:N]` then reaches
into the ghosts, where `residual_norms` reads `R(:,1-Ng)` (the undefined element
of the finding above), `resid_relnorm` reads a zero it filled itself, and the
marching loop's `mom_max/mom_min = maxval/minval(abs(mom(j_min:N)))` and the
`dtu` norms include ghost cells. The upper end is guarded and warned about
loudly; the lower end is not guarded at all. I did not find a place where
`r_esc` is validated against `r(1-Ng)` in `input_read.f90` (line 187 reads it
with no range test).

---

### The trust region's ray test uses two residuals whose admissibility it discards

Status: CONFIRMED
Severity: P2
Location: `src/modules/time_step/steady_newton.f90:1848-1860`

```
         Ytry = Y + D*(hray*s)
         call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,          &
                            admissible=okres)
         f2p = 0.5d0*sum((Ftry/Drow)**2)
         Ytry = Y - D*(hray*s)
         call eval_residual(Ytry, f_sp, fwork, Ftry, heat, cool,          &
                            admissible=okres)
         f2m = 0.5d0*sum((Ftry/Drow)**2)
         fslope = (f2p - f2m)/(2.0d0*hray)
```

`okres` is written twice and read nowhere. Everywhere else in this file an
inadmissible residual is refused rather than used (the header of
`eval_residual` states the rule explicitly), but here a probe whose chemistry
the sweep could not certify enters `fslope` and can reject a model that is
sound, shrinking the radius by 4x for a reason that is not the model's. The
backward probe `Y - D*(hray*s)` is also the one point in the solve that steps
AGAINST the descent direction, so it is the likeliest of the two to leave the
describable set. Reached only under `EXHALE_TRUST_REGION=1` with the carrier
unknown, so no default path and no regression case.

---

### The frozen-Jacobian finite-difference step keeps the flat `max(|Y|,1)` floor its replacement was written to remove

Status: CONFIRMED
Severity: P2 (diagnostic-only caller)
Location: `src/modules/time_step/steady_newton.f90:744` against `:788-791`

`build_banded_jac` uses `dYc(jcol) = sqeps*max(abs(Y(jcol)), 1.0d0)`, while
`build_banded_jac_full`, ten lines below, comments on exactly that expression:
"a flat max(|Y|,1) floor gives wind cells, ~1e-7 of the base in code units,
order-unity relative kicks and garbage columns", and uses
`cell_state_scales` instead. `build_banded_jac` has one caller,
`src/EXHALE_main.f90:674`, inside the `EXHALE_RESIDUAL` diagnostic block; the
production solves (`solve_steady_ptc:1271`, `solve_steady_jfnk:2130`) both call
the scaled version. So the columns the diagnostic prints are the ones the
comment calls garbage, and nothing says so at the diagnostic.

---

### The base boundary's reversal branch is wider than its comment says, and the ghost stratification uses the base composition on it

Status: CONFIRMED
Severity: P2 (documentation defect plus a small physical inconsistency in an
off-nominal branch)
Location: `src/modules/states/base_boundary.f90:346-349, 373-380, 446`

```
   w_rev = characteristic_branch_weight(M_i/base_face_mach_blend)
   ...
   rho_b   = (1.0d0 - w_rev)*rho_res + w_rev*rho_rev
   ...
   ! On the reversal branch the density came from the interior isentrope and
   ! cell 1's composition would be the consistent one; the difference is
   ! confined to |M_i| < base_face_mach_blend = 1e-6, ...
```

`characteristic_branch_weight` returns exactly 1 for `x <= -1`, so `w_rev = 1`
for every `M_i <= -1e-6`, i.e. over the WHOLE reversal branch and not only
inside the blend window. On that branch `rho_b` is the interior isentrope's
density while `c_b` (line 379) and, more importantly, `T_face` in
`base_ghost_averages` (line 446) are both formed with `base_reservoir_nhat` and
`jcomp = 0`, so the two ghost cells are stratified along the BASE
composition's caloric isentrope through a state whose density came from cell 1.
The comment defends this by a bound that does not hold. Severity is low because
the (rho, p) pair reproduces the face state exactly at `r_edg(0)` by
construction; what differs is the stratification across the two ghost cells.

---

### The Mixed-grid construction and one line of `Rec_BC` assume Ng = 2 without saying so

Status: CONFIRMED
Severity: P2
Location: `src/modules/init/define_grid.f90:77-83`; `src/modules/states/Apply_BC.f90:263-264`

`define_grid` seeds the Mixed grid with `r(1-Ng) = 1.0 - drc` and
`r(2-Ng) = 1.0` and then runs `do j = 1,N_low: r(j) = r(j-1) + drc`, which reads
`r(0)`. At Ng = 2 line 78 is exactly that seed; at Ng = 3 the loop reads an
`r(0)` nobody wrote and the whole grid is garbage, with no diagnostic.
`Rec_BC` likewise writes the literal `WL_out(:,N+2)`. `Ng` is a
`parameter = 2` in `parameters.f90:15`, so nothing is broken today; the lower
boundary of `Rec_BC` already carries a note about the Ng dependence and the
grid does not.

---

### HLLC returns an upwind face pressure where LLF and ROE return the mean

Status: SUSPECTED
Severity: P2
Location: `src/modules/flux/Num_Fluxes.f90:117,123,129,132` against `:232,297`

**What the code does.** The HLLC branch sets `p_out = pL` or `pR` according to
which wave region contains the face, while `lax_friedrichs_flux` and the ROE
branch both return `0.5*(pR + pL)`. `p_out` is consumed only under WENO3
(`RK_rhs.f90:119`, `dF(2,j) = dF(2,j) + (pR - pL)/dr`), where it IS the pressure
the momentum equation differences.

**Why the sum is nonetheless exact.** With `use_plm = .false.`, `Phys_flux`
omits `p` from the momentum component, so `NF(2) + p_out` reproduces the full
HLLC momentum flux exactly at every face. What the code does is split that flux
into a piece it differences with `(A_p . - A_m .)/dV` and a piece it
differences with `(. - .)/dr`; the split is exact but the two discrete
operators are not the same, so the choice of `p_out` is visible in the answer.

**What the face pressure should be.** The HLLC star pressure
`p* = pL + rhoL(SL - vL)(S* - vL)`, which at rest reduces to
`0.5(pL + pR) + O(jump^2)`, so the code's `pL` over-states the face pressure by
about half the reconstruction jump at that face. For a smooth WENO3 state the
jump is O(dr^3) and the resulting error in `(p_p - p_m)/dr` is O(dr^2), i.e. a
loss of one order in the momentum equation rather than an inconsistency.

**Why SUSPECTED.** I did not measure whether the asymmetry contributes
materially to the near-base momentum residual the `base_boundary` and
`low_mach_dissipation` headers are both about; closing that needs a run, and
this audit is read-only. `p*` is already computed on the way to `S_star` and
would cost nothing.

---

## Things checked and found sound

* The Roche potential and its radial derivative (`grav_field.f90`). Both match
  `Phi = -GM_p/r - GM_s/(a-r) - (G(M_p+M_s)/2a^3)(r - a M_s/(M_p+M_s))^2` with
  `b0 ~ GM_p`, `Mrapp = M_s/M_p`, `atilde = a/R_p`; `Dphi` is its exact
  derivative term by term. The `spherical_domain` branch drops the tidal and
  centrifugal terms consistently in both.
* The gravitational-work term in the energy flux (`RK_rhs.f90:121-123`). The
  expression is exactly the discrete identity
  `-rho v dphi/dr = -(1/dV)[A_p F_p (phi_p - phi_c) - A_m F_m (phi_m - phi_c)]`
  with `Gphi_i` at faces and `Gphi_c` at centres, which is why `S(3) = 0` in
  `Source.f90` is correct rather than a missing term. The duplicate in
  `positivity_limited_fluxes` matches it term for term.
* The Navier-Stokes stress and dissipation (`viscous_conduction.f90`). The
  header's `F_mu = (1/r^2)d(r^2 tau_rr)/dr + tau_rr/r` with
  `tau_tt = tau_pp = -tau_rr/2` and `q_mu = (4/3) mu (w' - w/r)^2` are the
  correct spherical Newtonian results, and `w F_mu + q_mu = div(tau.w)` closes;
  the tridiagonal coefficients reproduce them face by face, the base ghost
  column is correctly moved to the right-hand side of both Crank-Nicolson rows,
  and `mu = (4/15)(m_H/k_B) kappa` with the code-unit conversions
  `1/(n0 mu v0 R0)` and `T0/(n0 mu v0^3 R0)` is dimensionally right. The header
  is also right that the printed Koskinen+2022 (B5)/(B6) do not satisfy the
  energy identity.
* `eval_dt` carries no parabolic (diffusive) CFL limit, and does not need one:
  both transport operators are Crank-Nicolson. The low-Mach damping term IS
  explicit and its stability bound `eps4 < 1/(16 CFL)` is not enforced by
  `eval_dt`, but `input_read.f90:1114-1119` warns on it.
* The caloric EOS. `gamma_eff = 1 + 1/c_v` is simultaneously the energy index
  (`e = p/(gamma-1)`) and the isentropic index `c_p/c_v` for an ideal mixture
  with internal degrees of freedom, so using it in `sqrt(gamma p/rho)` is
  correct; the cubic Hermite table, its derivative, the low/high-T linear
  continuations and the safeguarded Newton inverse are all consistent.
  `continue_hydrostatic_isentrope`'s `d ln rho = c_v d ln T` and
  `isentropic_density_at_pressure`'s `d ln rho = d ln p / gamma` are the right
  relations for that EOS.
* The characteristic base boundary. The LODI relation
  `p_b - rho_i c_i v_b = p_i - rho_i c_i v_i` is the correct outgoing (v-c)
  invariant at a lower boundary; `v_b = v_i + (p_b - p_i)/(rho_i c_i)` follows;
  the weight-1 form quoted in the header is the correct simultaneous solution of
  that pair; `characteristic_branch_weight` is a C1 cubic smoothstep with the
  right endpoint values and zero endpoint slopes; the 8-point Gauss-Legendre
  nodes and weights are correct.
* The transonic IC. The Bernoulli integral
  `v^2/2 - c2 ln v - 2 c2 ln r + phi = B` and the critical condition
  `dphi/dr = 2 c2/r` are the correct isothermal-wind relations; the branch
  selection and the two bisections are sound.
* The HLLC star states (`usL`, `usR`) are Toro's, with `phi = rho(S-v)` used
  correctly in `p/phi`. The ROE wave-strength coefficients
  `a1 = (dp - rho a dv)/(2a^2)`, `a2 = drho - dp/a^2`,
  `a3 = (dp + rho a dv)/(2a^2)` are Toro eq. 11.68.
* Positivity plumbing: `positivity_limited_faces` and
  `positivity_limited_fluxes` index consistently, the `1-Ng..N+Ng` bounds of the
  ESWENO3 coefficient arrays `C1/C2/D1/D2/dV` are all covered at Ng = 2, and the
  four-cell stencil of `contact_mode_dissipation_flux` stays inside the
  `0:N+1` arrays it builds.
* Determinism of the OpenMP regions in this group. `RK_rhs` (two loops, one
  barrier, `face_flux`/`face_p` written by index and read by index), `PLM_rec`,
  the ESWENO3 loop and `U_to_W`/`W_to_U` are all cell-local with
  `schedule(static)` and no reductions, so the bitwise claims in their comments
  hold. `S0sav`/`S1sav` are allocated before the parallel region and written by
  owning index only.
* `Minmod_MC` is called with rank-1 actuals against `dimension(1,3)` dummies;
  this is legal sequence association (explicit-shape dummy) and the element
  order matches, so it is ugly rather than wrong.
* `assemble_residual` contains no time term and no `dt = 1e30` device; the
  residual is `dF - S - (heat - cool)` directly. Nothing in the group sets a
  large pseudo-`dt` to remove a time derivative.
* `define_grid`'s `real*8 :: tol = 1.0` is implicitly SAVEd, which would break
  the stretch-parameter Newton on a second call; `call define_grid` occurs
  exactly once (`init.f90:99`), so it is latent, not live.

## Files read completely

`src/modules/flux/Num_Fluxes.f90`, `src/modules/flux/speed_estimate_ROE.f90`,
`src/modules/flux/low_mach_dissipation.f90`,
`src/modules/states/Reconstruction.f90`, `src/modules/states/PLM_rec.f90`,
`src/modules/states/Source.f90`, `src/modules/states/Apply_BC.f90`,
`src/modules/states/base_boundary.f90`, `src/modules/states/caloric_eos.f90`,
`src/modules/time_step/RK_rhs.f90`, `src/modules/time_step/eval_dt.f90`,
`src/modules/time_step/energy_semi_implicit.f90`,
`src/modules/time_step/viscous_conduction.f90`,
`src/modules/time_step/steady_residual.f90`,
`src/modules/time_step/steady_newton.f90`, `src/modules/init/define_grid.f90`,
`src/modules/init/set_gravity_grid.f90`, `src/modules/init/init.f90`,
`src/modules/init/set_IC.f90`, `src/modules/functions/grav_field.f90`,
`src/modules/functions/UW_conversions.f90`, and `src/EXHALE_main.f90` lines
900-1470.

Read in part, for cross-checking only:
`src/modules/functions/utilities.f90:450-535` (the `dr_j` consumers),
`src/modules/files_IO/input_read.f90` (targeted lines for `r_esc`, `n0`, `v0`,
the low-Mach warning), `src/modules/init/parameters.f90` (targeted lines),
`src/modules/files_IO/write_setup_report.f90:300-320`,
`ATES/ATES-Code-main/src/modules/init/define_grid.f90` and
`ATES/ATES-Code-main/ATES_main.f90` (provenance of two findings).

## Not checked

* I did not run EXHALE, `make`, or the regression matrix. Every number above
  comes either from a standalone probe I compiled in the scratch directory
  (`grid_probe.f90`, `mixed_probe.f90`, `minloc_probe.f90`) or from arithmetic
  on the committed goldens' `r` column.
* The magnitude of the `dr_j` error on the mass-loss rate, the ionization front
  position and the heating profile: that needs a run with the corrected width,
  which is out of scope for a read-only audit.
* `IC_load` / `load_IC` and the restart interpolation of composition. The brief
  listed it, but the file is not in my group's list and I did not open it;
  `init.f90` only calls it. Whoever covers `files_IO` should take it.
* `wae_exhale_bridge`, `element_diffusion_step`,
  `photochemical_transport_step`, `excited_H_update`, `ioniz_eq`,
  `get_species_densities`, `comp_T_from_p` / `comp_p_from_T`: called from the
  marching loop but owned by other groups. I checked only their call ORDER
  against the residual's (`eval_residual` reproduces it: state, composition,
  ghosts, fluxes) and did not read them.
* `EXHALE_main.f90` outside lines 900-1470, including the diagnostic blocks at
  530-700 that are the only callers of `build_banded_jac`; I read the call
  sites but not the blocks around them.
* `binary_element_diffusion.f90` and `diffusive_photochemistry.f90` use
  `Dphi(r(j))*v0*v0/R0` as a physical `g`; the unit conversion is right
  (`phi` is in `v0^2`, `r` in `R0`), but those files belong to another group and
  I read only those two lines.

# Explicit viscosity and heat conduction

Molecular transport in the bulk (single-fluid) radial equations: the viscous
momentum force, its dissipation, and thermal conduction. This is the
Navier-Stokes level that CETIMB (Koskinen et al. 2013a, Icarus 226, 1678;
Koskinen et al. 2022, ApJ 929, 52) carries and that EXHALE's inviscid HLLC
scheme lacks. The motivation is item (A) of `TO_BE_DONE.md`: the near-base
momentum residual that floors the JFNK steady solve at `||R|| ~ 1e-3` in
full-physics runs.

Implementation: `src/modules/time_step/viscous_conduction.f90`.
Both switches default OFF.

**2026-08-10 (later the same day): the motivation stated above no longer
holds.** There is no near-base momentum wall. The `info = 2` residual floor
reported throughout Sec. 5 and Sec. 8 was traced to the JFNK diagonal scaling
and the stagnation watchdog, and the same cases now converge with `info = 0`;
see `docs/newton_scaling_and_base_wall.md`. The derivation, the discretization,
the verification of Sec. 7 and the physical A/B numbers of Sec. 8 are
unaffected — only the diagnosis they were measured against is. All `info = 2`
outcomes and every "worst residual cell `j = 1`" entry below were produced with
the old scaling and the old worst-cell print, both of which normalized momentum
by the base cell's `|rho v|`.

---

## 1. Equations

Koskinen et al. (2022) Appendix B writes the bulk momentum and energy
equations as their Eqs. (B2) and (B3),

```
d(rho w)/dt + (1/r^2) d/dr(r^2 rho w^2 + p) = -rho dU/dr + F_mu           (B2)
d(rho u)/dt + (1/r^2) d/dr(r^2 rho u w)
     = rho q - p (1/r^2) d/dr(r^2 w) + (1/r^2) d/dr(r^2 kappa dT/dr) + q_mu (B3)
```

with `u = c_v T`, and gives the viscous momentum force and the dissipation
functional as their Eqs. (B5) and (B6). As printed those read

```
F_mu = (4/3)(1/r^2) d/dr(r^2 mu dw/dr) - (dmu/dr)(dw/dr) - (16/3) mu w/r^2  (B5)
q_mu = (10/3) mu (dw/dr)^2 - (8/3) mu (w/r)(dw/dr) + (16/3) mu (w/r)^2      (B6)
```

**Neither expression is the Navier-Stokes result for a radial flow, and we do
not implement them as printed.** For the Newtonian stress
`tau_ij = mu (d_i w_j + d_j w_i) - (2/3) mu delta_ij div w` in spherical
symmetry the only non-zero components are

```
tau_rr = (4/3) mu (dw/dr - w/r),   tau_tt = tau_pp = -(1/2) tau_rr
```

so that, exactly,

```
F_mu = (div tau)_r = (1/r^2) d/dr(r^2 tau_rr) + tau_rr/r
     = (4/3)(1/r^2) d/dr(r^2 mu dw/dr) - (4/3)(dmu/dr)(w/r) - (8/3) mu w/r^2   (1)
q_mu = tau : grad w = (4/3) mu (dw/dr - w/r)^2                                 (2)
```

Comparison with the printed appendix:

| term | Koskinen (B5)/(B6) as printed | Navier-Stokes, Eqs. (1)/(2) |
|---|---|---|
| leading diffusion | `(4/3)(1/r^2) d/dr(r^2 mu dw/dr)` | same |
| viscosity-gradient | `-(dmu/dr)(dw/dr)` | `-(4/3)(dmu/dr)(w/r)` |
| geometric | `-(16/3) mu w/r^2` | `-(8/3) mu w/r^2` |
| dissipation, `(dw/dr)^2` | `10/3` | `4/3` |
| dissipation, `(w/r)(dw/dr)` | `-8/3` | `-8/3` |
| dissipation, `(w/r)^2` | `16/3` | `4/3` |

Two independent checks say the printed pair cannot be right. (2) is a perfect
square, hence non-negative, as a dissipation must be; (B6) is positive definite
but is not `tau : grad w` for any `tau`. And (1) and (2) are exact partners,
`w F_mu + q_mu = div(tau . w)`, so that the pair conserves total energy; the
printed pair does not satisfy that identity. The same appendix also prints the
pressure inside the spherical divergence in (B2), which likewise cannot be meant
literally (it would add a spurious `2p/r`). We read the appendix as loosely
transcribed rather than as a different physical model, and implement (1)-(2).
*This is an interpretation; it has not been checked against the CETIMB source.*

EXHALE evolves the **total** energy `E = rho w^2/2 + p/(gamma-1)`, so the
energy source that belongs with the momentum source `F_mu` is

```
S_E = w F_mu + q_mu + (1/r^2) d/dr(r^2 kappa dT/dr)                            (3)
```

(`w F_mu` is the viscous work that converts kinetic energy; `q_mu` is the
irreversible part that heats the gas.)

---

## 2. Transport coefficients

**Neither Koskinen et al. (2013a) nor (2022) prints `mu(T)` or `kappa(T)`.**
O'Neill & Chorlton (1989), *Viscous and Compressible Fluid Dynamics* (Ellis
Horwood), is cited in Koskinen et al. (2013a) Sec. 2 only as the textbook
source of the dissipation functional, not as a viscosity parameterization. So
there is no coefficient in those papers to reproduce, and we take the standard
atmospheric-escape choice instead.

**Heat conduction** — atomic hydrogen, Watson, Donahue & Walker (1981, Icarus
48, 150); quoted in this closed form as Eq. (A6) of Erkaev et al. (2016, MNRAS
460, 1300):

```
kappa(T) = 4.45e4 (T/1000 K)^0.7      erg cm^-1 s^-1 K^-1                      (4)
```

**Viscosity** — tied to (4) by the Chapman-Enskog/Eucken relation for a
monatomic gas, `kappa = (15/4)(k_B/m) mu` (Prandtl number 2/3):

```
mu(T) = (4/15)(m_H/k_B) kappa(T) = 1.44e-4 (T/1000 K)^0.7   g cm^-1 s^-1       (5)
```

The code evaluates (5) from (4) at run time rather than hard-coding
`1.44e-4`, so the two stay consistent by construction.

### Validity

(4) and (5) are **neutral atomic-hydrogen** values. Above the ionization front
the electron (Spitzer) conductivity `~ T^{5/2}` is far larger, and Coulomb
collisions raise the ion viscosity; that regime is *not* covered. The terms are
included for the dense, largely neutral base (`r <~ 1.1 R_p`), which is where
the momentum imbalance they are meant to damp lives. Koskinen et al. (2013a)
state (their footnote 4) that conduction and viscosity are not important in the
thermosphere of HD 209458 b, so the omission is expected to be numerically
small there — but it is an omission, and it is recorded here rather than hidden.

### Code units

Lengths in `R0`, velocities in `v0 = sqrt(k T0/m_H)`, densities in
`rho0 = n0 m_H`, temperatures in `T0`; a force density then scales with
`rho0 v0^2/R0` and an energy-density rate with `q0 = rho0 v0^3/R0`.
Substituting into (1)-(3) leaves the expressions form-invariant with

```
mu_code    = mu_cgs    /(rho0 v0 R0)             (= 1/Reynolds number)
kappa_code = kappa_cgs T0/(rho0 v0^3 R0) = (15/4) mu_code
```

---

## 3. Discretization

Finite volume on the existing grid, with the same face areas and cell volume
`RK_rhs` uses: `A_{j+1/2} = r_edg(j)^2`, `dV_j = (r_edg(j)^3 - r_edg(j-1)^3)/3`.

```
F_mu(j) = (A_p tau_p - A_m tau_m)/dV_j + tau_c(j)/r(j)
Q(j)    = (A_p kappa_p (T_{j+1}-T_j)/dr_p - A_m kappa_m (T_j-T_{j-1})/dr_m)/dV_j
q_mu(j) = (4/3) mu_j D_j^2 ,   D_j = dw/dr|_j - w_j/r_j
```

Face values of `mu` and `kappa` are linearly interpolated with the geometric
weight `(r_edg(j) - r_j)/(r_{j+1} - r_j)` (the grid is non-uniform). `D_j`
uses the centered difference `(w_{j+1}-w_{j-1})/(r_{j+1}-r_{j-1})`, one-sided in
the outermost cell.

Both operators are assembled **once** as tridiagonal coefficient triplets
(`viscous_momentum_coeffs`, `thermal_conduction_coeffs`) and those same
triplets are used for the residual source *and* as the matrix of the implicit
update — see Sec. 5.

### Boundary conditions

**Inner (base) face.** The flux is formed from ghost cell `j = 0`, i.e. from the
state `Apply_BC` anchors there: `rho_bc`, the valved or mass-flux base velocity,
and the BC pressure. During an implicit step that ghost is Dirichlet data and
its column moves to the right-hand side. The base cell therefore feels a viscous
stress and a heat flux *relative to its anchor* rather than only exchanging with
the cell above it, which is the entire point of the term.

Note what the existing BCs make of this in practice. With the legacy one-way
valve the ghost velocity is `max(w_1, 0)`, so for an outflowing base
`w_0 = w_1` and the *gradient* part of the base-face stress vanishes; what
survives is the geometric `-w/r` part of `tau_rr` and the stress against cell 2.
With `Base velocity: massflux` the ghost velocity is `F_c/(rho_bc r^2)`,
independent of `w_1`, and the base face becomes a genuine Dirichlet anchor. The
*thermal* anchor is unconditional: the legacy ghost pins `p = ntot_bc + dp_bc`
and `rho = rho_bc`, i.e. `T_0 = T0`, so conduction ties the base cell
temperature to the base temperature in every configuration.

**Outer face.** Zero diffusive flux. The upper ghosts are a zero-gradient copy
(PLM) or a linear extrapolation (WENO3) of the interior; the extrapolated ghost
carries a gradient that would inject an unphysical viscous stress and heat flux
into the outflow. A vanishing diffusive flux is the correct condition at a
supersonic outflow boundary, and it also removes the upper-ghost column from the
tridiagonal system entirely.

---

## 4. Time integration

Viscosity and conduction are diffusive: an explicit update is limited by
`dt < rho dr^2/mu` (and `dt < C dr^2/kappa`), which in the finely gridded base —
exactly where the terms matter — is far below the advective CFL step. They are
therefore integrated **semi-implicitly (Crank-Nicolson)** as an operator-split
stage of the marching loop, as CETIMB does:

```
(i)  rho (w* - w)/dt = (1/2)[F_mu(w) + F_mu(w*)]
(ii) C   (T* - T)/dt = (1/2)[Q(T) + Q(T*)] + q_mu(w*),   C = n_part/(gamma-1)
```

Each is a tridiagonal solve (Thomas algorithm). Stage (i) rebuilds the
conserved state at **fixed pressure**: the kinetic-energy change is then
`d(rho w^2/2) = w_avg rho (w* - w) = dt w_avg F_mu_avg`, which is precisely the
viscous work, so stage (i) leaves the internal energy untouched and the
`w F_mu` term of (3) is accounted for exactly. Stage (ii) adds `q_mu` and the
conduction to the internal energy.

At the resolutions tested these terms turn out *not* to be stiff: the explicit
diffusive limit at the WASP-121 b base is `~C dr^2/(2 kappa) ~ 1 t_s` against an
advective step of `~9e-5 t_s`, four orders of magnitude of margin. The implicit
treatment is therefore unconditional insurance rather than a necessity here; it
matters for a denser base, a coarser grid, or a larger coefficient.

The recursion runs along `r`, so the solves are inherently serial in `j`; they
are called from the serial part of the marching loop and are never placed inside
a cell-parallel region.

---

## 5. Why the marching loop and the Newton residual solve the same system

JFNK must solve *exactly* the system the marching relaxes: a term present in one
and absent from the other makes the two paths converge to different states.
(This section originally cited the `info=2` stagnation of
`docs/base_composition_handoff_plan.md` Sec. 11.8-11.9 as evidence for that
requirement. That attribution is withdrawn — the stagnation was the solver's
scaling and watchdog, Sec. 11.10 there — but the consistency requirement itself
stands on its own and is what the three properties below enforce.)

1. **One definition of the source.** `viscous_conduction_sources` returns the
   pair `(F_mu, w F_mu + q_mu + Q)`. `assemble_residual` subtracts exactly that
   pair from `R(2,:)` and `R(3,:)`; `viscous_conduction_step` relaxes exactly
   that pair. There is no second expression anywhere.
2. **One operator.** The residual evaluation and the implicit matrix are built
   from the same tridiagonal triplets, so they cannot differ in a face weight,
   an interpolation, or a boundary term.
3. **The split fixed point is the residual zero.** Over one step the conserved
   state changes by `-dt(dF - S)` (RK), `+dt(heat - cool)` (energy stage) and
   `+dt(F_mu, w F_mu + q_mu + Q)` (transport stage). A fixed point of the
   composite requires the sum to vanish, i.e.
   `dF - S - F_mu = 0` and `dF_E - S_E - (heat-cool) - (w F_mu + q_mu + Q) = 0`
   — which is the residual `assemble_residual` returns.

One consequence for the interface: `assemble_residual` now takes
`n_part = n_tot + n_e` instead of nothing, and forms `T = p/n_part` internally.
Passing the *particle count* rather than `T` keeps `T` a function of the Newton
unknowns, so the temperature dependence of the conduction operator is captured
by the residual's linearization — including in the frozen-radiation residual
that builds the banded preconditioner.

---

## 6. Input keys

| key | effect | default |
|---|---|---|
| `Viscosity: True` | calibrated `mu(T)`, Eq. (5), plus the dissipation `q_mu` | off |
| `Viscosity: <mu0> [<s>]` | diagnostic power-law override `mu = mu0 T^s` in **code** units | off (`mu0 = 0`) |
| `Conduction: True` | heat conduction with `kappa(T)`, Eq. (4) | off |

Both are one-word keys, so the value is word 2 of the line (as for `Solver:` and
`CFL:`). The pre-existing numeric `Viscosity:` branch read word 3, which no
input file ever exercised; the position is corrected here.

The resolved settings are echoed at startup by `input_read` and appear in
`EXHALE_setup.out` and in the `EXHALE_PARSE_DUMP=1` dump as `visc_on`,
`cond_on`, `visc_mu0`, `visc_s`. The two new dump lines make
`backup/regression/parse_golden/` stale; it has deliberately not been
re-snapshotted here (re-snapshotted 2026-08-15 together with the CODATA `k_B`
refresh; the snapshots now carry the `visc_on` / `cond_on` /
`visc_mu0` / `visc_s` lines).

---

## 7. Verification of the discretization

Four analytic checks on a geometrically stretched 500-cell grid
(`opcheck.f90`, run against the compiled module):

| check | expectation | measured |
|---|---|---|
| `w = c r` (uniform expansion), variable `mu(T)` | `D = 0`, so `F_mu = 0`, `q_mu = 0` | `max|F_mu| = 2.2e-16` (scale 0.37), `max|q_mu| = 6e-34` |
| `w = 1/r^2`, constant `mu` | `F_mu = 0` analytically | `max|F_mu|/(mu/r^4) = 3.3e-4` (2nd-order truncation) |
| `kappa = k0 T^0.7`, `T^1.7 = a + b/r` | conduction source `= 0` | `max|Q| = 2.3e-13` (flux scale 1.5e-7) |
| `S_E` vs `w F_mu + q_mu + Q` | identical by construction | difference `0`; `min q_mu = +7e-21 >= 0` |

The first row is the strongest: the discrete operator annihilates the exact
null space of the continuum operator, with a spatially varying viscosity, to
machine precision. That exercises the face interpolation, the geometric
`tau_rr/r` term and the strain-rate stencil simultaneously.

---

## 8. Gate results (2026-08-10)

**G1 — default-off byte identity.** `backup/regression/run_check.sh check` over
`wasp_full`, `wasp_he23off`, `mol_base_handoff`. All three carry no
`Viscosity`/`Conduction` key, all three reproduce their goldens bitwise.

**G2 — acceptance run** (`wasp_full` + `Solver: Newton`, `Viscosity: True`,
`Conduction: True`, no Shapiro filter, secondary ionization at its default
staging), against a control that is the identical binary and configuration
*without* the two keys:

| | transport ON | control (OFF) | flux-converged golden |
|---|---|---|---|
| JFNK | `info=2`, best `||R|| = 1.195` | `info=2`, best `||R|| = 1.201` | — |
| worst residual cell | `j=1, r=1.000, k=2` | `j=1, r=1.000, k=2` | — |
| `log10 Mdot` | 13.20 | 13.20 | 13.22 |
| negative densities | 20 entries / 5 cells, r = 1.0145-1.0154 | 20 / 5, same cells | 68 / 17, r = 1.0106-1.0140 |
| converged `T`, `rho` | — | — | on-vs-off differ by < 0.05% |

So: **the `||R|| < 1e-3` target is not met, and the terms change nothing.**
Mdot is 4.5% below the golden, inside the quoted 2-5% band, and identical to
the control — the shift is the known path dependence of the `du` stop, not the
new physics. The negative densities are the same cells, count and magnitude as
the control (the most negative entry differs in the 4th digit), so they are
untouched by molecular transport.

### Why viscosity cannot meet the target — the measurement

Evaluated on the flux-converged golden WASP-121 b state (`EXHALE_RESIDUAL=1`
with and without the keys, so the difference *is* the new source):

| quantity | value |
|---|---|
| volume-weighted `sum |R_mom| V` | 1.905e-3 |
| volume-weighted `sum |F_mu| V` | 1.485e-6 (0.078%) |
| worst cell `j=1`: `|R_mom|` | 1.16 |
| worst cell `j=1`: `|F_mu|` | 2.06e-4 (0.018%) |
| `mu` multiplier needed to cancel `R_mom` at `j=1` | **5.7e3** |
| volume-weighted energy source / `|R_ene|` | 0.16% |
| transport energy source / local radiative heating (base cells) | 0.3-1.3% |

A direct test with `Viscosity: 4.2e-5 0.7` (that 5.7e3 factor, i.e.
`mu ~ 0.8 g cm^-1 s^-1`) drives the `j=1` momentum residual from 1.16 to
4.1e-3 — and moves the imbalance to `j=2` (0.77), with the volume-weighted
norm rising from 0.107 to 3.68. So even an unphysically large viscosity does
not remove the base imbalance; it relocates it.

The same conclusion in cgs: at the WASP-121 b base
(`n ~ 3.1e12 cm^-3`, `dr ~ 3e6 cm`, `mu ~ 2.6e-4 poise`, `rho g ~ 4.3e-9
dyn cm^-3`) the viscous force `mu w/dr^2` reaches `rho g` only at
`w ~ 1.5e8 cm s^-1`, four orders above the actual base velocity. The cell
Reynolds number at the base is ~1e8.

*Interpretation.* The numbers above establish only that the viscous force is
too small by ~4 orders to cancel the base-cell momentum residual, and that
raising `mu` enough to cancel it at `j = 1` relocates the imbalance to `j = 2`.
The further reading recorded here originally — that the residual is therefore a
property of the lower boundary condition — was **withdrawn on 2026-08-10**:
the base momentum row was subsequently measured to be satisfiable (a 2.9 ppm
ghost-pressure change nulls it), and the `j = 1` "worst cell" was an artifact of
a diagnostic that normalized momentum by the base cell's own `|rho v|`. See
`docs/newton_scaling_and_base_wall.md`.

**G3 — the full-physics continuation**
(`vulcan_work/hd209_wind_response/photo_deep_secion_cont/`, HD 209458 b with a
`base.inp` handoff and a molecular base), re-run in a scratch copy with both
keys added:

| | transport ON | the 2026-08-09 baseline (no transport) |
|---|---|---|
| JFNK | `info=2`, best `||R|| = 1.662e-2` | `info=2`, best `||R|| = 2.802e-3` |
| final `du` | 0.117 at 28990 steps | 0.103 at 41375 steps |
| `log10 Mdot` | 9.67 | 9.66 |
| negative densities | none | none |

Still no Newton-converged full-physics state; against the secondary-ionization-off
reference of 9.703 both are ~7% low in Mdot. The JFNK *best* residual is worse
with the terms on, and the reason appears to be the hand-off timing: the JFNK
fires 2000 steps after the continuation starts, while the conduction operator's
relaxation time at the base is `C dr^2/kappa ~ 2 t_s`, i.e. ~2e4 base steps.
Unlike WASP-121 b, HD 209458 b's cold molecular base *is* sensitive to
conduction — base `T` moves by up to 4% below r = 1.02 and `rho` by up to 13%
near r = 1.07 — so the hand-off state is still adjusting. (Tentative; not
demonstrated by re-running with a later hand-off.)

**Marching/residual consistency, measured.** On that same transport-relaxed G3
state the volume-weighted residual is **3.61e-2 with the transport terms
included and 5.06e-2 with them excluded**. The state the marching produced is a
better steady state of the system *with* the terms than of the inviscid system,
which is the observable signature that the two paths solve the same equations.

**G4 — comparison with Koskinen et al. (2013a)** HD 209458 b. Only qualitative:
their C2 model uses a 1 microbar lower boundary at `T0 = 1300 K`, the average
solar spectrum (0.45 W m^-2 shortward of 912 A after the global-average
division by 4), and multi-species diffusive separation, none of which the
EXHALE benchmark reproduces. Their C2 peak temperature is 12,000 K near
1.5 R_p and their mass-loss rate is `4.1e7 kg s^-1` (`log10 Mdot = 10.61` in
g s^-1), with the sonic point above 5 R_p. The EXHALE HD 209458 b benchmark
(no transport) peaks at 6875 K at r = 1.51 R_p with
`log10 Mdot = 10.89`: the temperature-peak *radius* agrees, the peak
temperature is about half, and the mass-loss rate is a factor ~1.9 higher.
The transport terms do not move any of this (Sec. 8, G2 A/B).

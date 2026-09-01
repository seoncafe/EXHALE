# A molecular base that runs supersonic near He/H ~ 1-10: what the lower boundary condition imposes, and what it does not

Written 2026-08-31 as the standalone record of `TO_BE_DONE.md` item (P).
Background: `docs/Update_EXHALE.md` sections 82 (the third body), 93 (row
scaling and the NaN abort matrix), 95 (reconstruction positivity), 111 (ghost
positivity and the local flux correction); `docs/molecular_chemistry_audit_he_rich.md`;
`docs/EXHALE_BC_and_IC.tex` is the reference statement of the boundary
conditions and this document does not restate it, only the parts item (P)
turns on.

Everything in sections 2, 3.1, 3.3, 5 and 6 below was measured or derived
for this document from the code and from stored output. The numbers carried
over from item (P) -- the He/H = 1 profile rows and the section-93 abort
matrix -- are marked where they appear, because **the output directories of
the He/H = 1 and He/H = 10 molecular arms were not kept** and they cannot be
re-measured without re-running.

**Follow-up status, 2026-08-31.** Section 11 records new controlled runs with
the current executable and supersedes the earlier statement that the molecular
chemistry solver was excluded. The base failure is initiated by a molecular
equilibrium state that satisfies the element bounds but is not required to
satisfy the reaction residual. In several runs the first physical cell is
reported as fully ionized hydrogen at 28--244 K. The resulting recombination
and free-free cooling, not H3+ cooling, collapses the thermal state. The lower
boundary then amplifies or masks that failure depending on its velocity and
pressure closure; it is not the initiating defect.

---

## The judgement first

**The lower boundary condition of EXHALE imposes two of the three primitive
variables and lets the third float, and that is the correct count only while
the base flow is subsonic.** `BC_component_constrho` pins `rho = rho_bc` and
`p = ntot_bc + dp_bc`, and sets the ghost velocity to `max(v_1, 0)` -- a copy
of the first interior cell. For gas entering the domain from below at
`0 < v < c` the boundary has two incoming characteristics (`v`, `v+c`) and one
outgoing (`v-c`), so exactly two conditions must be imposed and one component
must be left free: the code's choice is the well-posed one. Once `v > c` all
three characteristics enter the domain, the outgoing path through which the
interior state used to reach the boundary is gone, and the copied velocity is
no longer a boundary condition at all -- it is a feedback with no restoring
term. **Nothing in the code tests for that crossing.** This is the answer to
item (P)'s question "should the base BC be able to reach a supersonic inflow
at all": it should not, and it is not prevented from doing so.

What the original investigation did not settle was why the crossing happens in
a window near He/H ~ 1-10 and nowhere else. The follow-up in section 11 changes
one of its exclusions:

* **the chemistry solver is not excluded.** Section 93 showed that row scaling
  changed the MINPACK path without moving the old abort. It did not check the
  reaction residual of the state retained after a nonconverged solve or after
  projection onto the element bounds. The controlled runs in section 11 expose
  that unchecked path in the first physical cell;
* update positivity (sections 95 and 111: the same arm now completes 12000
  steps with 0 retaken steps and 0 dt bisections, and still marches
  `du = 13.5`, so positivity was never what made the state supersonic);
* **the base mean molecular weight on its own**, which item (P) named as a
  candidate ("a mean molecular weight that has moved by 5x"). Computed from
  the code's own base composition functions, the base sound speed that the
  boundary condition fixes moves from 2.71e5 to 1.98e5 cm/s over He/H =
  0.0793 to 1000 -- a factor 1.37, monotonic and saturating. It cannot make a
  window at 1-10, and it puts the sonic threshold at essentially the same
  velocity at every composition tested.

The initiating defect is therefore no longer open. The exact constrained
chemical root and the width of the composition window still need to be
computed after the solver is corrected. Section 11 gives the required solver,
boundary-composition, and validation changes.

**Scope, measured and sharper than item (P) recorded it.** No configuration
in production is inside the window, but not because the compositions avoid it:
646 of the 942 `input.inp` files under `LHS1140b/` sit at He/H between 1 and
10, squarely inside it. They are safe because **none of the 942 has molecular
chemistry on**. Every `input.inp` in the repository that does have it on (six:
two under `docs/lower_atmosphere_figs/`, four under `examples/`) is at
He/H <= 0.0833. The exposure is therefore a switch, not a composition: turning
molecular chemistry on for the He-rich LHS 1140 b ladder -- a natural thing to
want for a sub-Neptune -- lands directly in the window.

**A side finding recorded in section 7**: the base row of item (P)'s own
profile table is internally inconsistent. Its pressure, 19.674 erg/cm^3, is
not a value the base pressure BC can produce at that composition (the pinned
value is `ntot_bc` x `p0` = 11.805), and the particle count per nucleus it
implies, 1.500, requires every hydrogen nucleus to be ionized at 950 K in a
shielded base. The deeper rows of the same table are consistent (0.750 =
fully molecular, and 1.000 = atomic and neutral in the shocked cell). The
supersonic conclusion does not depend on it -- the recorded base velocity of
1.08e7 cm/s is 50 times the sound speed the BC pins, whichever reading is
taken -- but the snapshot should be re-measured and kept.

---

## 1. The symptom, row by row

Item (P) records this profile, from `backup/regression/mol_diffusion` re-run
at `He/H number ratio: 1.0` with `HeH_base 1.0`, at step 3405 of a 12000-step
relaxation. **Carried over from item (P); not re-measured (the run directory
was not kept).**

```
 idx     r        rho[mH/cm3]      v[cm/s]      p[cgs]       T[K]      Mach
   1   1.00000   2.5000e+14   1.0808e+07   1.9674e+01     950.00    38.60
   2   1.00019   2.4993e+14   1.0808e+07   4.0382e+00     195.05    85.18
   3   1.00039   2.4989e+14   1.0807e+07   7.8128e+00     754.84    61.23
 ...
 233   1.13895   2.9009e+14   8.2214e+06   5.5775e+00     464.20    59.40
 234   1.14072   6.1568e+14   7.9549e+06   7.7944e-01      30.56   223.98
 235   1.14251   3.9703e+15   6.0585e+06   1.8708e+03    8531.80     8.84
```

Reading it against the code:

**`idx` is the row index of `output/Hydro_ioniz.txt`, and row 1 is a ghost
cell, not the first interior cell.** The file is written over `j = 1-Ng .. N+Ng`
with `Ng = 2`, so rows 0 and 1 are the two lower ghosts and row 2 is `j = 1`.
Checked against the stored golden of the same case: rows 0 and 1 carry the
same `rho`, `v` and `p`, and rows 0, 1, 2 carry the *same velocity*
(120.50092137458105 cm/s, against 241.19 in row 3) -- the signature of
`W_in(2,index) = max(W_in(2,1), 0)`, the ghost velocity copied from the first
interior cell. Item (P)'s table starts at row 1, so its `idx 1` is the ghost
`j = 0` and its `idx 2` is the first interior cell.

**`rho` at the base is the pinned reservoir value and nothing else.**
2.5000e14 mH/cm^3 = `rho_bc` x `n0` with `rho_bc = (1 + 4 x 1)/(1 + 1) = 2.5`
and `n0 = 1e14 cm^-3`, exact. The density at the base does not evolve; it is a
Dirichlet condition.

**The velocity is not decaying upward -- the Mach number is rising because the
sound speed is collapsing.** From row 1 to row 234, `v` *falls* from 1.081e7
to 7.955e6 cm/s while the Mach number rises from 38.6 to 224.0. The whole rise
is in `c = sqrt(gamma p / (rho m_H))`, which falls from 2.80e5 to 3.55e4 cm/s
as the temperature drops from 950 K to 30.6 K (`c` is proportional to the
square root of `T` at fixed mean mass, and `T` falls by a factor 31, i.e. `c`
by 5.6, while `v` falls by only 1.36). The layer is not accelerating; it is
freezing.

**The last row is a shock front, not a boundary artifact.** Across rows 234 to
235 the density rises by 6.4x (and by 20x above the cells just below it), the
temperature jumps 30.6 -> 8532 K, the velocity drops, and the Mach number
falls from 224 to 8.8. Density up, temperature up, velocity down, Mach down
across one interface is a compression front. The cold hypersonic stream below
`r = 1.1425` is running into it.

**The particle count says what the gas is at each level.** From `p = (n_tot +
n_e) k_B T`, the particle count per (H+He) nucleus is `p/(rho k_B T)` x
`rho_bc`:

| row | particle count per nucleus | reading |
|---|---|---|
| 1 (ghost) | 1.500 | *not physical here*: needs `n_e` = 0.5 per nucleus, i.e. hydrogen fully ionized, at 950 K in a shielded base |
| 2 (cell 1) | 1.500 | the same, at 195 K |
| 3 | 0.750 | exactly `ntot_bc` at He/H = 1: hydrogen entirely in H2, `n_e` ~ 0 |
| 233 | 0.750 | the same |
| 234 | 0.750 | the same |
| 235 | 1.000 | atomic and neutral: every nucleus a free particle |

So the deep molecular layer and the shocked cell are internally consistent and
physically readable; the ghost and the first interior cell are not. Section 7
takes that up.

**The state is not a wind.** Item (P) records `du = 12.8` at the time of the
abort and section 111 records `du = 13.5` for the completed 12000-step run --
`du` is the fractional radial spread of `rho v r^2`, and a converged wind has
it at 1e-3. Independently: the recorded base row carries a mass flux
`4 pi R0^2 rho v` = 6.7e17 g/s, or log10 = 17.52 in the run's own
`Rate/2 + Mdot/2` convention, against the log10 Mdot = 8.89 the same run
reports. The base is draining the reservoir 8.6 decades faster than the wind
above it carries mass away.

---

## 2. What the lower boundary condition imposes

### 2.1 The three components

`BC_component_constrho` (`src/modules/states/Apply_BC.f90`, lines 41-101),
called for both lower ghosts by `Apply_BC_W` and again for the reconstructed
face states by `Rec_BC`:

* **Density.** `W_in(1,index) = rho_bc`, unconditionally. There is no branch.
  `rho_bc = comp_rho_bc() = comp_mass_per_H()/(1 + HeH)` (`composition.f90`),
  the gas mass per (H+He) nucleus in units of m_H. It is a constant of the
  run, computed once in `input_read`.
* **Velocity.** Three branches. The default -- and the one every molecular
  case uses, since none of them sets the keys -- is the legacy one-way valve
  `W_in(2,index) = max(W_in(2,1), 0.0)`: the ghost copies the first interior
  cell when the gas moves outward and is clamped at zero when it moves inward.
  `valve_eps > 0` replaces the kink by the smooth
  `0.5*(v_1 + sqrt(v_1^2 + eps^2))`, which is the same copy for `v_1 >> eps`.
  `base_v_massflux` (`Base velocity: massflux`) is the only branch that does
  *not* read `v_1`: it sets `v = F_c/(rho_bc r^2)` from the mass-flux constant
  measured over the constant-momentum region `[j_min:N]`.
* **Pressure.** Three branches. Default (legacy, `Base ghost temperature:
  isothermal`): `W_in(3,index) = ntot_bc + dp_bc`, a constant, which at the
  pinned `rho_bc` is the statement `T_base = T0`. `hydrostatic_base`
  extrapolates the interior gradient from cells 1 and 2 so that `dp/dr` is
  continuous through the base face (with an admissibility guard added in
  section 111). `base_ghost_T_continuous` sets
  `(ntot_bc + dp_bc) * p_1/n_part_cell1`, i.e. `dT/dr = 0`.

The regression case this item lives on, `backup/regression/mol_diffusion`,
sets none of the optional keys, so it runs the legacy triple: **`rho` pinned,
`p` pinned, `v` copied from the interior with a floor at zero.**

Verified against stored output: the golden of that case has
`p_ghost = 8.99394` erg/cm^3, and `ntot_bc x n0 k_B T0` with the run's own
`ntot_bc = 0.571428571` (read from `EXHALE_resolved.out`), `n0 = 1e14`,
`T0 = 1140 K` is 8.9939 -- six digits. The pressure at the base is the pinned
constant, not an evolved quantity.

### 2.2 The characteristic count, and where it stops being right

The lower boundary is an *inflow* boundary whenever `v > 0`: gas enters the
domain from below. The three characteristic speeds of the 1-D Euler system are
`v - c`, `v`, `v + c`, and a boundary condition is well posed when the number
of conditions imposed equals the number of characteristics entering the
domain.

* **`0 < v < c` (subsonic inflow).** `v` and `v + c` enter; `v - c` leaves.
  Two conditions must be imposed and one component must be determined by the
  interior. The code imposes `rho` and `p` and takes `v` from cell 1. **This
  is the correct count**, and it is why the base works: the outgoing `v - c`
  characteristic is the path by which the interior tells the boundary what it
  is doing, and the pinned pressure is the restoring force acting against it.
  The "breathing base" limit cycle described in `EXHALE_BC_and_IC.tex`
  section 2.4 is that restoring loop in action -- an oscillation about `v = 0`,
  bounded.
* **`v > c` (supersonic inflow).** All three characteristics enter. The
  boundary state must now be *fully* specified, and no component may be taken
  from the interior, because no information travels from the interior to the
  boundary at all. The code nevertheless keeps `v_ghost = v_1`. The copy is
  then not a closure but a feedback: the reservoir supplies mass at the fixed
  density `rho_bc` at whatever speed the interior happens to have, and the
  supply rate `rho_bc v_1 r^2` grows with `v_1` with nothing acting against
  it. The pinned pressure, which is the restoring term in the subsonic case,
  can no longer reach the interior.

**Nothing in the code tests the base Mach number.** A search of
`src/modules/states/`, `src/modules/time_step/` and `parameters.f90` finds
`Mach` only in the optional low-Mach contact dissipation (`lowmach_damp_*`,
a *low*-Mach device, default off), in a comment in `Reconstruction.f90`, and
in `viscous_conduction.f90`'s outer boundary. There is no diagnostic, no
warning and no guard at the base. That is the concrete form of item (P)'s
question, and section 9 recommends the counter.

This section states a well-posedness property of the boundary condition, which
is a property of the equations and holds independently of the composition. It
does **not** say that the boundary condition *causes* the crossing -- only that
once the crossing happens the boundary has no mechanism to undo it. Which
closure lets `v_1` get there in the first place is section 5.

### 2.3 What the boundary condition fixes and the composition cannot move

Because both `rho` and `p` at the base are constants of the run, so is the
base sound speed:

    c_bc = sqrt(gamma * p_bc / rho_bc_cgs)
         = sqrt(gamma * k_B * T0 * ntot_bc / (rho_bc * m_H))

with `gamma = g = 1.666666666667` and `m_H = mu = 1.67353284e-24 g`
(`parameters.f90` lines 431-432). Evaluated from the code's own base
composition functions for the hot-Uranus gate (`T0 = 1140 K`, `n0 = 1e14`,
`q_H2_base = 0.75`, no metals):

| He/H | `ntot_bc` | `rho_bc` | mean mass [m_H] | `p_bc` [erg/cm^3] | `c_bc` [cm/s] | `v` for M = 1 |
|---|---|---|---|---|---|---|
| 0.0793 | 0.57143 | 1.2204 | 2.136 | 8.994 | 2.709e5 | 2.709e5 |
| 1 | 0.75000 | 2.5000 | 3.333 | 11.805 | 2.169e5 | 2.169e5 |
| 10 | 0.95455 | 3.7273 | 3.905 | 15.024 | 2.004e5 | 2.004e5 |
| 100 | 0.99505 | 3.9703 | 3.990 | 15.662 | 1.982e5 | 1.982e5 |
| 1000 | 0.99950 | 3.9970 | 3.999 | 15.732 | 1.980e5 | 1.980e5 |

(The `p_bc` column reproduces the stored golden at He/H = 0.0793 to six
digits, which is the check that the formula is the code's.)

Two consequences.

1. **The base mean molecular weight cannot make a window.** It rises
   monotonically, 2.14 -> 4.00, saturating at pure helium, and the sound speed
   it sets falls monotonically by a factor 1.37 across four decades of
   composition. Item (P) named "the base pressure/density BC pair under a mean
   molecular weight that has moved by 5x" as a candidate; the measured move is
   1.87x in the mean mass and 1.37x in the threshold velocity, and both are
   monotonic. A monotonic factor of 1.37 does not produce a failure at 1 and
   10 and a clean run at 100.
2. **The sonic threshold is about 2e5 cm/s at every composition tested.** The
   recorded base velocity 1.08e7 cm/s exceeds it by a factor 50. Whatever
   drives the base, it does not need to work harder at He/H = 1 than at 100 to
   cross the line.

For contrast, the same quantity measured directly from the eight stored
regression goldens (recomputed for this document from their own `p` and
`rho`):

| case | base Mach | max Mach | `c_bc` [cm/s] | `v_ghost` [cm/s] |
|---|---|---|---|---|
| `lower_profile` | 0.0000 | 0.941 | 4.008e5 | 0.000 |
| `mol_base_handoff` | 0.0004 | 2.048 | 2.709e5 | 120.5 |
| `mol_diffusion` | 0.0004 | 2.048 | 2.709e5 | 120.5 |
| `mol_ir_bands` | 0.0005 | 2.112 | 2.696e5 | 134.4 |
| `mol_lyman_werner` | 0.0000 | 1.898 | 2.709e5 | 0.000 |
| `mol_metals` | 0.0005 | 2.113 | 2.696e5 | 134.9 |
| `wasp_full` | 0.0021 | 1.010 | 5.098e5 | 1065 |
| `wasp_he23off` | 0.0020 | 1.012 | 5.098e5 | 1022 |

A healthy base sits four orders of magnitude below its own sound speed. The
`mol_lyman_werner` and `lower_profile` bases are exactly at rest, which is the
valve holding a momentary inflow at zero.

---

## 3. The composition dependence

### 3.1 The atomic ladder, re-measured for this document

Item (P) quotes an atomic ladder as evidence that helium alone is innocent.
Those six runs are stored, so the table was recomputed here from their own
`output/Hydro_ioniz.txt` (`c = sqrt(gamma p / (rho m_H))`, `gamma = 5/3`, the
code's `g`; physical cells and the base only):

| arm | He/H | base Mach | max Mach | at | min T [K] | `v_base` [cm/s] | `c_bc` [cm/s] |
|---|---|---|---|---|---|---|---|
| `heh0p55` | 0.55 | 0.000 | 0.821 | 29.05 R_p | 133.0 | 0.169 | 1.227e5 |
| `heh2p13_diff_kzz1e9` | 2.13 | 0.000 | 0.813 | 29.05 R_p | 226.0 | 0.072 | 1.011e5 |
| `heh10` | 10.0 | 0.000 | 0.776 | 29.05 R_p | 184.8 | 0.097 | 9.13e4 |
| `heh100` | 100.0 | 0.000 | 0.793 | 29.05 R_p | 226.0 | 0.068 | 8.85e4 |
| `heh1000` | 1000.0 | 0.000 | 0.793 | 29.05 R_p | 226.0 | 0.063 | 8.82e4 |
| `heh10000` | 10000.0 | 0.000 | 0.801 | 29.05 R_p | 226.0 | 0.067 | 8.82e4 |

The values reproduce item (P)'s table exactly (base Mach 0.000, max 0.8 at
29.0 R_p, min T 133/226/185/226/226/226). The base velocities are of order
0.1 cm/s -- six orders below the sound speed -- and the flow stays subsonic to
the outer boundary. Raising He/H by four decades does not perturb the base.

### 3.2 The molecular matrix

From section 93 (the same gate, the same binary, only the composition
changed; 12000-step cap, `OMP_NUM_THREADS=1`). **Carried over; the run
directories were not kept.**

| He/H | NaN abort | log10 Mdot |
|---|---|---|
| 1 | yes, step 3406 | 8.91 |
| 10 | yes, steps 5474-8827 | 8.87 |
| 100 | no, 12000 steps | 11.26 |
| 1000 | no, 12000 steps | 11.09 |

Two groups separated by more than two decades in `Mdot`. Since section 111 the
He/H = 1 arm completes its 12000 steps (0 steps retaken, 0 dt bisections,
8156 interface fluxes dropped to first order, log10 Mdot 8.89), so the *abort*
is gone; the *state* is not, and `du` stays at 13.5.

### 3.3 What the non-monotonicity excludes, and what it does not

It excludes anything monotonic in He/H as the sole cause. Measured or derived
above, the following are all monotonic across the window and therefore cannot,
alone, produce a failure at 1 and 10 with a clean run at 100:

* base mean molecular weight (2.14 -> 4.00) and base sound speed (2.71e5 ->
  1.98e5 cm/s), section 2.3;
* the H2 fraction of the base particle count, `f_H2 = 0.750, 0.333, 0.048,
  0.005, 0.0005` at He/H = 0.0793, 1, 10, 100, 1000 (derived from
  `h2_bound_fraction`), which falls monotonically;
* the third-body correction of section 82 itself: it divides R12/R13/R15 by
  the base mean mass, 2.14 -> 4.00, again monotone.

A window therefore requires a **product of two opposing trends**: something
that needs hydrogen (so it dies above He/H ~ 10) multiplied by something that
needs helium (so it is absent below He/H ~ 1). The obvious candidate class is
"a molecular layer that is still thick enough to matter, in a gas whose
thermodynamics is already helium's". Item (P) says the same thing in words;
nothing measured here identifies which quantity it is.

Note also what *triggers* the abort, from section 93: the tree **as found**
(third body = `rho/m_H`) completed at every He/H tested; only with the
physically correct third body (`n_tot`) do He/H = 1 and 10 abort. The
correction is right and stays; it exposed the state rather than creating it.

### 3.4 The two ladders are not the same experiment

This is a limitation of the evidence, not a result. The atomic ladder of
section 3.1 is **LHS 1140 b** -- `T_eq = 226 K`, `M = 0.0176 M_J`,
`R = 0.158 R_J`, PLM only, no `base.inp` handoff, `2D approximate method:
Mdot`. The molecular matrix of section 3.2 is the **hot-Uranus gate** --
`T_eq = 1140 K`, `M = 0.0457 M_J`, `R = 0.49 R_J`, PLM+WENO3, a `base.inp`
handoff, `Rate/2 + Mdot/2`. Comparing them changes the planet, the base
temperature by a factor 5, the reconstruction schedule, the handoff and the
molecular switch all at once.

The controlled comparison that would isolate the molecular switch -- the
hot-Uranus gate at He/H = 1 with `Molecular chemistry: False` -- **has not
been run**. It is diagnostic D1 in section 9 and it is the cheapest thing on
the list.

---

## 4. Causes excluded, and by what measurement

| candidate | status | the measurement |
|---|---|---|
| the molecular chemistry solver | **identified initiating defect** | Row scaling changed the MINPACK path, but `ionization_equilibrium.f90` still retains a physical-bounds iterate when `info /= 1`, without checking its final reaction residual. If all iterates lie outside the bounds, it projects one onto the element budget and again does not recompute the residual. Section 11 measures fully ionized H at 28--244 K in the affected first cell. |
| the reconstruction producing a negative pressure | **excluded as the cause of the abort** | Section 95: the negative pressure was already in the *cell average* (`p = -3.64985e-2`, the identical number in the cell and both of its faces, the slope limited to zero), so the reconstruction only passed it on. The real violation was update positivity at CFL 0.6 in a Mach 38-224 layer. `positivity_limited_faces` was added anyway and is not idle -- 1388 non-positive face states in 9000 steps of `wasp_full`, which HLLC had been swallowing -- but it never fires on this arm. |
| update positivity | **excluded as the cause of the supersonic state** | Section 111: with the local flux correction the same arm completes 12000 steps with 0 retaken steps, 0 dt bisections, no NaN, and 8156 corrected interface fluxes. It marches. `du` is still 13.5. Positivity was keeping the run alive, not making it supersonic. |
| helium-rich composition as such | **excluded** | Section 3.1, re-measured: six atomic arms from He/H = 0.55 to 10000, base Mach 0.000 in all six, maximum Mach 0.78-0.82. |
| the base mean molecular weight alone | **excluded** | Section 2.3, derived from the code's base composition functions: monotone, 1.37x in the sonic threshold over four decades. |
| the base boundary condition | **remaining** | Sections 2 and 5. |

---

## 5. Hypotheses for the window, and how each is decided

All four are hypotheses. None is supported by a measurement in this document
beyond what is stated, and they are not mutually exclusive.

### H1. The free velocity has no restoring path above M = 1

*Statement.* The subsonic base is stable because the outgoing `v - c`
characteristic carries the interior's back-pressure to a boundary whose
pressure is pinned; that is the breathing limit cycle. Some perturbation
(molecular cooling of the first cells, the front at 1.14 R_p, the cold start)
pushes `v_1` past 2e5 cm/s once, and from there the copied ghost velocity is a
feedback with no term acting against it: the reservoir feeds harder the faster
the interior runs.

*How to decide it.* `Base velocity: massflux` replaces the copy by
`v = F_c/(rho_bc r^2)`, which is set by the wind's own flux constant and does
not read `v_1` at all. If the same arm run with that key stays subsonic at the
base, the free-velocity feedback is the mechanism; if it still runs
supersonic, H1 is wrong and the driver is above the boundary. This is a
one-run test on an existing case (D3).

*Secondary signature.* Instrumenting `v_1` and the base Mach every step would
show whether the crossing is a single event early in the run followed by
monotone growth (H1) or a slow drift (H1 is then not the whole story).

### H2. The pinned base state and the interior chemistry disagree about what the base gas is

*Statement.* The pressure BC pins `p = ntot_bc + dp_bc` at `rho = rho_bc`,
which is the statement `T_ghost = T0` **for a gas whose particle count is
`ntot_bc`**. `ntot_bc` is computed once, in `input_read`, from
`q_H2_base = 0.75`; the interior cells get their particle count from the
molecular solver every step. Where the two disagree, the boundary carries a
contact discontinuity: the same `rho`, a different temperature. That is
exactly the failure `base_ghost_T_continuous` was introduced for -- on
HD 189733 b the pin sat at 1183 K against a first cell at 490 K, and removing
the jump cut the alternating base density amplitude by 18x
(`parameters.f90` lines 553-573, `docs/hd189_base_checkerboard.md` section
4.2 E2). In item (P)'s snapshot the first interior cell is at 195 K against a
`T0` of 1140 K, a factor 5.8 -- larger than the HD 189733 b case that
motivated the key.

*Why this could produce a window.* `ntot_bc` at He/H = 1 is 0.750 (hydrogen
entirely in H2) and at He/H = 100 it is 0.995 (essentially atomic, because
there is almost no hydrogen left to be molecular). The mismatch between "what
the pin assumes" and "what the chemistry returns" therefore has room to be
large only where there is still a lot of hydrogen *and* the pin has been
driven fully molecular -- which is the window.

*How to decide it.* Two runs. (a) He/H = 1 with `Base ghost temperature:
continuous`, which removes the `T0` pin while keeping the density anchor;
(b) He/H = 1 with `Hydrostatic base: True`, which replaces the pin by the
extrapolated interior gradient. If either keeps the base subsonic, the pinned
pressure closure is implicated. Read `EXHALE_resolved.out` for `ntot_bc` and
the first-cell particle count in each.

### H2b. `q_H2_base = 0.75` is out of range at He/H >= 0.167 and the code caps it silently

*This one is a fact about the experiment, not a hypothesis, and it weakens the
composition scan as a controlled comparison.* `h2_bound_fraction`
(`composition.f90` lines 284-296) converts the mixture mixing ratio `q` into
the fraction of H nuclei bound into H2,

    x2 = 2 q (1 + HeH)/(1 + q),   capped at x2 = 1,

and the cap is reached at `HeH = (1 - q)/(2 q)`, which for `q = 0.75` is
**0.1667**. The requested `x2` is 0.925 at He/H = 0.0793 (no cap), 1.71 at
He/H = 1, 9.43 at 10, 86.6 at 100 and 858 at 1000. The cap is right -- `q` and
`HeH` are not independent, and at He/H = 1 the largest attainable `q` is
0.5/(0.5 + HeH) = 0.333 -- but it is applied **without a message**: the
resolved-setup report prints the resulting `ntot_bc`, and the startup line
prints `q_H2(base, photochemical) = 0.750`, which at that composition is not
a value the base can have.

The consequence for the scan: raising `He/H` from 0.0793 to 1 changed *two*
things, the composition and the effective base molecular content (from 92.5 %
of H nuclei in H2 to 100 %, saturated by the cap). Every arm at He/H >= 0.167
in this matrix ran with a saturated base.

*How to decide its relevance.* Re-run He/H = 1 with `q_H2_base` removed from
`base.inp`, which selects the Visscher/Koskinen chemical-equilibrium fit at
`(p_base, T0)` instead, or with a value inside the attainable range. If the
base stays subsonic, the out-of-range handoff was doing the work.

*Recommendation (not applied -- `src/` untouched by this document):*
`h2_bound_fraction` should say when it caps, and the setup report should carry
the effective `x2` alongside `q_H2_base`. A silent clamp on an out-of-range
input is how a scan changes two variables while its author believes it changed
one.

### H3. Molecular cooling collapses the first cells and the sound speed with them

*Statement.* The Mach rise in section 1 is a sound-speed collapse:
`T` = 950 -> 195 -> 30.6 K over the layer. H3+ cooling (`h3p_cooling.f90`), the
H2 rovibrational bands and, when `Molecular IR bands` is on, H2O/CO
(section 110) all act at the base and all scale with the molecular density.
If at He/H ~ 1-10 the base is cold enough for a molecular layer but hot enough
for the layer to still radiate strongly, the first cells cool, `c` collapses,
and a base velocity that was marginally subsonic becomes strongly supersonic
without any change in `v`.

*Why a window.* Cooling of this kind needs H2 (so it dies above He/H ~ 10 with
the hydrogen) and needs the helium-loaded mean mass to make `c` small (so it
is weaker below He/H ~ 1). That is the required product of opposing trends
from section 3.3.

*How to decide it.* `output/Cooling_breakdown.txt` (written by every run) at
the first ten cells, for the arms of the D-list below: read the H3+ and H2
channels as a fraction of the total, and the resulting `T(1)/T0`. If the base
temperature ratio is non-monotonic in He/H with a minimum in 1-10, H3 is
supported. This costs nothing beyond the runs already needed for D1/D2.

### H4. The shock at 1.14 R_p is the object, and the supersonic layer is its inflow

*Statement.* The alternative causal order: a dense hot wall forms somewhere
above the base (row 235: 20x density, 8532 K), the column below it drains
onto it, and the base velocity is the *consequence* of that drain rather than
its cause. The boundary condition would then be an accomplice (it happily
supplies the mass) rather than the driver.

*How to decide it.* Time resolution, which the stored snapshot cannot give:
write `Hydro_ioniz.txt` every 100 steps for the first 3500 steps of the
He/H = 1 arm and find which forms first -- the base crossing M = 1, or the
front. If the front is first, H4; if the base is first, H1/H2/H3.

---

## 6. Scope: what is exposed, measured

**Nothing in production is inside the window.** But the reason is not the
composition:

* **LHS 1140 b.** 942 `input.inp` files under `LHS1140b/`. **646 of them
  (69 %) are at He/H between 1 and 10** -- the values 2.09 (61 files), 2.13
  (40), 9.7 (34), 11.1 (35, just outside), 3.0 (22), 8.0 (20), 9.0 (18),
  8.45 (17) dominate the ladder. **Zero of the 942 have `Molecular chemistry:
  True`.** This confirms, by direct count, what
  `docs/molecular_chemistry_audit_he_rich.md` section 4b states.
* **Every molecular configuration in the repository.** Six `input.inp` files
  carry `Molecular chemistry: True`: `examples/15_molecular`,
  `examples/16_molecular_metals`, `examples/18_oxygen_chemistry`,
  `examples/19_molecular_ir_bands`, `docs/lower_atmosphere_figs/data_g1m`,
  `docs/lower_atmosphere_figs/data_g2`. All are at He/H <= 0.083333.
* **The regression matrix.** All five molecular cases (`mol_base_handoff`,
  `mol_metals`, `mol_lyman_werner`, `mol_diffusion`, `mol_ir_bands`) are at
  He/H = 0.0793, base Mach 0.0000-0.0005 (section 2.3).
* **The A2 oxygen work** (sections 108-111) is HD 189733 b at He/H =
  0.083333333.

**What the item blocks.** The helium-rich molecular regime itself: a
sub-Neptune with He/H of order unity and a molecular lower atmosphere cannot
be run today, because that combination is exactly the window. The LHS 1140 b
ladder is the obvious place this would be wanted -- 646 of its runs are
already at those compositions and would need only the switch. Until this item
is closed, molecular chemistry and He/H in 1-10 must not be combined without
checking the base Mach number of the result.

---

## 7. A defect in the recorded snapshot

Recorded here because the arm has to be re-run anyway and this should be
checked when it is.

**The base row's pressure is not a value the base pressure BC can produce.**
The `mol_diffusion` configuration uses the legacy pin, so
`p_ghost = (ntot_bc + dp_bc) * p0` with `dp_bc = 1e-10` (`input_read.f90`
line 1277). At He/H = 1 with `q_H2_base = 0.75`, `ntot_bc = 0.750` and
`p0 = n0 k_B T0 = 15.7394`, so `p_ghost = 11.805` erg/cm^3. Item (P)'s table
records **19.674**, which is `1.2501 * p0`. `comp_ntot_bc` cannot return
1.2501: it starts at 1 (or `(1 + HeH + sum(ab))/(1 + HeH)` with metals, which
is 1.0005 here) and `h2_bound_fraction` only ever subtracts. The same table's
density, 2.5000e14, *is* exactly `rho_bc * n0`, so the composition and the
normalization are the expected ones.

**The particle count of that row and the next is not physical.** Section 1's
table: rows 1 and 2 carry 1.500 particles per nucleus, which at He/H = 1 needs
`n_e` = 0.5 per nucleus -- every hydrogen nucleus ionized -- at 950 K and
195 K in a base the run shields. Rows 3, 233 and 234 carry exactly 0.750
(= `ntot_bc`, hydrogen entirely molecular, `n_e` ~ 0) and row 235 exactly
1.000 (atomic, neutral). The two anomalous rows are precisely the ghost and
the first interior cell.

Three readings, none of which can be tested without the run:

1. **The dump is not a state after `Apply_BC`.** The run aborted on its NaN
   detector mid-step; if the arrays were written from a partial update the
   ghost need not carry its boundary value. Under this reading the recorded
   Mach numbers are correct as printed and item (P)'s table stands.
2. **The printed pressure carries a factor `gamma`.** 19.674/11.805 = 1.66665,
   which is `g` to five digits. Under this reading the true base pressure is
   the pinned 11.805, the particle counts of rows 1 and 2 become 0.900
   (10 % of nuclei bound in H2, `n_e` ~ 0, physical), and every recorded Mach
   number is low by `sqrt(gamma)` = 1.291: base 49.8 instead of 38.6, and
   289 instead of 224 at row 234. Against this reading, row 235's own
   particle count is exactly 1.000 with the pressure *as printed*, and would
   become 0.600 -- below the fully molecular limit 0.750 and therefore
   impossible -- if the factor were global.
3. **The arm did not run the legacy pin.** `hydrostatic_base` extrapolates and
   is not bounded by `ntot_bc`; but extrapolating item (P)'s own rows 2 and 3
   to `r = 1.0` gives 0.264, not 19.674, so this reading is not supported by
   the table itself.

Reading 1 is the most economical. Either way, **the conclusion of item (P) is
unaffected**: the base sound speed that the boundary condition pins at this
composition is 2.169e5 cm/s, and the recorded base velocity of 1.081e7 cm/s
is 50 times it.

*What to do about it:* the re-run of section 9 should keep its output
directory, and should be stopped cleanly (a step cap) rather than on an abort,
so the dumped ghost is a post-`Apply_BC` state.

---

## 8. What is settled, and what is not

**Settled.**

* The base boundary condition imposes `rho` and `p` and copies `v` from the
  first interior cell. Read from `Apply_BC.f90` and confirmed against stored
  output (three consecutive rows sharing one velocity; the pinned pressure
  reproducing `ntot_bc * p0` to six digits).
* That is the well-posed condition count for subsonic inflow and not for
  supersonic inflow, and no code path tests which regime the base is in.
* The base sound speed is a constant of the run, 2.71e5 -> 1.98e5 cm/s over
  He/H = 0.0793 -> 1000 for this planet; the sonic threshold is ~2e5 cm/s
  throughout.
* The chemistry solver, the reconstruction and update positivity are all
  excluded as causes of the supersonic state (section 4).
* Helium alone is excluded: six atomic arms from He/H = 0.55 to 10000 have a
  base Mach of 0.000 (re-measured here).
* The base mean molecular weight alone is excluded: monotone, 1.37x in the
  threshold velocity over four decades.
* Scope: 0 of 942 LHS 1140 b runs and 0 of the 6 molecular configurations in
  the repository are inside the window; 646 of the 942 are at the window's
  *compositions* and are outside it only because molecular chemistry is off.

**Not settled.**

* Which closure of the boundary condition -- the copied velocity, the pinned
  pressure, or neither -- first lets `v_1` reach 2e5 cm/s.
* Whether the base crosses M = 1 before or after the front at 1.14 R_p forms.
* Whether the window has anything to do with the silently capped
  `q_H2_base = 0.75` (H2b), which made every arm at He/H >= 0.167 run with a
  saturated molecular base.
* The window's edges. Only He/H = 1, 10, 100, 1000 were run; the boundary is
  somewhere in (10, 100) on one side and below 1 on the other, unmeasured.
* Whether the same thing happens on a different planet with molecules on and
  He/H ~ 1 -- the two ladders compared in item (P) differ by more than the
  molecular switch (section 3.4).
* The internal inconsistency of the recorded base row (section 7).

---

## 9. Next steps, cheapest first

Each of D1-D4 is one 12000-step marching run of an existing case with one key
changed, single-threaded, and each is decisive for one hypothesis. Keep the
output directory of every one; write `Cooling_breakdown.txt` and
`EXHALE_resolved.out` with it.

**D1 -- the missing control (decides section 3.4).** The hot-Uranus gate at
He/H = 1 with `Molecular chemistry: False`. This is the comparison item (P)
does not have: same planet, same base level, same grid, molecules the only
difference. If it stays subsonic, the molecular layer is necessary and the
LHS 1140 b ladder is irrelevant to the diagnosis; if it goes supersonic, the
problem is helium plus this planet and the molecular network is a bystander.
Do this first -- it halves the hypothesis space either way.

**D2 -- remove the out-of-range handoff (decides H2b).** He/H = 1 with
`q_H2_base` deleted from `base.inp` (falling back on the equilibrium fit) or
set to a value below 0.5/(0.5 + HeH) = 0.333. Report the resulting `ntot_bc`
from `EXHALE_resolved.out` next to the base Mach.

**D3 -- remove the velocity feedback (decides H1).** He/H = 1 with
`Base velocity: massflux`. The ghost velocity then comes from the wind's flux
constant and never reads `v_1`. Note the recorded cost: this key slows
convergence about 5x (`parameters.f90` line 640 and its note), so it is a
diagnostic, not a proposed default.

**D4 -- remove the temperature jump (decides H2).** He/H = 1 with
`Base ghost temperature: continuous`, and separately with
`Hydrostatic base: True`. The first cell in item (P)'s snapshot is at 195 K
against `T0 = 1140 K`, a larger jump than the HD 189733 b case that motivated
the key.

**D5 -- the causal order (decides H4).** Whichever of D1-D4 still goes
supersonic: dump the profile every 100 steps to step 3500 and record when the
base first crosses M = 1 and when the front at ~1.14 R_p first appears.

**D6 -- the window edges.** He/H = 0.3, 3, 30 added to the existing 1, 10,
100, 1000. Cheap and it constrains section 3.3's "product of two opposing
trends" more than any single run does.

**Two code changes to consider afterwards** (neither made here; `src/` was not
touched):

* a base Mach diagnostic in the same style as the positivity counters --
  reported at the end of a run, and a warning the first time the base face
  goes supersonic. A run that silently marches a Mach 224 base for 12000 steps
  should not be able to do so quietly. This is a diagnostic, not a guard: what
  the boundary condition should *do* in that regime is the open physics
  question, and a counter does not prejudge it.
* a message from `h2_bound_fraction` when `x2` is capped, and the effective
  `x2` in the setup report (H2b).

---

## 10. Files and reproduction

**Code read for this document**

| File | What it settles |
|---|---|
| `src/modules/states/Apply_BC.f90` | `BC_component_constrho`: `rho = rho_bc`, `v = max(v_1,0)` (or valve / mass flux), `p = ntot_bc + dp_bc` (or hydrostatic / continuous-T). `Apply_BC_W`, `Rec_BC` |
| `src/modules/functions/composition.f90` | `comp_rho_bc`, `comp_ntot_bc`, `h2_mixing_ratio_base`, `h2_bound_fraction` (and its cap) |
| `src/modules/files_IO/input_read.f90` | lines 1224-1282: `mass_per_H`, `ntot_bc`, `rho_bc`, `dp_bc = 1e-10`, `p0 = n0 mu v0^2`, the `thereis_mol -> molecular_base` coupling |
| `src/modules/init/parameters.f90` | `g = 1.666666666667`, `mu = 1.67353284e-24`, `RJ = 6.9911e9`; the base-BC option block (lines 534-645) |
| `src/modules/flux/Num_Fluxes.f90` | `a = sqrt(g p/rho)`, the sound speed the code itself uses |
| `src/modules/files_IO/write_output.f90` | the `Hydro_ioniz.txt` schema: column 2 is a mass density in units of m_H per cm^3 |

**Stored output measured**

| Path | Used for |
|---|---|
| `LHS1140b/exhale/{heh0p55, heh2p13_diff_kzz1e9, heh10, heh100, heh1000, heh10000}/output/Hydro_ioniz.txt` | the atomic ladder, section 3.1 |
| `backup/regression/golden/*/Hydro_ioniz.txt` | the eight base Mach numbers, section 2.3 |
| `backup/regression/mol_diffusion/{input.inp, base.inp, output/, EXHALE_resolved.out}` | the case definition, the ghost layout, `ntot_bc = 0.571428571` |
| `LHS1140b/**/input.inp` (942 files) | the scope count, section 6 |

**Reproducing the Mach tables**

```python
import numpy as np
m_H = 1.67353284e-24          # parameters.f90 'mu'
g   = 1.666666666667          # parameters.f90 'g'
d = np.loadtxt('output/Hydro_ioniz.txt')      # '#' lines are comments
r, rho, v, p, T = d[:,0], d[:,1], d[:,2], d[:,3], d[:,4]
cs = np.sqrt(g*p/(rho*m_H))   # rho is in m_H per cm^3
M  = abs(v)/cs
# rows 0,1 are the lower ghosts (Ng = 2); row 1 is the base face state
```

and the base scalars, which need no run at all:

```python
q, HeH, T0, n0 = 0.75, 1.0, 1140.0, 1e14
x2   = min(2*q*(1+HeH)/(1+q), 1.0)            # composition.f90 h2_bound_fraction
ntot = 1.0 - 0.5*x2/(1+HeH)                   # comp_ntot_bc
rhob = (1 + 4*HeH)/(1 + HeH)                  # comp_rho_bc
c_bc = np.sqrt(g*1.380649e-16*T0*ntot/(rhob*m_H))
p_bc = ntot*n0*1.380649e-16*T0
```

**Related records**

* `TO_BE_DONE.md` items (O), (P), (Q) -- (O) and (Q) resolved, (P) is this
  document.
* `docs/Update_EXHALE.md` sections 82, 93, 95, 110, 111.
* `docs/molecular_chemistry_audit_he_rich.md` sections 4.1 (what M is), 4b
  (no LHS 1140 b run has molecules on), 7.
* `docs/EXHALE_BC_and_IC.tex` section 2 (the lower boundary), 2.4 (the
  breathing base), 2.6 (the optional stabilizers and the selectable base
  mode).
* `docs/hd189_base_checkerboard.md` section 4.2 E2 (the `T0` pin as a contact
  discontinuity, the precedent for H2).

---

## 11. Follow-up code audit and controlled runs

### 11.1 Revised physical judgment

The molecular equilibrium solve is the initiating defect. The affected cell
is allowed to use a bounded but nonconverged MINPACK iterate as though it were
chemical equilibrium. In the measured cases this makes hydrogen fully ionized
in a cold, dense, shielded cell. Recombination and free-free cooling calculated
from that state then remove energy at an enormous rate, the pressure and sound
speed collapse, and the hydrodynamic boundary condition either copies the
resulting velocity into the reservoir or hides it behind a different ghost
closure.

The lower boundary analysis in sections 2 and 5 remains correct: copying an
interior velocity into a supersonic inflow is not a well-posed boundary
condition. It is, however, the amplifier rather than the first cause in this
case. Changing only the hydrodynamic boundary cannot turn a non-root chemical
state into a physical atmosphere.

### 11.2 Verification setup

The repository `EXHALE.x` and
`backup/regression/mol_diffusion/EXHALE.x` had the same ELF BuildID
`14ae663c1ec59691762bc0466db38a49bff370c5`; the repository executable was
built on 2026-08-31 from the current source tree. Every controlled case was a
copy of `backup/regression/mol_diffusion` under `/tmp`, with
`OMP_NUM_THREADS=1` and `EXHALE_MAXSTEPS=12000`. No stored regression output
or production result was overwritten.

All rows below use the code's own sound speed,

    c = sqrt(gamma p/(rho m_H)),  gamma = 5/3,

from the final `Hydro_ioniz.txt`. The `r=1` row is the second lower ghost and
the first physical cell is at `r=1.000193 R_p`.

| case | isolated change at He/H = 1 | base-ghost Mach | first-cell Mach | first-cell T [K] | maximum Mach | minimum T [K] | result |
|---|---|---:|---:|---:|---:|---:|---|
| D1 | molecular network off; passive `Molecular base: True` retained | 0 | 9.15e-4 | 853.2 | 1.80 at the outer wind | 819.4 | the molecular network is necessary for the base failure |
| D2 | attainable `q_H2_base=0.30` instead of 0.75 | 73.1 | 153.7 | 132.1 | 530.6 | 11.4 | the silent H2 cap is not the cause |
| D3 | `Base velocity: massflux` | 3.34e-6 | 0.0115 | 11.4 | 1.81 at the outer wind | 11.4 | suppresses the copied-velocity feedback, but leaves a nonphysical cold base |
| D4a | `Base ghost temperature: continuous` | 35.5 | 38.9 | 11242 | 164.0 | 388.7 | does not cure the failure and overheats the base |
| D4b | `Hydrostatic base: True` | 0 | 1.58 | 27.95 | 2.19 | 0.044 in the ghost | removes outward ghost inflow but leaves a supersonic, nonphysical first cell |

D1 is a controlled molecular-network switch, not an atomic-EOS switch:
`Molecular base: True` keeps `ntot_bc=0.75`, the same base particle count as
the He/H=1 molecular case. Its stable base therefore cannot be attributed to
changing the lower pressure anchor back to an atomic value.

D2 also corrects an ambiguity in the original diagnostic plan. Removing
`q_H2_base` would select the equilibrium fit, whose value can still exceed
the attainable mixture limit at He/H=1 and be capped. Setting `q_H2_base=0.30`
is below the maximum 1/3 and gives `ntot_bc=0.769230769`. The failure becomes
more extreme, so the out-of-range value in the original scan is a real input
validation problem but not this runaway's cause.

### 11.3 Direct evidence in the retained state

The first physical cell in `Ion_species.txt` has the following hydrogen-nucleus
fractions:

| case | H I/H | H II/H | H in H2/H | first-cell T [K] |
|---|---:|---:|---:|---:|
| D2, valid H2 input | 0 | 1.000000 | 8.3e-55 | 132.1 |
| D3, mass-flux velocity | 0 | 1.000000 | 0 | 11.4 |
| D4b, hydrostatic pressure | 0 | 1.000000 | 7.6e-24 | 27.95 |

These are inside the element simplex, because no fraction is negative and the
hydrogen fractions sum to one. That only proves conservation. It does not
prove that the reaction balances vanish. Complete hydrogen ionization at
11--132 K in the dense base is not a physically acceptable equilibrium for
this shielded configuration.

The cooling ledger shows the consequence rather than an independent coolant
problem. In D4b at 27.95 K the first-cell total cooling is
`3.8494e-3 erg cm^-3 s^-1`: recombination contributes `3.3502e-3`, free-free
emission contributes `4.9918e-4`, and H3+ contributes only `2.19e-47`. In D3
at the 11.4 K temperature floor, recombination and free-free cooling are
`3.4988e5` and `4.6354e4 erg cm^-3 s^-1`, respectively. These rates are large
only because the accepted composition supplies an unphysical electron and
proton reservoir.

An additional scratch build set `h3p_cooling_rate=0` and was sampled through
step 5525. The saved 5001-step profile still had base-ghost Mach 46.4,
first-cell Mach 91.6, and a minimum temperature of 96.7 K. This directly
excludes H3+ emission as the initiating cause. The experiment used a modified
temporary executable and is diagnostic only; no source file in the repository
was changed for it.

### 11.4 The code path that admits the state

The relevant selection block is in
`src/modules/radiation/ionization_equilibrium.f90`:

1. Three molecular starts are tried: the prior state, a molecular balance
   estimate, and an atomic balance estimate.
2. `conv_ieq` is true only for `info == 1`.
3. Any iterate inside the element bounds receives `ok_rank=1`, even when
   `info /= 1`. A converged bounded iterate receives rank 2.
4. If no rank-2 result exists, the best rank-1 iterate is copied into
   `sys_x`. Its final reaction residual is not tested.
5. If every iterate is outside the bounds, one is projected by
   `clamp_fractions_to_element_budget`. The reaction residual is not
   recomputed after that nonlinear projection either.

The comments call both retained objects roots, but the implementation does
not establish that claim. Bounds and conservation are necessary conditions;
they are not reaction equilibrium. Row scaling cannot correct this acceptance
rule. It can reduce the number of `info=4` exits while leaving one repeatedly
failing base cell, and one cell is enough to launch the pressure wave.

The thermodynamic coupling completes the causal chain:

1. `ioniz_eq` installs the accepted molecular/ionization fractions.
2. `get_species_densities` obtains `n_tot` and `n_e` from them.
3. `comp_p_from_T` rebuilds pressure as `(n_tot+n_e) k_B T`.
4. A jump from molecular H to H II roughly doubles the particle count at
   fixed H/He nuclei and temperature.
5. The energy solve then evaluates recombination and free-free cooling from
   the same false ion density, collapses the temperature, and reaches the
   11.4 K numerical floor in several variants.
6. The legacy velocity boundary copies the resulting interior velocity. Once
   it is supersonic, the characteristic mismatch described in section 2
   amplifies the state.

`q_H2_base` currently affects only `comp_ntot_bc`; it does not impose the
molecular partition in the inflowing ghost. Thus the hydrodynamic pressure
anchor may describe a molecular reservoir while `ioniz_eq` independently
turns that same ghost into ionized atomic gas. The D2 ghost illustrates the
contradiction: the resolved base particle count is molecular, but the retained
species state contains H II rather than H2.

### 11.5 Required changes

#### A. Refuse non-roots immediately

This is the smallest safe correction and should precede any attempt to make
the run converge:

* accept a molecular state only when it is inside the element bounds **and**
  its normalized reaction residual meets a stated tolerance;
* do not treat `info /= 1` as sufficient evidence of failure or success by
  itself, but never accept the iterate without an explicit residual check;
* after any projection, recompute every reaction row and reject the state if
  the residual changed beyond tolerance;
* if all three starts fail, stop with cell index, radius, temperature, density,
  `info`, maximum element violation, maximum normalized reaction residual, and
  the candidate fractions. Continuing with a known non-root produces a
  physically meaningless wind.

This change will initially make the He/H=1 case fail early and clearly. That
is the correct behavior until the constrained solver below exists.

#### B. Solve the constrained chemical problem, rather than projecting it

The durable solution is a nonlinear solve whose iterates cannot leave the
physical composition domain. A clean formulation uses nonnegative species
number densities as the unknowns, the H and He conservation equations as
explicit residual rows, charge neutrality for the electron density, and the
independent reaction-balance rows. Positive transformed variables or a
bounded trust-region method can enforce nonnegativity. A continuation from the
dense molecular limit toward the local radiation field is preferable to three
unconnected guesses when multiple mathematical basins exist.

The acceptance test must be applied to the original, untransformed reaction
and conservation equations. A small transformed step or a successful library
status is not by itself evidence of chemical equilibrium.

The existing fraction layout can be retained temporarily only if the solver
enforces the coupled H and He simplex throughout the step, including the H and
He simultaneously carried by HeH+. Post-solve rescaling is not equivalent.

#### C. Give an inflowing molecular reservoir a composition boundary

When the lower boundary supplies mass, its elemental composition and molecular
partition are incoming data. If `base.inp` supplies `q_H2_base`, the ghost
species should be constructed from that value and `HeH_base`, with the
remaining stages obtained from a lower-atmosphere state or an explicitly
documented local equilibrium closure. The hydrodynamic EOS and the ghost
species must use the same particle count. A lower-atmosphere profile should
take precedence when present.

This changes `q_H2_base` from an EOS-only scalar to a physical composition
boundary and therefore changes input semantics. It should be implemented as a
deliberate interface revision, not as another hidden branch. If the intended
quantity remains EOS-only, the code must instead refuse a mismatch between
that EOS state and the species state used at the same ghost.

**IMPLEMENTED 2026-09-01 (section 117 of `docs/Update_EXHALE.md`).** C and
E were carried out together as one base-boundary block. Measured: the
equation of state had believed 92.51% of the base H nuclei were bound into
H2 while the chemistry delivered 99.97%, a 6.4% disagreement in the base
particle count which showed up as a base ghost at 1213.35 K where
T0 = 1140 K had been asked for. After C the ghost carries the stated
partition (x2 = 0.925114) and sits at 1140.0000 K exactly. E refuses an
unattainable request instead of capping it, and refuses nothing that was
running: 8 of 8 regression cases byte-identical, 27 of 27 tracked
`base.inp` files inside the ceiling. The scope recorded below was measured
before the work and held. Measured scope,
recorded here so the implementation starts from numbers rather than from a
survey:

* the ceiling is `q_H2,max = 0.5/(0.5 + He/H)`: **0.8631** at He/H = 0.0793,
  **0.3333** at He/H = 1, **0.0476** at He/H = 10;
* the five molecular regression cases (`mol_base_handoff`, `mol_metals`,
  `mol_lyman_werner`, `mol_diffusion`, `mol_ir_bands`) all carry
  `q_H2_base 0.75` at He/H = 0.0793, below the 0.8631 ceiling, so **E's
  refusal cannot break their goldens**. C's composition promotion may still
  move them, and that is to be measured when it is implemented;
* the tracked example `LHS1140b/examples/scalar_base_cno` sits at
  `q_H2_base 0.1927` with `HeH_base 2.0921`, whose ceiling is 0.1929 -- it
  passes at **99.9% of the limit**. The comparison at the boundary
  (inequality sense and round-off) must be written so that this example is
  not refused;
* what E *does* refuse is our own diagnostic arm: the He/H = 1 arm carries
  `q_H2_base 0.75` against a ceiling of 0.3333, **2.25x over**, and the code
  has been silently capping it.

**A question this raises for the window re-measurement, recorded as a
question and not as a claim.** If the He/H = 1 arm's base composition has all
along been a silently capped, physically impossible value, then that arm was
not a clean probe of the window. D2 above shows the cap is not the cause of
the runaway, and that stands. What is *not* established is how much the
phenomenon itself -- "a molecular base runs supersonic near He/H ~ 1-10" --
is entangled with an over-specified base H2 fraction. Answering it requires
re-specifying the arm with a physically admissible `q_H2_base` and
re-measuring, and that is the first item of the window re-measurement block.

#### D. Keep the hydrodynamic guard, but do not call it the chemistry fix

Add the base-face Mach diagnostic proposed in section 9 and warn on the first
supersonic inflow. `Base velocity: massflux` is useful as a diagnostic and
prevents the legacy copied-velocity feedback, but D3 proves that it can coexist
with a first cell at the temperature floor and enormous false recombination
cooling. A Mach guard must therefore be an independent validity check.

#### E. Validate the base H2 input

For a mixture value `q_H2=n_H2/(n_H2+n_H+n_He)`, the maximum at fixed He/H is

    q_H2,max = 0.5/(0.5 + He/H).

Refuse `q_H2_base > q_H2,max` instead of silently capping it, and write both
the requested `q_H2` and the implied fraction of H nuclei in H2 to
`EXHALE_resolved.out`. The original `q_H2_base=0.75` is impossible at He/H=1,
although D2 shows that this is not the runaway's cause.

### 11.6 Separate defect in `Base IR field`

Turning on `Base IR field` is physically motivated, but the current H3+ net
exchange implementation cannot be used as a workaround. The test reached NaN
at step 3. `h3p_emission_lte` holds the emission fit at its 30 K edge, while
`h3p_net_cooling_rate` continues to evaluate the absorption multiplier
`exp(E_k/T)` below 30 K, down to the 11.4 K trial floor. The intended low-
temperature detailed-balance cancellation is then lost and the heating term
diverges. In the measured evaluation at 11.4 K the dimensionless bracket is
about `-3.38e136`.

This should be repaired with a spectrally consistent absorption model or a
bounded net-exchange closure whose low-temperature and `T=T_rad, W=1` limits
are both finite and exact. Clamping only the exponential is not physically
sufficient because the total Miller emission fit contains many transitions
and cannot be paired with one 3--4 micron Boltzmann factor at all temperatures.

### 11.7 Focused acceptance tests after the solver change

The following checks follow the affected path and are sufficient before a
broader regression run:

1. Re-run the He/H=1 molecular hot-Uranus case and require zero accepted
   states whose normalized reaction residual exceeds tolerance.
2. At the base ghost and first ten physical cells, record the H/He species,
   particle count, temperature, cooling channels, and Mach number. No cell may
   be accepted merely because it lies inside the element bounds.
3. Require base inflow Mach below one throughout relaxation. Report the first
   crossing as a failed physical validation, not just a numerical message.
4. Require that recombination and free-free cooling vanish with the ion
   reservoir in the cold molecular limit; the 28 K, fully ionized state above
   must not recur.
5. Repeat D2 with `q_H2_base=0.30` to show that the result is independent of
   the old silent cap.
6. Repeat D3 only after the chemistry passes, to determine whether the
   mass-flux boundary improves convergence of a physical state.
7. Re-run the existing low-He/H molecular regression because the solver is a
   shared module. Atomic-only cases do not need to be regenerated unless the
   implementation changes their branch.

The He/H window and any production mass-loss rate should be remeasured only
after these conditions pass. Values from the present He/H=1 molecular runs are
diagnostics of an invalid state, not physical predictions.

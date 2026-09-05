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

**Status, 2026-09-02: the window is closed and it was not physics.** Section 12
identifies the two code defects that produced it -- an on-the-spot He
recombination photon budget with no competing absorber, and an initial
condition that clipped the H2 fraction at 1 instead of at the element-ratio
ceiling -- and shows the He/H = 1 arm relaxing to a subsonic, warm, fully
molecular base once either one is cut. Read section 12 first; sections 1-11
are the investigation that led to it and their conclusions about *causes* are
superseded there, while their measurements of the boundary condition itself
stand.

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
term. This is the answer to item (P)'s question "should the base BC be able to
reach a supersonic inflow at all": it should not. Nothing in the code tested
for that crossing when this was written; the base-face Mach diagnostic of
11.5-D now does, and reports the maximum inflow Mach of every run. **What has
changed since is that the arm no longer produces a crossing at all**: section
12 shows the supersonic base was driven by two code defects, and with those
corrected the same configuration keeps the base at Mach 0.086. The
characteristic count below is still the reason a crossing would be
inadmissible; it is no longer a description of what this configuration does.

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

**Superseded by section 12 (2026-09-02).** Section 11 was right that the
chemistry solver must not be excluded and that the retained state of the first
physical cell is where the failure starts. It was wrong about *why* that state
is retained: the fully ionized cell is an exact root of the hydrogen balance
row, not a non-root that slipped through acceptance, because the H I
photoionization rate fed into that row had been inflated to 1e+94 s^-1 by a
division by a vanishing H I density. There is no composition window to
compute: at He/H = 1 the corrected code gives a subsonic molecular base.

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

---

## 12. The window was two code defects

**Judgment first: there is no physical window. "A molecular base runs
supersonic near He/H ~ 1-10" is not a property of helium-rich molecular
atmospheres; it is what two defects in this code do to each other, and with
both corrected the same arm relaxes to a subsonic, warm, fully molecular base.
Nothing about the boundary condition needed to change to get there.** The
characteristic count of section 2 stands as written -- a supersonic inflow
would still be an over-specified boundary -- but the arm no longer produces
one, so the crossing that section 11 was trying to explain never happens.

The two defects are independent of each other, they are in different files,
and cutting either one removes the runaway. That is the strongest statement
this document can make about causation and it is what the controlled runs
below show.

### 12.1 Defect 1: He recombination photons with no competing absorber

`he_rec_coupling` (`src/modules/radiation/util_ion_eq.f90`) implements the
Draine (2011) on-the-spot treatment of the photons emitted when He II
recombines: a photon above 13.6 eV emitted inside the cell is assumed to be
absorbed inside the cell, and it ionizes hydrogen. The rate it returns is the
photon production rate divided by the H I density, because a rate per H I atom
is what the ionization balance wants:

    dP_HI(j) = P_add/max(nhi(j), 1.0d-99)

One channel of that sum -- the >= 24.6 eV ground-capture continuum -- carried a
competition factor,

    y = 1/(1 + R n_HeI/n_HI),    R = sigma_HeI(24.6)/sigma_HI(24.6) ~ 6,

which is the fraction of those photons that H I wins from He I. It vanishes
correctly as `n_HI -> 0`. **Every other channel had no competition factor at
all**: in the He 2^3S branch, which is the one the arm runs, those are the
singlet-excited captures (the 584 A resonance and the 2^1S two-photon
continuum), the 19.8 eV 2^3S decay line, and the 2^3S collisionally converted
into the singlets; in the atomic branch it is the case-B cascade. All were
assumed to deposit *every* photon into H I. With `n_HI = 0` the floor `1.0d-99` then converts a perfectly ordinary
photon production rate into an enormous rate coefficient. In the He/H = 1 arm
`P_add ~ 1.5e-5 cm^-3 s^-1` became `dP_HI ~ 1.5e94 s^-1`, and `x(H II) = 1`
became the exact root of the hydrogen balance row of the first physical cell.

Three things follow, and each of them is why this was hard to see:

* the state is a *root*, not a failed solve. The residual test added in
  section 113 passes it, because the normalizing turnover of the hydrogen row
  (`s(1)` in `System_HeH_mol.f90`) is built from the same contaminated `P_HI`.
  The normalized residual is 2.7e-17 while the raw residual is -1.3e+16;
* the heating of the same channels is a volumetric rate, `n_HeII n_e alpha`
  times an energy, and stays finite (2.6e-10 erg cm^-3 s^-1). The rate and the
  heat of one channel were structurally inconsistent, which is the signature
  to look for;
* it needs `n_HI` to reach *exactly* zero. Defect 2 is what produces that.

**Why this is wrong physics and not a missing clamp.** In a molecular layer the
species that absorbs a 19.8-24.6 eV photon is H2, whose photoionization cross
section in that band is 1.2 to 3.7 times that of H I (`sigma_H2`, Yan,
Sadeghpour & Dalgarno 1998, already in `cross_sec.f90`); in a helium-rich gas
the >= 24.6 eV photons go to He I. The photon must be divided among the
absorbers in the ratio of their absorption coefficients, and the share H2 takes
must arrive as **H2 photoionization**, not as hydrogen ionization and not as
nothing. Clamping `dP_HI` would have hidden the runaway and kept giving H2's
photons to H I.

The routine now computes, for every channel c at its own photon energy E_c,

    w_s(E_c) = n_s sigma_s(E_c) / sum_s' n_s' sigma_s'(E_c),

over s in {H I, H2} and, for the >= 24.6 eV channel, He I; returns `dP_HI` and
a new `dP_H2` from the H I and H2 shares; excludes both from the effective He II
recombination coefficient (`alpha_B + (w_HI + w_H2) alpha_1`); and deposits the
photoelectron of each share at its own threshold. Where no absorber is present
the photon leaves the cell and every returned rate is zero. Measured in the
molecular base of `mol_diffusion`, where `n_H2/n_HI = 1.7e3`, H2 takes 99.98%
of these photons; the earlier code gave all of them to H I.

### 12.2 Defect 2: the IC clipped x2 at 1 instead of at the element ceiling

`set_IC.f90` seeded the molecular column from the chemical-equilibrium fit and
capped the result:

    x2_ic = 2 q (1+He/H)/(1+q)
    if (x2_ic .gt. 1.0d0) x2_ic = 1.0d0
    dfHI_ic = min(x2_ic/mass_per_H, f_sp(j,isp_HI))

This is the same silent cap section 117 removed from the base particle count,
left behind in the initial condition. `q_h2_equilibrium` is a
solar-composition fit (item P26) and exceeds the attainable ceiling
`0.5/(0.5 + He/H)` for any He/H > 0.0964, so above that helium abundance the
cap engages, `x2_ic` is exactly 1, and `dfHI_ic` equals `f_sp(j,isp_HI)`. The
subtraction on the next line then leaves `f_sp(j,isp_HI) = 0` **exactly** --
the input defect 1 needs. At the He/H = 0.0793 of the five molecular
regression cases `x2_ic = 0.984 < 1`, the cap never engages, and none of them
can show the failure.

The IC now takes its H2 fraction from the same source the inflowing ghost
does: `base_h2_nuclei_fraction()` when a lower-atmosphere handoff states the
partition, and otherwise the fit clipped **at the ceiling** with a one-line
startup notice. Clipping at the ceiling still gives `x2 = 1`, fully molecular
hydrogen, which is the physical statement the ceiling makes -- and with defect 1
corrected that is harmless.

### 12.3 The controlled runs

All rows are `backup/regression/armD_D2` (the hot Uranus at He/H = 1 with
`q_H2_base 0.300861`, i.e. the x2-matched arm of section 117), 12000 steps,
`OMP_NUM_THREADS=1`.

| run | change | max base Mach | energy-floor activations | log10 Mdot |
|---|---|---:|---:|---:|
| baseline (stored `armD_D2/run.log`) | as it stood | 73.09, supersonic on 11994 of 12000 steps, first at step 6 | 224494 in 335 cells | 8.88 |
| control A | `He_rec_coupling: False` (defect 1 disabled) | 0.0876 | 0 | 10.89 |
| control B | IC cap moved to the ceiling (defect 2 fixed) | 0.0864 | 0 | 10.91 |
| both corrected | the code as it now stands | 0.0864 | 0 | 10.91 |

Each control cuts a different link and they reach the same solution, which is
what "two defects feeding each other" predicts. In the corrected run the first
physical cell at `r = 1.000193 R_p` carries
`n(H2) = 2.49e13`, `n(H I) = 2.28e10`, `n(H II) = 1.4e3 cm^-3` at 1165.7 K:
warm, subsonic, and molecular, against `n(H I) = 0`, `n(H II) = 5.0e13`,
132.1 K in the D2 row of section 11.2.

### 12.4 Why D4a escaped

D4a (`Base ghost temperature: continuous`) was the one variant of section 11.2
whose first cell was not fully ionized, and its first physical cell is why: it
retained `n(H I) = 3.96e9 cm^-3` against `n(H II) = 7.7e-4`, so the division
by `n_HI` never approached the `1.0d-99` floor and defect 1 never fired at
full strength. Its base failed for a different reason -- the continuous ghost
temperature put the first cell at 11242 K -- which is why its numbers look
unlike every other row of that table. D4a is not evidence against the
mechanism; it is the row in which the mechanism was switched off by an
accident of the state.

### 12.5 Status of the section 11.5 items

* **A (refuse non-roots immediately)** -- done in section 113, and it could not
  have caught this. The acceptance test is a *normalized* residual, and the
  normalization is the row turnover built from the same `P_HI` the defect
  inflated. Both numerator and denominator carried the 1e+99, so the test
  reported 2.7e-17 on a state whose raw residual was -1.3e+16. A residual test
  cannot detect a defect in the coefficient it normalizes by.
* **B (constrained solve rather than projection)** -- done in section 114. Same
  remark: the constrained solver solves the equations it is given, and it was
  given the inflated rate.
* **C, E (composition boundary, validate the base H2 input)** -- done in
  section 117. C is what makes the corrected IC of 12.2 possible: without a
  single definition of the base H2 fraction there is nothing for the IC to
  agree with.
* **D (Mach diagnostic)** -- done; it is the `base boundary:` line quoted in
  12.3, and it is now reporting a subsonic base rather than warning.

The other half of C, the one this document asked for as an alternative --
*"If the intended quantity remains EOS-only, the code must instead refuse a
mismatch between that EOS state and the species state used at the same ghost"*
-- is **also done, 2026-09-02 (item P35)**. With `Molecular base: True` and
`Molecular chemistry: False` the base particle count was molecular while the
species state stayed atomic, and the `ntot_bc` recomputation that reconciles
them is gated on `thereis_mol`, so the isothermal base boundary silently sat at
`ntot_bc x T0`. Measured on `backup/regression/armD_D1`: 876.34 K where
`Equilibrium temperature: 1140.0` was asked for, exactly `0.768722 x 1140`.
`input_read` now stops on the combination and prints `q_H2`, `x2`, `ntot_bc`,
the implied ghost temperature and the two ways to fix it.

The refusal invalidates 31 existing run directories, which is why it was put to
the user before being written: `armD_D1` itself, and 30 HD 209458 b VULCAN
handoff runs under `vulcan_work/` whose base ghosts sat between 0.555 and 0.994
of their requested 1450 K. None of them has been re-run. Each carries an
`INVALID_BASE_TEMPERATURE.md` with its own `ntot_bc` and ghost temperature, the
two parent directories carry a summary table, and the documents that quote them
now carry a stale marker at the quoting passage.

The boundary condition itself was not changed, and on this evidence it did not
need to be. What section 11.5 called "the boundary-condition problem" in item
P1 was the *downstream* half of a chain whose source was the photon budget of
one cell.

### 12.6 The arm converges: the first quotable He/H = 1 rate

`backup/regression/arm_heh1_x2matched` -- the same configuration without the
12000-step cap, PLM+WENO3, `du_th 0.5 1e-3`, `Solver: Newton`,
`OMP_NUM_THREADS=1` -- was run to see whether the corrected arm would reach the
`du` threshold on its own and hand over to the JFNK Newton finish. Marching
alone does not get there in any reasonable time (below), but the state it is
descending towards **does** solve under the JFNK when the solver is handed it
directly, and that is where the rate below comes from.

The run is well behaved and it is converging, just slowly. It crosses
PLM -> WENO3 at step 17632, and from about step 26000 onward both the mass-flux
spread and the steady residual fall monotonically with no oscillation:

| step | flux spread | steady residual (ref) |
|---:|---:|---:|
| 26000  | 4.32e-1 | 1.156 |
| 45500  | 8.76e-2 | 4.33e-2 |
| 106500 | 3.07e-2 | 1.61e-2 |
| 152500 | 1.96e-2 | 1.04e-2 |
| 244500 | 1.03e-2 | 5.58e-3 |

Over that range the decay is geometric and its rate is itself falling -- about
0.8% per 1000 steps at step 1.5e5 and 0.6% per 1000 at step 2.4e5 -- which put
the 1e-3 threshold of order 4e5 further steps away.

**It does not get there: the marching turns around (measured 2026-09-02, run
still going).** `du` reaches its minimum **9.96e-3 at step 251,287** and then
rises again -- 0.737 at step 260,000, 7.1e-2 at 270,000, 6.6e-2 at step
272,816 -- so this trajectory does not reach `du_th = 1e-3` at all, and the
extrapolation above is superseded by the measurement. That is the `du` floor
item (S) of `TO_BE_DONE.md` records for the molecular A2 wind, here met under
WENO3 rather than PLM, and it is a separate problem from the one this document
is about: the state being descended towards is subsonic, warm and molecular at
the base, with zero energy-floor activations and no non-root acceptances.
Nothing in it resembles the failure of sections 1-11.

**The arm does not have to march there: the JFNK solves it (2026-09-02).** The
step-168,550 state of that run (`du = 1.735e-2`, still falling monotonically)
was taken as a restart -- `Load IC? True`, `du_th [PLM,WENO3]: 1.0e9 1.0e-3`,
`Solver: Newton 100.0`, so the hand-off is armed and the solver is entered
without waiting for the `du` floor. It marches 2002 steps to `du = 3.1094e-2`,
enters the JFNK **three times, every one of them returning `info = 0`**, and
closes as *converged: JFNK steady solution* at `||R|| = 6.572e-4`, with

> **`log10 Mdot = 10.32` g/s**, and the same at 1 and at 8 threads

**This is the first quotable quantitative He/H = 1 value.** Everything this
document quoted for the arm before it -- including the 12000-step snapshot
10.91 of section 12.3 -- is a relaxation snapshot, and item (P17) of
`TO_BE_DONE.md` says why such a number is not a rate.

The solve is clean by every counter the code keeps: over the whole run the
marching phase accepts 1,011,494 converged roots, 16 roots without solver
convergence and 18 constrained-continuation roots and **no non-root state at
all**, with no non-finite reaction residual and no cell left outside the
element simplex; the steady solver's own 26 iterate sweeps and 490 probe sweeps
are equally clean, and not one admissibility event of section 121 of
`Update_EXHALE` fires (no truncated GMRES cycle, no zeroed Jacobian color, no
rejected trial). The reference state is kept at
`.../scratchpad/jfnkprobe/fixsched/run/heh1_jfnk_t1/`.

**The structure it converges to, three layers.** Read from that run's
`output/Hydro_ioniz.txt` and `output/Ion_species.txt`. The grid is
**500 physical cells from `r = 1.000193` to `4.7253 R_p`** (corrected
2026-09-02: the earlier "1.0000 .. 4.814 R_p in 503 physical cells" counted
the two outer ghost rows and the inner ghost at `r = 1.0000`, which the output
file carries alongside the physical cells; the setup report says
`Grid cells: 500 computational cells (+ 2 ghost cells on each side)`). The
outer-boundary values below are quoted on the physical range.

| layer | radial range | what marks its end | T [K] | n_e [cm^-3] |
|---|---|---|---:|---:|
| molecular  | 1.0000 - 1.0468 | `f_H2 = 0.5` at 1.0468, T = 623 | 1140 -> 623 | 2.7e5 -> 3.5e7 |
| atomic     | 1.0468 - 1.5859 | `n(H II) = n(H I)` at 1.5859, T = 4827 | 623 -> 4827 | 3.5e7 -> 7.2e7 |
| ionized    | 1.5859 - 4.7253 | outer boundary, `x_HII = 0.948` | 4827 -> 2146 | 7.2e7 -> 3.9e6 |

The base is `f_H2 = 0.9251` of the H nuclei, as the handoff states it. The H3+
peak sits at the very bottom of the molecular layer, `r = 1.0014 R_p`,
`n(H3+) = 3.16e5 cm^-3` at 442 K -- well below the H2 front, not at it. The
temperature minimum of the whole domain, 440 K, sits at `r = 1.0665` just above
the H2 front; the maximum, 4967 K, is at `r = 1.4207` in the atomic layer, so
the H ionization front sits just outside the temperature peak, and the gas
cools again to 2146 K at the outer boundary. `n_e` spans 2.7e5 to 1.4e8 cm^-3
over the domain, and the outflow leaves the outer boundary at 7.34 km/s.

**The marching run is still going, as a cross-check** that the state the JFNK
solved is the one marching was descending towards; its output is at
`.../scratchpad/hrcfix/runs/arm_heh1_x2matched/`. Note that the restart was
taken at step 168,550, i.e. **before** the turnaround above, so the JFNK solved
the state on the descending branch; what the branch after step 251,287 settles
to is a separate question and is not answered here.

**The thread dependence this arm used to show was a defect of the solution's
definition, and it is repaired (item (P37), 2026-09-02).** Before the repair
the same case converged to `log10 Mdot` 10.32, 10.33 and 10.35 depending on the
thread count and the OpenMP schedule. The cause was in the constrained
continuation: its SEED floored the species of `n_unknown` /
`species_of_unknown`, which are rung-partition state built later in the solve,
so on entry to a cell they still described whichever cell that thread had
solved last. The seed now uses the cell's own species set, and the rung context
is cleared on entry. The case is byte-identical at 1 thread and at 8 threads
under `static`, `static,1` and `dynamic,8`, and so are `mol_sec_ion` and
`mol_metals`. **The numbers in this section are from the repaired code**; the
serial answer moved too, because it was order-contaminated in the same way, and
the cells that moved are exactly the constrained-continuation (class-5) cells.


## 13. The ladder, converged

**Judgment first. Every rung of the He/H ladder now has a Newton-converged
solution, the JFNK returns `info = 0` on all of them, and the ladder is
monotonic and undramatic: raising He/H from 0.0793 to 30 thins the molecular
layer, pulls the hydrogen ionization front inward, heats the wind, and leaves
the mass-loss rate almost unchanged. There is no window and no crossing
anywhere on the ladder, which is what section 12 predicted. The JFNK solution
of the He/H = 1 arm is also stable under pure time marching -- 50,000 marching
steps leave it where it was and lower its steady residual -- so on this
evidence the arm has one solution, not two. What the campaign did turn up is a
new code defect, in the same He recombination coupling section 120 repaired,
which fires only for He/H >= 3 and only in the fully ionized outer wind; it
moves the ionization front of those rungs and does not touch their molecular
layer or their rate.**

Measured 2026-09-02 on a build made from a copy of the working tree taken at
14:24 KST -- the tree as it then stood, with the section-117 to section-121
changes in it and uncommitted (`seed_species_of_cell` is present in
`constrained_chemical_equilibrium.f90`, so the section-121 constrained-solver
seed repair is in; HEAD was `6d07d48`). The copy was built and run outside the
repository so that a concurrent change to the tree could not enter the
campaign, and none of these arms sets `Stellar LW flux`, which is what that
change touched. The runs are in the campaign scratch directory under
`.../scratchpad/ladder/runs/`. No metals, no `opacity.inp`;
`Secondary_ionization` at its STAGED default.

### 13.1 How each rung was converged

Each rung is `arm_heh1_x2matched` with `He/H number ratio` and the `base.inp`
`HeH_base` / `q_H2_base` changed, and nothing else. `q_H2_base` is fixed by
holding the H2 nucleus fraction at `x2 = 0.9251143` -- the value the molecular
matrix cases carry -- so that the ladder varies He/H alone:

| He/H | `q_H2_base` | x2 | ceiling `0.5/(0.5+He/H)` | q/ceiling |
|---:|---|---:|---:|---:|
| 0.0793 | 0.75 | 0.9251143 | 0.8631107 | 0.869 |
| 0.3 | 0.552344723 | 0.9251143 | 0.6250000 | 0.884 |
| 1 | 0.300861 | 0.9251134 | 0.3333333 | 0.903 |
| 3 | 0.130760315 | 0.9251143 | 0.1428571 | 0.915 |
| 10 | 0.0438965268 | 0.9251143 | 0.0476190 | 0.922 |
| 30 | 0.0151472127 | 0.9251143 | 0.0163934 | 0.924 |

Every value is admissible: `q_H2_base` is below its element-ratio ceiling on
every rung, by 8 to 13%.

The He/H = 0.0793 rung is `backup/regression/mol_base_handoff` without its
`maxsteps` pin, which differs from the arm in one further respect -- it does
not set `He_diffusion` -- so a **second 0.0793 rung was run with the arm's
`He_diffusion: True` and `He_Kzz: 1.0e9`** to keep that variable from riding
along with He/H. It is reported as its own row.

The recipe of section 12.6 was used with one change, and the change is the
reason all six converged. Each rung was first marched with
`du_th [PLM,WENO3]: 0.5 1.0e-3` and `Solver: Newton 2.0e-2`, a step budget of
250,000 (`EXHALE_MAXSTEPS`) and `OMP_NUM_THREADS=8`. Every rung reaches
`du = 2e-2` on its own, and the staged secondary ionization fires there first
(`secondary ionization activated ahead of the Newton finish`, at step 37,857 to
227,624 depending on the rung), so the wind re-relaxes under the full physics
and has to cross the switch a second time. Three rungs (0.0793 with and without
diffusion, and 0.3) got there inside the budget and **the JFNK entered at that
second crossing failed on all three, `info = 2` at `||R|| = 1.253` to
`1.258e-2`**, after which the code falls back to marching as it is written to.
The other three (3, 10, 30) never made a second crossing inside 250,000 steps.
The two 0.0793 marching runs were stopped by hand at steps 198,638 and 199,597
once their in-run attempt had failed and the restart below had been started
from their output; the other four ran the full 250,000.

**What converges all six is the section-12.6 restart.** Take the state the
marching left, restart it with `Load IC? True`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-3` and `Solver: Newton 100.0` -- which arms the
hand-off immediately, so the solver is given the state rather than made to
march to it -- and the run marches 2002 steps (the secondary-ionization hold)
and then enters the JFNK. Five of the six rungs converged on the FIRST such
restart: 0.0793-with-diffusion from its step 199,597 state, and 0.3, 3, 10 and
30 from their step 250,000 states. The 0.0793 rung without diffusion needed
TWO: the first restart, from its step 198,638 state, reached
`||R|| = 1.440e-3` and still returned `info = 2`, then marched itself down to
`du = 2.9e-3` over 20,000 steps, and the restart from THAT converged. Each
restart takes under two minutes; the marching that precedes it takes hours.

### 13.2 The rungs

Seven `info = 0` states, `OMP_NUM_THREADS=8`: the six converged here, plus the
He/H = 1 reference state of section 12.6 re-read on the physical cell range.
`n(JFNK)` is the number of steady solves each run made, every one of them
`info = 0`.

| He/H | 0.0793 | 0.0793 (diffusion) | 0.3 | 1 | 3 | 10 | 30 |
|---|---:|---:|---:|---:|---:|---:|---:|
| n(JFNK), all info=0 | 1 | 2 | 3 | 3 | 1 | 1 | 1 |
| final residual (JFNK norm) | 9.89e-4 | 9.16e-4 | 7.97e-4 | 6.57e-4 | 6.02e-4 | 9.21e-4 | 4.95e-4 |
| `du` of the solution | 1.61e-3 | 1.50e-3 | 1.06e-3 | 3.71e-4 | 1.04e-3 | 1.72e-3 | 8.66e-4 |
| **log10 Mdot [g/s]** | **10.42** | **10.49** | **10.35** | **10.32** | **10.36** | **10.42** | **10.39** |
| **log10 Mdot, whole-column measure** | **10.28** | -- | -- | **10.29** | -- | -- | **10.37** |
| flux spread `r > 1.2 R_p`, old / new measure | 9.27e-3 / 2.39e-3 | -- | -- | 2.58e-3 / 2.59e-3 | -- | -- | 1.19e-2 / 2.86e-3 |
| H2 front `f_H2 = 0.5` [R_p] | 1.1012 | 1.1200 | 1.0681 | 1.0469 | 1.0467 | 1.0357 | 1.0250 |
| T at the H2 front [K] | 765 | 974 | 626 | 622 | 1140 | 1496 | 1287 |
| `f_H2` down 4 decades by [R_p] | 1.273 | 1.354 | 1.154 | 1.100 | 1.068 | 1.057 | 1.044 |
| H3+ peak [R_p] | 1.0226 | 1.0200 | 1.0167 | 1.0014 | 1.0017 | 1.0015 | 1.0027 |
| n(H3+) at the peak [cm^-3] | 3.09e5 | 2.97e5 | 3.20e5 | 3.16e5 | 2.42e5 | 1.95e5 | 1.48e5 |
| T at the H3+ peak [K] | 922 | 1110 | 796 | 442 | 930 | 980 | 842 |
| n(H3+) > 1 cm^-3 out to [R_p] | 1.181 | 1.177 | 1.118 | 1.081 | 1.058 | 1.049 | 1.035 |
| ionization front `x_HII = 0.5` [R_p] | 2.573 | 2.755 | 2.040 | 1.593 | 1.396* | 1.333* | 1.241* |
| T minimum [K] | 502 | 671 | 433 | 440 | 758 | 979 | 842 |
| T maximum [K] (at [R_p]) | 2751 (1.815) | 2698 (1.967) | 3642 (1.551) | 4967 (1.421) | 5864 (1.464) | 6328 (1.558) | 6660 (1.544) |
| n_e at the base [cm^-3] | 2.95e5 | 2.92e5 | 2.90e5 | 2.72e5 | 2.59e5 | 2.77e5 | 3.92e5 |
| n_e maximum [cm^-3] | 7.77e7 | 7.14e7 | 1.02e8 | 1.39e8 | 1.60e8 | 1.67e8 | 1.70e8 |
| n_e at the outer boundary [cm^-3] | 6.80e6 | 7.87e6 | 5.02e6 | 3.90e6 | 3.64e6 | 3.73e6 | 3.43e6 |
| `x_HII` at the outer boundary | 0.794 | 0.774 | 0.869 | 0.948 | 1.000* | 1.000* | 1.000* |
| outflow speed at the boundary [km/s] | 8.97 | 8.97 | 8.28 | 7.34 | 6.85 | 6.74 | 6.61 |

\* contaminated by the defect of section 13.4; the coupling-off control values
are 1.646, 1.593 and 1.506 R_p for the front and 0.908, 0.910 and 0.917 for
`x_HII`.

**The two added rows are a 2026-09-02 re-measurement, and where they carry a
number they supersede the row above them.** Every rate in the `log10 Mdot` row
was obtained with a convergence measure that looked only at `r > r_esc`, i.e.
at 109 of the 500 cells, so the column between the base and 2 R_p was
determined by nothing; item (P42) and `Update_EXHALE.md` section 127 replace it
by the larger of the wind's residual and the below-escape layer's, each on its
own physical scale. Three rungs were re-run on the current tree under the new
measure, all `info = 0`. The signature of the repair is the flux-spread row:
under the new measure the spread of `rho v r^2` below 2 R_p equals the spread
above it on all three rungs, where before the inner column was 4 to 5 times
worse on two of them. The rates move by at most 0.05 dex; the structure below
`r_esc` moves by 4 to 27% in the L2 norm of (rho, v, p). The remaining four
rungs have not been re-run and their entries in the `log10 Mdot` row are still
old-measure values. Part of the 0.05-0.14 dex difference between the two rate
rows is not the measure at all: the re-runs are on the current tree, which
carries sections 125 and the `q_H2_base` ghost-composition change, and the same
three rungs re-run on the current tree under the OLD measure give 10.33, 10.29
and 10.37.

Read across, the ladder is monotonic in everything that describes the
structure. **More helium means a thinner molecular layer** (the H2 front moves
from 1.101 to 1.025 R_p), **a hotter wind** (T maximum 2751 to 6660 K), **an
ionization front closer to the planet** (2.573 to 1.241 R_p) and **a slower
outflow** (8.97 to 6.61 km/s). The trend has the sign helium's extra opacity
and its 24.6 eV photoelectrons predict: each H nucleus is accompanied by more
absorbing and heating helium, the layer heats, and the H2 front is pushed down.

**The mass-loss rate is the exception: it is flat.** `log10 Mdot` spans 10.32
to 10.49 over a factor 378 in He/H, and it is not even monotonic -- it has a
shallow minimum at He/H = 1. The 0.17 dex spread includes the 0.07 dex the
`He_diffusion` switch alone moves the 0.0793 rung, so the He/H dependence of
the rate over this range is at most a factor 1.5 and probably less.

**He diffusion is a small effect at He/H = 0.0793.** With it on the H2 front
moves outward by 0.019 R_p, the ionization front by 0.18 R_p, and `log10 Mdot`
from 10.42 to 10.49. It is not zero, and the two 0.0793 rows should not be read
as one.

**Ledger.** None of the six converging runs admitted any base inflow at all --
the `base boundary:` line, which prints whenever a run has one, is absent from
all six logs (the He/H = 1 reference of section 12.6 did have one, subsonic at
max Mach 1.067e-2) -- and none activated the energy floor. Root acceptance is
clean for He/H <= 10: zero non-root acceptances, zero non-finite reaction
residuals, zero projected/handback roots. **He/H = 30 is not clean**: 34
non-root states accepted under the relaxation amnesty (12 in the marching
phase, 22 in the steady solver's probe sweeps), 153 + 6 + 52
projected/handback roots, and 165 cells with no in-simplex starting point.
Every one of the printed non-root cells sits at `r = 2.12 - 2.22 R_p` at
5700-5820 K -- the ionized outer wind, which is where section 13.4's defect
lives. The largest element-budget violation grows along the ladder, with one
step out of order: none reported at 0.0793, 1.2e-3 at 0.3, 2.1e-1 at 1,
1.5e-1 at 3, 2.8 at 10, 1.1e1 at 30.

**Quotable.** Item (P17) may be lifted for the six converged states listed
above, with the asterisked entries excluded and the He/H = 30 row carrying its
ledger.

### 13.3 The converged solution is stable under time marching

**Judgment: it stays.** The reference state of section 12.6 was restarted with
`Load IC? True`, `Secondary_ionization: Immediate` (so the coupling is on from
step 0, as it is in the state being tested), `du_th [PLM,WENO3]: 1.0e9 1.0e-3`
and `Solver: Newton 1.0e-30`, which selects the marching branch and puts the
JFNK hand-off out of reach, and marched **50,000 steps with no steady solve at
all**.

| step | `du` | `dtu` | steady residual (ref) |
|---:|---:|---:|---:|
| 2 | 3.81e-4 | 9.2e-5 | 6.01e-4 (step 500) |
| 1,035 | **2.582e-4** (minimum) | | |
| 3,968 | **5.094e-3** (maximum) | 3.8e-5 | 1.12e-2 (step 5,500) |
| 10,002 | 6.65e-4 | 2.0e-6 | 2.40e-3 (step 10,500) |
| 20,002 | 2.60e-4 | 8.4e-8 | 3.32e-4 (step 20,500) |
| 35,002 | 3.88e-4 | 5.3e-8 | 3.99e-4 (step 35,500) |
| 50,000 | 2.625e-4 | | **2.070e-4** |

The trajectory has one excursion and then settles. `du` sits at its minimum for
about 3000 steps, rises by a factor 20 over the next 1000, decays over the next
15,000, and then holds between 2.58e-4 and 2.63e-4 for the last 30,000 steps
(mean 2.603e-4) with `dtu` at 3e-8. **The steady residual ends below where the
JFNK left it**: 2.07e-4 against the solver's own 6.572e-4. `log10 Mdot` is
10.32 at step 50,000, the same value to two decimals. Handing the marched state
back to the JFNK closes it again immediately: `info = 0`, `||R|| = 3.532e-4`,
`log10 Mdot = 10.32`.

**Where the excursion lives.** Comparing the periodic output at steps ~3001 and
~4001, the cells whose mass flux moves are `r = 2.12 - 2.19 R_p` -- inside the
`[j_min:N]` window over which `du` is measured, which begins at the escape
radius `r = 2.0 R_p` -- and they move by 0.3%. It is not a base event: over the
same interval density, temperature, pressure and every species change by less
than 0.5% below 1.1 R_p. The one quantity that changes by order unity anywhere
is the velocity in the base cells `r = 1.006 - 1.017`, where it changes sign
and `|v|` is of order 1 cm/s against 7.5 km/s in the wind. That is the base
breathing this document has described throughout; `du` does not see it, because
`du` is measured above the escape radius.

**Where the marched state differs from the JFNK state, after 50,000 steps.**
The differences are confined to the base and to the temperature minimum:
`|dT/T|` reaches 5.8% at `r = 1.0015` and about 2.1% over `r = 1.075 - 1.083`;
density and n(He I) move 6.2% at `r = 1.0015`; n(H II), n(He II) and n(HeH+)
move 32-34% at the same cell, where they are trace (n(H II) goes from 58 to
76 cm^-3 against n(H2) = 2.5e13). Above `r = 2 R_p` nothing moves at the
precision of the rate.

**What this answers, and what it does not.** The two-solution question of item
(P17) asked whether the marching branch and the JFNK branch are two different
steady solutions. On this arm they are not: marching does not leave the JFNK
solution, it improves its residual, and it returns the same rate; the two
states differ only in a base layer whose velocity is four orders below the
wind. **It does not show that no second solution exists elsewhere** -- it shows
that this one is stable to this perturbation. In particular the marching branch
of section 12.6 after its turnaround at step 251,287 was not tested here, and
whatever it settles to is still unmeasured.

Run ledger for the 50,000 steps: 25,199,496 converged roots, no root without
solver convergence, no projected/handback root, no non-root acceptance, no
energy-floor activation, no base inflow.

### 13.4 A defect the ladder found: the He recombination coupling divides by a vanishing absorber in the ionized wind

**Judgment: physically wrong, and it is the same failure mode section 120
repaired at the base, surviving in the one place that section did not
reach -- a cell with no absorber for the sub-24.6 eV channels at all.**

The diagnostic section 120 added, `n_cells_he_rec_photoionization_dominant`
(cell-steps in which the coupling's `dP_HI` exceeds the direct stellar
photoionization rate by more than 1000), was reported there as zero on every
case measured. **It is not zero on the helium-rich half of this ladder**, and
the threshold is sharp:

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 |
|---|---:|---:|---:|---:|---:|---:|
| cell-steps flagged (converging restart) | 0 | 0 | 0 | 1.50e5 | 2.35e5 | 3.03e5 |
| largest `dP_HI/P_HI` | -- | -- | -- | 4.25e105 | 4.20e106 | 5.60e106 |

Instrumenting the routine (a diagnostic build of the same source in
`.../scratchpad/ladder/diag/`, printing the arguments each time a new maximum
is set) locates it exactly. In the He/H = 30 solution:

| j | r [R_p] | `dP_HI` [s^-1] | `P_HI` [s^-1] | n(H I) | n(H2) | n(He I) |
|---:|---:|---:|---:|---:|---:|---:|
| 348 | 1.5859 | 6.21e5 | 2.40e-5 | 1.70e-3 | 3.99e-27 | 1.57e8 |
| 354 | 1.6305 | 1.39e12 | 2.41e-5 | 6.65e-10 | 1.44e-23 | 1.26e8 |
| 355 | 1.6383 | 1.50e27 | 2.41e-5 | 0 | 3.81e-25 | 1.21e8 |
| 460 | 3.2929 | 2.60e100 | 2.51e-5 | 0 | 2.10e-151 | 2.04e6 |

`P_HI` is perfectly ordinary throughout, so the ratio is not a denominator
artifact: it is `dP_HI` that diverges, exactly as in section 120.

**Why the section-120 form does not protect this case.** The correction writes
the H I share as `w_HI = 1/(1 + R_He n_HeI/n_HI + R_H2 n_H2/n_HI)` and the rate
as `P_channel * w_HI / max(n_HI, 1e-99)`, so the returned rate is really
`P_channel/(n_HI + R_He n_HeI + R_H2 n_H2)` and stays finite **as long as one
absorber remains**. For the >= 24.6 eV ground-capture continuum one always
does here, because He I is abundant; that channel is fine. But the other four
channels -- the 584 A resonance at 21.2 eV, the 19.8 eV 2^3S line, the 2^1S
two-photon continuum at 16.1 eV, and the singlet-excited captures -- are
**below the He I threshold**, so He I is correctly not among their absorbers,
and their denominator is only `n_HI + R_H2 n_H2`. In the fully ionized outer
wind of a helium-rich rung n(H I) is exactly zero and n(H2) is numerical
residue, 1e-25 falling to 1e-151. The guard that is supposed to catch this,
`absorb_low = (nhi(j) + n2) .gt. 0.0d0`, tests for *exactly* zero, and 1e-151
passes it. The code then divides one cell's photon production by 1e-151.

**The physics the code should state instead** is the one section 120 already
wrote down for a cell with no absorber: the photon leaves and the rate is zero.
The on-the-spot assumption is a statement that the cell is optically thick to
these photons, and a cell holding 1e-151 H2 molecules and no hydrogen is
transparent to them. Section 120 recorded the untested optical depth as an
approximation; what was not recorded is that the untested case does not merely
overestimate the rate, it produces a rate coefficient of 1e100 s^-1.

**What it does to the solutions, measured.** Four controls with
`He_rec_coupling: False`, each restarted from the corresponding converged rung
and re-converged (`info = 0` in all four):

| He/H | log10 Mdot on / off | ionization front [R_p] on / off | `x_HII` at the boundary on / off | H2 front [R_p] on / off |
|---:|---|---|---|---|
| 0.3 | 10.35 / 10.31 | 2.040 / 2.119 | 0.869 / 0.852 | 1.0681 / 1.0628 |
| 3 | 10.36 / 10.34 | 1.396 / 1.646 | 1.000 / 0.908 | 1.0467 / 1.0453 |
| 10 | 10.42 / 10.40 | 1.333 / 1.593 | 1.000 / 0.910 | 1.0357 / 1.0354 |
| 30 | 10.39 / 10.36 | 1.241 / 1.506 | 1.000 / 0.917 | 1.0250 / 1.0249 |

At He/H = 0.3, where the counter never fires, switching the coupling off moves
the ionization front by 0.08 R_p and `x_HII` by 0.017: the ordinary size of a
real physical channel. At He/H >= 3, where it does fire, the front moves by
0.25 R_p and `x_HII` is driven to **exactly 1.000** with the coupling on. That
saturation is the signature: cells just inside the front, where n(H I) is small
but the gas is not yet fully ionized, are handed `dP_HI` between 6e5 and 1e12
against a stellar `P_HI` of 2.4e-5, and `x(H II) = 1` becomes the exact root --
which is section 120's mechanism, at the ionization front rather than the base.

**What it does not touch** is the molecular layer and the rate. The H2 front
moves by at most 0.0014 R_p between the two controls at He/H >= 3, and
`log10 Mdot` by at most 0.03 dex. The structural statements of section 13.2
below the H2 front, and the comparison of section 13.5, therefore stand as
measured; the ionization-front and outer-`x_HII` entries for He/H >= 3 do not,
and are asterisked.

**Not repaired here.** The source tree was in use by another change while this
campaign ran, so the fix was not written. It is recorded as its own item in
`TO_BE_DONE.md`. The shape of the fix is not in doubt -- the sub-24.6 eV
channels need an absorber test with a physical floor rather than `> 0.0`, and
below it the channel must return zero -- but choosing that floor is a physics
decision (it is the density at which the cell stops being optically thick to a
16-21 eV photon over its own width) and should be made deliberately.

**Fixed 2026-09-02, `Update_EXHALE.md` section 125**, and not with a floor: the
cell now keeps the fraction `1 - exp(-tau_c)` of each channel's photons and the
rest leave, which removes the divergence and the untested on-the-spot
assumption together (items P33 and P40 closed). The three contaminated rungs
were re-converged with the correction; the corrected ionization fronts are
1.618 (He/H = 3) and 1.444 R_p (He/H = 30) with `x_HII` 0.914 and 0.924, and
the counter is zero on both. The He/H = 1 rung, which never fired the counter,
moves as well -- front 1.586 to 1.758 R_p -- because the ordinary overestimate
of a thin cell is not a runaway and never tripped the threshold.

### 13.5 The P23 gate: compared against the published hot Uranus and super-Earth

**Judgment: the gate does not pass, and the explanation P23 proposed for the
disagreement is refuted by this measurement.** P23's step (a) was to converge
our metals-off molecular case and see whether H3+ becomes extended rather than
base-confined, on the reasoning that the difference from Frelikh &
Murray-Clay's H-and-He-only model is most likely the metal electrons our
metals-on gates carry. The metals-off case is now converged, and **H3+ is still
base-confined**: it falls from 3.09e5 cm^-3 at 1.023 R_p to 7.9e-4 cm^-3 at
1.5 R_p and 1.7e-5 at 2.0 R_p, nine decades over half a planetary radius.
Composition is therefore not the explanation, and by P23's own instruction the
next places to look are the planet parameters, the secondary-ionization gap and
the base composition.

The comparison is with Frelikh & Murray-Clay (2026, ApJ 996, 96), their Table 3
and sections 4 and 5.4, and with Koskinen et al. (2022, ApJ 929, 52) Model A,
their section 3.2.1 and Figures 7 and 8, both read from the publisher PDFs in
`references/`. **Koskinen's hot Uranus, not Frelikh's super-Earth, is our near
twin**, and it is worth writing down how near:

| | this work (He/H = 0.0793) | Koskinen+2022 Model A | Frelikh & Murray-Clay 2026 |
|---|---|---|---|
| planet | hot Uranus | Uranus-like | super-Earth |
| mass | 0.0457 M_J = 14.5 M_E | ~14.5 M_E | 7.2 M_E |
| radius | 0.49 R_J = 5.38 R_E (R_J = 6.9911e9 cm, the code constant; the escape base) | escape base at 1.34 R_p, R_p the 1 bar radius | 1.9 R_E |
| gravity at the base | 493 cm s^-2 | ~494 cm s^-2 (14.5 M_E at 1.34 x 4.007 R_E) | 1959 cm s^-2 |
| orbital distance | 0.048 au | 0.05 au | 0.05 au |
| base temperature | 1140 K | 1140 K | 1000 K |
| base pressure | 9.0e-6 bar | 1e-6 bar | 4.7e-7 bar |
| XUV flux at the planet | 1.56e3 erg cm^-2 s^-1 (13.6 eV - 1.24 keV) | 1.6 W m^-2 = 1.6e3 erg cm^-2 s^-1 | ~2.6e3 erg cm^-2 s^-1 (12-248 eV) |
| species | H, He (metals off) | H, He only | H, He only |
| base H2 | x2 = 0.925 (`q_H2_base = 0.75`, representative) | q(H2) = 0.84 at the lower boundary | H2 and He the major constituents |

Mass, gravity at the escape base, orbital distance, base temperature and XUV
flux all agree with Koskinen's Model A to within a few percent. The one input
that differs materially is the base level: ours is nine times deeper in
pressure, which if anything should give us *more* H2 at altitude, not less.

**Where we agree with them.**

* the **rate**: our `log10 Mdot = 10.42` g/s against Koskinen's
  1.9e7 kg s^-1 = 1.9e10 g s^-1, i.e. 10.28. A factor 1.4;
* the **outflow speed**: 8.97 km/s at our outer boundary (4.73 R_p) against
  their 9.8 km/s at theirs, and Frelikh's 5.8 km/s at 4.4 R_p;
* the **wind is mostly neutral**, which is the property both papers single out.
  Our `x_HII` is 0.078 at 1.5 R_p, 0.31 at 2.0 R_p and 0.73 at 4.0 R_p, so
  n(H I) exceeds n(H II) out to 2.57 R_p. Koskinen: "the density of neutral H
  is higher than the proton density at all radii". Frelikh: "The wind is mostly
  neutral H at the sonic point."

**Where we do not.**

* **H2 is truncated and theirs is not.** Our `f(H2)` falls from 0.995 at the
  base to 2.8e-5 at 1.5 R_p and 6.3e-6 at 2.0 R_p; four decades between 1.10
  and 1.27 R_p. Koskinen: "the density of H2 remains significant at all
  altitudes" and "the density of H2 is significant at all radii of the hot
  Uranus model". Frelikh: "The H2 is not cut off sharply as in the hot Jupiter
  case; its gradual drop off allows H3+ to be present throughout";
* **H3+ is base-confined and theirs is not.** Above; Koskinen report "the
  appearance of the molecular ions H3+ and HeH+ at high altitudes", Frelikh
  "the presence of H3+ throughout the outflow, instead of being confined to a
  narrow region at the base";
* **we are colder.** Our 502-2751 K against Koskinen's rise to about 4400 K at
  3 R_p and 5540 K at their upper boundary, and Frelikh's peak of 4200 K
  (their Figure 20 caption; their section 5.4 quotes the range as
  1000-3800 K).

**The mechanism Frelikh give for their result does not apply to ours, and that
is the sharpest thing this comparison says.** Their section 4 explains the
survival of H2 by temperature: "As the temperature threshold for significant
thermal dissociation is never reached in a super-Earth, molecular hydrogen is
present throughout, instead of being confined to a layer at the base." Our H2
front sits at **765 K**, and our whole domain peaks at 2751 K, against the
48,350 K barrier of the thermal-dissociation reaction R12 that both codes take
from Baulch et al. (1992). At 765 K that rate coefficient is `1.5e-9
exp(-48350/765)`, which is zero for any purpose. **Our H2 is not being
thermally dissociated; it is being removed by something else, at a temperature
where their argument says it should survive.** Identifying which channel is a
measurement this campaign did not make -- the candidates are H2
photoionization (the 15.4 eV channel, which both codes carry), the
H2+ -> H3+ -> dissociative-recombination route, and charge exchange with H+ --
and it is the obvious next step for this gate. Lyman-Werner photodissociation
is not a candidate: it is off in these runs, and Koskinen's Table 1 does not
carry it either. **That measurement has since been made; section 13.6 has the
answer, and it is not one of the three candidates listed here.**

**A caveat that belongs in the verdict.** Our base H2 partition is
`q_H2_base = 0.75`, which its own `base.inp` header calls a representative
photochemical value rather than a computed one, and P23 records that a
structural claim should not rest on a representative base. The claim made here
is comparative and one-sided: our base is *more* molecular than Koskinen's
(x2 = 0.925 against q(H2) = 0.84) and nine times deeper, so a base-composition
error would have to act in the direction of *suppressing* H2 at altitude to
explain the disagreement, and it is not obvious how it would.

**Verdict, stated as P23 asks.** Qualitatively **different**, and the
difference is **not explained** by the composition axis P23 proposed: our
metals-off, Newton-converged hot Uranus keeps a sharp H2 truncation and a
base-confined H3+ at a temperature far below thermal dissociation, on a planet
whose mass, base gravity, orbital distance, base temperature and XUV flux match
Koskinen's Model A. The rate, the outflow speed and the neutral-dominated wind
do agree. The gate does not pass.

The figure is `docs/lower_atmosphere_figs/fig_ladder_three_layer.png`:
`f(H2)`, `n(H3+)`, `n_e` and `T` against radius, one curve per rung. It is
drawn by `make_fig.py` in the campaign scratch directory beside the runs, from
the `output/` of the six converged rungs; it is not regenerated by
`docs/lower_atmosphere_figs/make_figures.py`.

### 13.6 What removes the H2

**Judgment: nothing in the chemistry is doing the work that section 13.5
attributed to it. The H2 truncation is very likely a property of the CLOSURE,
not of a reaction: EXHALE solves H2 in LOCAL chemical equilibrium, and in the
region where `f(H2)` falls its four decades the omitted advection term of the
H2 continuity equation is 4 to 500 times larger than the whole net chemical
rate.** Koskinen et al. (2022) say the same thing from the other side about
their own model: "The prevalence of neutral H and H2 at high radii in the hot
Uranus model is due to a combination of a relatively low temperature and high
escape flux that replenishes H2 and H at high altitudes where they would
otherwise be predominantly dissociated and ionized". Their H2-at-all-altitudes
is a transport result, and the state they say the gas would reach without the
transport is the state we compute.

**How it was measured.** The converged He/H = 0.0793 state of section 13.2 was
reloaded with `Load IC? True` (and `Secondary_ionization: Immediate`, so that
the coupling the converged state carries is on from step 0 -- with the STAGED
default the reloaded state is a different one, and n(H I) at the base moves by
a factor 4.6) into a scratch build of this tree carrying one output-only
addition: every term of the `fvec(4)` H2 balance row of
`System_HeH_mol.f90`, and the branch rates of H2+, H3+ and HeH+, written per
cell from the same arrays as `Ion_species.txt`. Below 1.5 R_p the reloaded
state reproduces the reference to 7e-4 in n(H2), 4e-5 in n(H I) and 8e-5 in
n(H II); `log10 Mdot` is 10.42 either way.

**Gross rates are not the answer, because the network cycles.** In a cell at
chemical equilibrium the production and destruction sums of the H2 row are
equal by construction (measured: they agree to 1 part in 1e6), and most of
both is a null cycle. R20 (He+ + H2 -> HeH+ + H) followed by R19
(HeH+ + H -> H2+ + He) and R9 (H2+ + H -> H+ + H2) returns the H2 and leaves
only a He+ + H -> He + H+ charge exchange behind. The budget below is
therefore **cycle-resolved**: each molecular ion is assigned its expected net
H2 yield from its own branch rates (H2+ from R5/R8/R9, H3+ from R6/R7/R11,
HeH+ from R16/R18/R19, solved as a 3-state linear system), and each primary
channel is charged the net H2 it actually removes. The check that the
bookkeeping is right is that the net loss then equals the external sources --
three-body association R15 plus the He(2^3S) + H associative branch -- to
1e-4 at every radius.

| | 1.05 R_p | 1.15 R_p | 1.25 R_p |
|---|---:|---:|---:|
| T [K] | 795 | 502 | 1220 |
| v [cm/s] | -3.9 | 3.87e2 | 1.09e4 |
| `f(H2)` | 0.966 | 1.43e-2 | 1.47e-4 |
| n(He+)/n(He) | 6.3e-11 | 8.9e-5 | 1.8e-2 |
| **net H2 loss, total [cm^-3 s^-1]** | **6.75e3** | **1.33e2** | **8.51e1** |
| net, H2 + hv (P_H2, 15.4 eV) | 5.85e3 (87%) | 5.30 (4.0%) | 1.16e-1 (0.1%) |
| net, R20 He+ + H2 -> HeH+ + H | 2.62e2 (3.9%) | 8.93e1 (67%) | 6.39e-1 (0.8%) |
| net, R17 He+ + H2 -> He + H + H+ | 2.82e2 (4.2%) | 3.70e1 (28%) | 8.43e1 (99%) |
| net, R13 H+ + H2 + M -> H3+ + M | 3.56e2 (5.3%) | 7.19e-1 (0.5%) | 9.1e-6 |
| net, R23 He+ + H2 -> H2+ + He | 4.39 (0.07%) | 1.04 (0.8%) | 9.0e-3 |
| net, R12 H2 + M -> H + H + M | 1.1e-11 | 3.4e-33 | 1.4e-11 |
| net, R14 H2 + e -> H + H + e | 5.9e-19 | 5.3e-38 | 2.5e-13 |
| net, R10 H+ + H2(v>=4) | 3.9e-3 | 6.1e-13 | 9.0e-5 |
| net, He(2^3S) + H2 ionization | 6.1e-5 | 7.9e-3 | 3.2e-3 |
| net, Lyman-Werner | 0 | 0 | 0 |
| source R15 (3-body) | 6.75e3 | 1.09e2 | 3.4e-3 |
| source He(2^3S) + H -> HeH+ route | -1e-6 | 2.44e1 | 8.50e1 |
| `tau_chem`(H2) = n(H2)/net [s] | 8.8e8 | 1.3e7 | 8.0e3 |
| `H_p/v` [s] | 1.8e7 | 2.4e5 | 2.5e4 |
| `r/v` [s] | 9.5e8 | 1.0e7 | 4.0e5 |
| **&#124;r^-2 d(r^2 v n_H2)/dr&#124; [cm^-3 s^-1]** | **2.96e5** | **1.55e4** | **4.74e1** |
| **that term / net chemical rate** | **44** | **117** | **0.6** |

**Four things the table settles.**

* **It is not thermal dissociation, and section 13.5 was right about that.**
  R12 runs at 1e-11 to 1e-33 cm^-3 s^-1 -- between fourteen and thirty-five
  decades below the net rate. R14, the electron-impact channel, is smaller
  still. Lyman-Werner is identically zero everywhere in the domain, as the run
  configuration requires.
* **Below ~1.10 R_p it is H2 photoionization**, 87% of the net at 1.05 R_p,
  and it removes more than one H2 per event: the H2+ it makes is consumed by
  R8 (H2+ + H2 -> H3+ + H), which eats a second H2, and the H3+ then returns
  an H2 only through R6, whose branching against R7 (3H) is 2.16/(2.16+5.04)
  = 0.30. The expected net yield per H2+ created is -0.66 at 1.05 R_p, i.e.
  each photoionization costs 1.66 H2. The He recombination photon coupling
  `dP_H2_hrc` of section 120 contributes nothing to P_H2 below 1.10 R_p and
  2-16% of it at 1.15-1.25 R_p (measured by an A/B on `He_rec_coupling`), so
  it is not a driver of the truncation.
* **Through the front itself, 1.10-1.27 R_p, it is helium ions.** At 1.15 R_p
  R20 and R17 carry 95% of the net between them; at 1.25 R_p R17 alone carries
  99%. Both are He+ + H2 reactions, and the He+ fraction rises from 6e-11 at
  1.05 R_p to 1.8e-2 at 1.25 R_p. R20's gross rate is ten times the whole net
  loss and only 6.7% of it survives the HeH+ -> H2+ -> H2 return, which is why
  the gross table names the wrong reaction.
* **Above ~1.2 R_p the surviving H2 is not left over from below, it is made in
  place.** The three-body source R15 has fallen to 3.4e-3 cm^-3 s^-1 at
  1.25 R_p while the whole net loss is 85, and what balances it is the 10%
  associative branch of He(2^3S) + H -> HeH+ + e (Garcia Munoz 2025, network
  rows 199 and 203, as cited at row (7) of `mol_heh_rows`), routed to H2
  through R19 and R9. The `f(H2)` = 1e-4 to 1e-6 floor
  of the outer wind is therefore a helium-metastable product and depends
  on that branching ratio.

**Why this is read as a closure problem and not a rate problem.** The last
two rows of the table evaluate, on the converged solution, the one term the
molecular network does not carry: `f(H2)` is an equilibrium root solved cell
by cell inside `ioniz_eq`, and H2 is transported only when `thereis_oxychem`
is on and the carrier operator of `diffusive_photochemistry.f90` runs -- which
requires metals, and this run has none. Evaluated on our own solution, that
omitted divergence is 44x the net chemical rate at 1.05 R_p, 117x at 1.15 R_p,
and only at 1.25 R_p does the chemistry finally win (0.6x). The equilibrium
profile is therefore not a steady state of the full continuity equation
anywhere below ~1.25 R_p; a self-consistent solution there must be shallower,
which is the direction the disagreement asks for. The measurement is one-sided
and should be read that way: it shows the closure is violated on its own
solution, not how much H2 a transported solution would keep.

Both comparison papers carry the term we omit. Koskinen et al. (2022)
Appendix B equation (B1) is a continuity equation for each species with an
advection and a diffusion velocity, and their Appendix B states of Figure 16
that "Chemical loss of H2 (dissociation and dissociative ionization
reactions) dominates over production but outflow brings fresh H2 up from
deeper layers of the atmosphere and primarily balances the total loss rate."
Frelikh & Murray-Clay (2026) plot "advection of H2" as a term of their
Figure 7 and note that its sign "represents an addition of H2 from the
boundary". Neither paper balances H2 chemically; we do, because that is all
our closure can do.

**The reaction networks compared.** Ours is transcribed from Koskinen's
Table 1, so against Koskinen the network is not a candidate at all. Against
Frelikh & Murray-Clay's Table 1 there are three differences that touch the
front, and one of them is large:

| channel | EXHALE (= Koskinen Table 1) | Frelikh & Murray-Clay Table 1 | ratio at 1220 K |
|---|---|---|---|
| He+ + H2 -> He + H + H+ | R17, `1e-9 exp(-5700/T)` (Moses & Bass 2000) | k16, `8.8e-14`, T-independent (Schauer et al. 1989) | ours 106x larger (657x at 2000 K, 0.13x at 502 K) |
| H3+ + H -> H2+ + H2 | R11, `2.1e-9 exp(-20000/T)` (Harada et al. 2010) | k11, `2.0e-9`, no barrier (Yelle 2004) | ours 1.2e7x smaller; this is an H2 SOURCE |
| H2 + hv -> H + H | opt-in Lyman-Werner (`lyman_werner.f90`), off here | k22, always on, from their section 2.10 cross sections | a sink they carry and we switched off |
| He+ + H2 -> HeH+ + H | R20, `4.2e-13` | k15, `4.2e-13` | identical |
| HeH+ + H2, HeH+ + H, HeH+ + e | R18/R19/R16 | k17/k18/k20 | identical |
| H3+ + e -> H2 + H : -> 3H | R6/R7, `2.16e-8`/`5.04e-8` (Larsson 2008) | k6/k7, `2.9e-8`/`8.6e-8` (Sundstrom 1994 / Datz 1995) | H2 return branch 0.30 vs 0.25 |
| H + H + M -> H2 + M | R15, `2.8e-31 T^-0.6` (Cohen & Westberg 1983) | k13, `8e-33 (300/T)^0.6` = `2.45e-31 T^-0.6` | ours 1.14x larger, documented at R15 |

The R17 line is the one worth carrying forward. It is our dominant net sink at
the top of the front (99% at 1.25 R_p), it is the same reaction and the same
primary measurement (Schauer et al. 1989) in both networks, and the two
networks differ only in how that measurement is apportioned -- Table 1 gives
the dissociative channel the Moses & Bass (2000) Jovian-ionosphere fit and the
HeH+ channel 4.2e-13 on top of it, while Frelikh (and Garcia Munoz 2025) give
the dissociative channel a T-independent value of order 1e-13 and nothing
else. The header of `mol_rates.f90` already records this disagreement and
already says "R17 is the dominant H2 sink and H+ source once He+ is present";
this measurement confirms that it is, at exactly the radii the gate is about.
**It cannot explain the disagreement with Koskinen**, who use the same R17 and
still keep H2 at all altitudes -- which is itself an argument that the
transport term, not the rate, is what separates us from them.

**Base depth moves the wrong way.** Section 13.5 recorded, as a caveat, that
our base is nine times deeper than Koskinen's 1 microbar and argued that this
should give us more H2, not less. That is now measured rather than argued. Two
12,000-step snapshots from a cold start, identical but for the base level
(`Log10 lower boundary number density` 14.00 -> 13.046, which is what actually
sets the level; the `p_base` key of `base.inp` only chooses where the
chemical-equilibrium fit is evaluated and is inert when `q_H2_base` is given):

| | base 8.99e-6 bar | base 1.00e-6 bar |
|---|---:|---:|
| H2 front, `f(H2) = 0.5` [R_p] | 1.154 | 1.074 |
| `f(H2)` down 4 decades by [R_p] | 1.230 | 1.152 |
| n(H3+) peak [cm^-3] (at [R_p]) | 6.8e4 (1.030) | 4.6e4 (1.000) |
| n(H3+) > 1 cm^-3 out to [R_p] | 1.190 | 1.100 |
| `log10 Mdot` [g/s] | 10.57 | 10.42 |
| T range [K] | 847-2828 | 697-2456 |

Sign and order of magnitude only, per item (P17): these are relaxation
snapshots at `du` = 9.9e-3 and 4.3e-3, not converged states, and the
8.99e-6 bar snapshot's front sits 0.05 R_p outside the converged reference's
1.101 R_p, which is the size of the snapshot error. **The direction is
unambiguous: a shallower base makes the truncation sharper and moves it
inward.** Raising our base to Koskinen's level would make the disagreement
larger, so base depth is not the explanation.

**The energy side, for the record.** Section 13.5 also recorded that we are
1700-2800 K colder than Koskinen. That is not a cooling-channel difference in
any obvious way: over 1.2-2.0 R_p the total radiative cooling is 21 to 156
times below the total heating, so the wind's temperature is set by expansion
and advection, not by a radiative channel.

| r [R_p] | heating (share) | cooling (share) | cool/heat |
|---|---|---|---:|
| 1.20 | He I 0.478, H I 0.415, He(2^3S) Penning 0.060, He recomb. 0.046 | recomb. 0.784, bremsstrahlung 0.216 | 6.4e-3 |
| 1.50 | H I 0.699, He I 0.209, He recomb. 0.048, He(2^3S) Penning 0.039 | recomb. 0.657, bremsstrahlung 0.260, He I lines 0.083 | 2.3e-2 |
| 2.00 | H I 0.765, He I 0.159, He recomb. 0.051, He(2^3S) Penning 0.019 | recomb. 0.596, bremsstrahlung 0.249, He I lines 0.156 | 4.7e-2 |

Ly-alpha cooling is below 0.1% of the total everywhere in this range -- at
under 3000 K it is not excited -- and every metal column is identically zero
(metals off). Whether the temperature difference and the H2 difference have
the same cause is not settled here.

**What this does not show.** It does not show that a transported H2 would
survive to the radii Koskinen and Frelikh report; that requires solving the
transported problem, and the advection estimate above is made on a profile
that the closure itself produced. It does not rank the R17 rate against the
missing transport term -- both act at the top of the front and only one of
them has been evaluated. And it says nothing about the He/H = 30 arm or about
the metals-on gates, where the electron density and hence every dissociative
recombination branch is different.

**Suspected code defects: none found in this measurement.** The H2 balance row
of `System_HeH_mol.f90` was read term by term against the network it claims to
implement and every Koskinen Table-1 channel that touches H2 is present with
the published coefficient; the one deliberate departure (R15 from Cohen &
Westberg 1983 rather than Ham et al. 1970, and R12 by detailed balance from
it) is argued in place at the function. The cycle bookkeeping closes on the
external sources to 1e-4, which is an independent check that no channel is
missing from the row. The one structural gap is the one this section is about,
and it is a closure choice rather than a coding error: H2 has no transport
operator unless the oxygen carriers are switched on, and no comment in
`System_HeH_mol.f90` or `ionization_equilibrium.f90` states the range of
validity of the local-equilibrium assumption for H2. That belongs in the code
at the H2 row.

The figure is `docs/lower_atmosphere_figs/fig_h2_budget_heh0793.png`: the
net H2 removal rate in each channel and `f(H2)` against radius (top), and the
same channels as time scales against `H_p/v` and `r/v` (bottom), with the
omitted advection term drawn on both. It is drawn by `fig2.py` in the campaign
scratch directory; it is not regenerated by
`docs/lower_atmosphere_figs/make_figures.py`.

### 13.7 The closure is now removable, and the transported problem has a moving front

**Judgment: the term section 13.6 measured to be missing is now available as an
option, `Molecular carrier transport: True`, on the molecular network alone
instead of only with the oxygen chemistry. It moves the H2 front and the H3+
extent in the direction the gate asks for. It does NOT yet give a steady state
on this planet: the transported front is still propagating outward at the gas
speed when the steady solver is handed it, so the option is default off and the
P23 gate is still open.** The full account, with every measurement, is section
128 of `docs/Update_EXHALE.md`; this section records what it means for the
ladder.

**What moves, on the marching path.** From the converged He/H = 0.0793 state of
section 13.2 (reloaded with `Secondary_ionization: Immediate`), 20,000 marching
steps with the transport on against the same 20,000 with it off:

| | `f(H2)=0.5` | `1e-2` | `1e-4` | `n(H3+)>1` to | `f(H2)` at 1.20 | at 1.25 |
|---|---:|---:|---:|---:|---:|---:|
| transport off | 1.1008 | 1.1542 | 1.2698 | 1.1810 | 1.01e-3 | 1.42e-4 |
| transport on | 1.1016 | 1.1602 | 1.3072 | 1.1976 | 3.26e-3 | 5.80e-4 |

and after a further 32,000 steps with the transport on the front is at 1.1731 /
1.3375 and `n(H3+) > 1` reaches 1.2157, i.e. **still moving outward and not
saturated**. The control is what makes this readable: the same case with the
transport OFF, run to 141,000 steps, moves the other way (`f(H2)` at 1.25 R_p
from 1.42e-4 to 1.25e-4, the four-decade point from 1.2698 to 1.2635 R_p). So
the outward motion is the transport, not marching drift.

**What the state after one steady pass looks like.** With the wind
JFNK-converged (`||R|| = 8.7e-4`) and the carriers relaxed once at that wind,
`f(H2)` at 1.20 R_p is 212 times its untransported value, `n(H3+)` at the same
radius goes from 0.16 to 157 cm^-3, and the cold dip at the front is gone
(T 768 -> 1677 K at 1.20 R_p) -- the direction section 13.5 needs, where we sit
1700-2800 K below Koskinen. Beyond about 1.35 R_p the two solutions rejoin
(ratio 1.8, then 0.84 at 1.5 R_p and 0.99 at 2.0 R_p), which is the consistency
check that matters: where the Damkohler number is large the transported and the
local solutions agree, and the transport changes the front rather than the
wind. That state is not converged and nothing from it is quoted as a rate.

**Why the gate is still open.** The front propagates 0.082 R_p per outer pass
of 0.1 flow times, about 2050 cm/s against a local gas speed of 3000 cm/s, and
it had not stopped after two passes. A real feedback pushes it: where H2
arrives, the He-ionizing continuum is absorbed and n(He II) falls by a factor
50-100 at 1.15-1.20 R_p -- exactly where section 13.6 measured He+ + H2 to
carry 95% of the whole net H2 loss -- so the sink that would stop the front is
removed by the front itself. Whether that ends in one steady front or in two
states is not settled. Until a marching run long enough to stop the front says
where it stops, and a steady solve from that state returns `info = 0`, the
comparison of section 13.5 stands as written and the transported numbers above
are a direction, not a result.

**A separate finding about this planet's base, and a correction to a premise
this campaign briefly held.** The converged He/H = 0.0793 state's base cell
carries `v = -250 cm/s` and a cell-centred mass flux 196 times the wind's with
the opposite sign, alternating and decaying over the next four cells. The first
reading of that was that the base face is an outflow, so escape could not
replenish H2 from below on this planet whatever the closure. **That reading was
wrong.** Integrating the scheme's own steady mass residual down from the smooth
wind recovers the face flux it actually transports, and at the base face it is
+0.36 to +0.98 times the wind flux, OUTWARD -- i.e. the base is supplying the
wind at about the right rate, in every converged configuration measured. The
cell-centred velocity there is a collocated odd-even mode, not infall: `rho`,
`p` and `T` are monotone through the same cells and the layer carries 288 cells
per scale height (`docs/p44_base_sawtooth.md`, `TO_BE_DONE.md` item (V)).

The correction has a consequence for this section: the carrier operator's base
composition condition is gated on the wind's mass flux, not on the base cell's
velocity, and it therefore DOES act on this planet. Over 20,000 marching steps
it moves the base cell from 0.9985 to 0.9574 of the H nuclei in H2, against the
0.9251 the handoff states -- 56% of the jump -- while leaving the front where it
was (1.1016 -> 1.1017 R_p). Koskinen's replenishment mechanism is therefore not
excluded here by the discretization; whether our wind supplies enough of it is
the question the front trajectory above is still answering.

## 14. The ladder at 1 microbar

**Judgment first. The ladder repeats at the base level the gate's own radius
implies, and it is easier and flatter than it was at 9 microbar. Five of the six
rungs converge from a COLD START with no restart at all -- the marching reaches
the hand-off, the JFNK enters by itself and returns `info = 0` on both gates of
section 133 -- where at 9 microbar the in-run hand-off failed on every rung that
reached it and all six needed the section-12.6 restart. `log10 Mdot` is 10.28 or
10.29 on every rung, and 10.28 to 10.30 counting both convergence cycles, across
a factor 378 in He/H -- a spread of 0.01 to 0.02 dex against 0.10 dex over the
same six rungs at 9 microbar (0.17 dex counting section 13.2's extra
`He_diffusion` rung). The structure keeps the direction section 13.2 measured --
more helium means a thinner molecular layer, a hotter wind, an ionization front
closer in and a slower outflow -- but every feature has moved inward with the
lighter base, and at the helium-rich end the molecular layer is gone before the
first cell: at He/H = 30 the base cell holds 4.6e-4 of its H nuclei in H2 where
the handoff states 0.985, and n(H3+) never exceeds 0.36 cm^-3 anywhere. The one
rung that did not converge on its own is the gate itself, He/H = 0.0793, and it
did not fail a gate -- it ran out of its 200,000-step budget at `du = 1.5e-2`,
never reaching the `du < 1e-2` hand-off. Its restart closed in 3.7 seconds. On
the P23 comparison nothing changes: the H2 front is now at 1.043 r_base instead
of 1.101, so the disagreement with Koskinen's H2-at-all-altitudes is larger, not
smaller.**

Measured 2026-09-03 with a copy of the tree's `EXHALE.x` as built at 07:41 KST
(md5 `b07ed10c0976303be82d9b6c5ae1e447`, copied to the campaign scratch at
08:02 KST), the binary that passes the block-F goldens 10/10. The runs are in
that campaign directory under `.../scratchpad/ladder2/runs/`, one directory per
rung, `OMP_NUM_THREADS=8`, all six launched at once. No metals, no
`opacity.inp`, no carrier transport; `Secondary_ionization` at its STAGED
default for the cold start.

### 14.1 The six rungs, and how each was converged

Each rung is the corresponding retargeted input of section 136 taken as it
stands: the He/H = 0.0793 rung is `backup/regression/mol_base_handoff` with its
`maxsteps` pin removed, and the other five are `armHeH_0p3`,
`arm_heh1_x2matched`, `armHeH_3`, `armHeH_10` and `armHeH_30`. All six carry
`p_base 1.000e-06` and no `Log10 lower boundary number density`, so the base
level comes from the handoff, and their `q_H2_base` values hold the H2 nucleus
fraction at `x2 = 0.985447826` while He/H alone varies. The startup report
confirms the level on every rung, for instance

```
 - Base level: p =  1.0000E-06 bar (from base.inp p_base) -> n0 =  1.1690E+13 cm^-3
```

with `n0` running 1.1690e13 (He/H = 0.0793) to 6.4561e12 cm^-3 (He/H = 30). As
in section 13.1 the gate rung differs from the five arms in one further respect
-- it does not set `He_diffusion` -- and that difference rides along with He/H
between the first and second rung of the table below.

Every rung was marched from a cold start with `PLM+WENO3`,
`du_th [PLM,WENO3]: 0.5 1.0e-3`, `Solver: Newton` (hand-off at the default
`du < 1e-2`) and a budget of 200,000 steps, and then left to the automatic
hand-off of sections 126 and 133.

**Five of the six converged that way, with no intervention.** The staged
secondary ionization fires at the first `du = 1e-2` crossing (step 42,878 to
58,552 depending on the rung), the wind re-relaxes under the full physics,
crosses again, and the JFNK enters and returns `info = 0` on both gates. One of
them (He/H = 0.3) took three solves inside the same run and two (He/H = 1 and 3)
took two, separated by the diffusion outer passes; the other two took one each.
**No solve in this campaign returned anything but `info = 0`, and the line
"no descent direction exists" is not printed anywhere** -- the failure mode
sections 126 and 138 were opened on does not occur on this ladder. The section
126 escape is not idle, though: on the He/H = 1 rung the damped Gauss-Newton
step is taken 8 times in the first solve, each time cutting the merit by a
factor 1.37 to 4.01 (`||Fs||2` 9.04e1 -> 2.51e1 at iteration 3,
2.51e1 -> 6.27e0 at iteration 4, and so on).

**The gate rung is the exception, and it is a budget exception rather than a
gate failure.** He/H = 0.0793 crossed `du = 1e-2` at step 79,208, the staged
switch fired, `du` rose to 1.03e-1 and then decayed monotonically but slowly,
reaching only 1.51e-2 at step 200,000 -- 1 h 24 m of wall time, and no second
crossing, so the JFNK was never armed. Restarted the section 12.6 way
(`Load IC? True`, `Secondary_ionization: Immediate`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-3`, `Solver: Newton 100.0`) it marched 2 steps,
entered the solver and met both gates at `||R|| = 7.742e-6` and flux spread
4.975e-3, in **3.7 seconds**.

### 14.2 The rungs

> **PRE-P53 (marked 2026-09-03).** These six states predate the
> composition-dependent caloric EOS and the molecular reaction heat, and they
> were finished on the residual measure that preceded P54. Section 14.8 has the
> same six rungs re-converged and re-finished, with a before/after column.

Six `info = 0` states, both gates met, `OMP_NUM_THREADS=8`. `log10 Mdot` is the
number the run prints (the `Rate/2 + Mdot/2` convention, as in section 13.2);
`||R||` is on the section-133 scale, where converged states sit at 1e-6 to
1e-5, so it is NOT comparable with the section 13.2 row, which is on the old
scale and larger by a factor 50 to 1000. All radii and profile values are read
on the physical cells only, the file's two ghost rows at each end dropped
(section 133.6).

> **This table is one grid, and the ladder has since been run on three
> (section 14.7).** All six rungs were re-converged on a 2x refined base grid
> (`Base grid [dr,cells]: 1.0e-4 100`, `Grid cells: 550`) and the gate rung and
> He/H = 0.3 also on a 4x grid, every state `info = 0` on both gates. Three
> things follow for the column values below.
>
> **(a) The He/H = 0.3 column carries a first-cell artifact, and the quantities
> that cell sets move when it is refined away**: `n(H3+)` at the peak
> **4.053e5 -> 7.535e4 cm^-3** (the peak moving 1.0002 -> 1.0022 R_p),
> `T minimum` **312 K -> 341 K** (moving from cell 1 to 1.0418 R_p),
> `log10 Mdot` 10.29 -> 10.26, the ionization front 2.0788 -> 2.0296 R_p and the
> `f(H2) = 0.5` front 1.0330 -> 1.0308.
>
> **(b) The artifact is not confined to that rung -- refining moves it to
> He/H = 10**, whose 2x state has the non-monotone first cells (977, 939, 958,
> 1016, 1057, 1075 K) and whose peak `n(H3+)` rises 874 -> 1556 cm^-3 and
> `log10 Mdot` 10.29 -> **10.34**, the largest rate move on the ladder and
> reproduced by a second restart cycle. The other four rungs move by 0.01 to
> 0.6 percent in `x2` at cell 1 and 0.00 to 0.04 dex in rate.
>
> **(c) None of this is a discretization error settling.** The base layer is
> time dependent: restarted under pure time marching, these states swing
> `rho v r^2` at 1.03 R_p by 16 to 50 percent of the wind value over 5000 steps
> (section 14.7.4). The rows below, and their refined counterparts, are phase
> samples of that oscillation, and on the gate rung the 2x to 4x change exceeds
> the 1x to 2x change on every quantity but `x2` at cell 1. **Read the values
> below as one grid's sample, not as the converged ladder**; section 14.7.2 has
> all six rungs on two grids side by side.

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 |
|---|---:|---:|---:|---:|---:|---:|
| n(JFNK), all info=0 | 1 | 3 | 2 | 2 | 1 | 1 |
| final `\|\|R\|\|` (section 133 scale) | 7.74e-6 | 6.78e-6 | 9.22e-7 | 1.90e-6 | 5.94e-6 | 8.64e-6 |
| flux spread, the gate's own | 4.98e-3 | 2.58e-3 | 1.78e-3 | 2.43e-3 | 4.27e-3 | 3.59e-3 |
| `du` of the solution | 1.51e-2 | 5.94e-3 | 5.61e-3 | 5.58e-3 | 5.33e-3 | 6.24e-3 |
| **log10 Mdot [g/s]** | **10.29** | **10.29** | **10.29** | **10.28** | **10.29** | **10.29** |
| marching steps | 200,000 + 2 | 60,552 | 56,226 | 50,426 | 55,335 | 44,878 |
| wall clock, 8 threads | 1 h 24 m + 3.7 s | 27 m 41 s | 25 m 42 s | 23 m 39 s | 27 m 02 s | 21 m 53 s |
| x2 at the ghost (the handoff) | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 |
| x2 at cell 1 | 0.9588 | 0.9883 | 0.8869 | 0.7453 | 0.3538 | **4.56e-4** |
| H2 front `f(H2) = 0.5` [R_p] | 1.0431 | 1.0330 | 1.0202 | 1.0099 | below cell 1 | below cell 1 |
| `f(H2) = 1e-2` [R_p] | 1.0799 | 1.0643 | 1.0397 | 1.0224 | 1.0080 | below cell 1 |
| `f(H2) = 1e-4` [R_p] | 1.1690 | 1.1138 | 1.0882 | 1.0637 | 1.0518 | 1.0465 |
| T at the H2 front [K] | 562 | 681 | 921 | 915 | -- | -- |
| H3+ peak [R_p] | 1.0044 | 1.0002 | 1.0002 | 1.0002 | 1.0002 | 1.0002 |
| n(H3+) at the peak [cm^-3] | 1.76e5 | 4.05e5 | 2.26e4 | 5.97e3 | 8.74e2 | 3.64e-1 |
| n(H3+) > 1 cm^-3 out to [R_p] | 1.0997 | 1.0810 | 1.0508 | 1.0475 | 1.0122 | never > 1 |
| ionization front `x_HII = 0.5` [R_p] | 2.376 | 2.079 | 1.757 | 1.531 | 1.405 | 1.354 |
| `x_HII` at the outer boundary | 0.816 | 0.858 | 0.897 | 0.919 | 0.927 | 0.930 |
| n_e at cell 1 [cm^-3] | 7.63e5 | 4.71e5 | 5.52e6 | 1.83e7 | 4.31e7 | 5.33e7 |
| n_e maximum [cm^-3] (at [R_p]) | 8.42e7 (1.093) | 9.49e7 (1.200) | 1.23e8 (1.143) | 1.54e8 (1.081) | 1.78e8 (1.064) | 1.73e8 (1.052) |
| n_e at the outer boundary [cm^-3] | 5.34e6 | 4.42e6 | 3.52e6 | 2.97e6 | 2.83e6 | 2.76e6 |
| T minimum [K] (at [R_p]) | 424 (1.066) | **312 (1.0002)** | 629 (1.057) | 636 (1.042) | 747 (1.032) | 946 (1.025) |
| T maximum [K] (at [R_p]) | 2810 (1.646) | 3685 (1.500) | 5064 (1.431) | 6344 (1.396) | 7070 (1.426) | 7231 (1.453) |
| outflow speed at the boundary [km/s] | 8.84 | 8.08 | 7.28 | 6.80 | 6.62 | 6.57 |
| energy-floor hits | 0 | 0 | 0 | 0 | 0 | 0 |
| positivity limiter: face states | 0 | 0 | 0 | 0 | 0 | 0 |
| positivity: ghost cells / face fluxes | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| non-root acceptances (all ledgers) | 52 | 54 | 64 | 111 | 756 | 184 |
| projected/handback roots | 0 | 0 | 0 | 0 | 0 | 0 |
| cells with no in-simplex start | 0 | 0 | 0 | 0 | 0 | 0 |
| largest element-budget violation | none reported | 9.7e-2 | 1.6e-1 | 2.4e-1 | 4.2e-3 | 1.6e-1 |
| trial states refused (section 138.3) | 0 | 0 | 0 | 0 | 0 | 0 |
| base inflow, max Mach (cold start) | 3.54e-2 | 4.04e-2 | 4.11e-2 | 5.95e-2 | 7.28e-2 | 7.83e-2 |

**The 9 microbar ladder beside it.** Section 13.2's states, on the same rows
that mean the same thing. The `f(H2)` row there is "down 4 decades from the
base", not the absolute `1e-4` of the table above, so the two are only
comparable where the base fraction is near 1; the "4 decades below cell 1"
measure on the new states is 1.172, 1.115, 1.090, 1.067, 1.056 and 1.216 R_p,
the last of which is meaningless because that rung's cell 1 already holds
4.6e-4. The section 13.2 ionization fronts for He/H >= 3 are the asterisked
ones, superseded by the section 125 re-convergence (1.618 at He/H = 3, 1.758 at
He/H = 1, 1.444 at He/H = 30).

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 |
|---|---:|---:|---:|---:|---:|---:|
| log10 Mdot, 9 microbar (section 13.2) | 10.42 | 10.35 | 10.32 | 10.36 | 10.42 | 10.39 |
| log10 Mdot, 1 microbar | 10.29 | 10.29 | 10.29 | 10.28 | 10.29 | 10.29 |
| H2 front `f = 0.5`, 9 microbar | 1.1012 | 1.0681 | 1.0469 | 1.0467 | 1.0357 | 1.0250 |
| H2 front `f = 0.5`, 1 microbar | 1.0431 | 1.0330 | 1.0202 | 1.0099 | below cell 1 | below cell 1 |
| n(H3+) peak, 9 microbar | 3.09e5 | 3.20e5 | 3.16e5 | 2.42e5 | 1.95e5 | 1.48e5 |
| n(H3+) peak, 1 microbar | 1.76e5 | 4.05e5 | 2.26e4 | 5.97e3 | 8.74e2 | 3.64e-1 |
| n(H3+) > 1 out to, 9 microbar | 1.181 | 1.118 | 1.081 | 1.058 | 1.049 | 1.035 |
| n(H3+) > 1 out to, 1 microbar | 1.100 | 1.081 | 1.051 | 1.048 | 1.012 | never |
| ionization front, 9 microbar | 2.573 | 2.040 | 1.593 | 1.396\* | 1.333\* | 1.241\* |
| ionization front, 1 microbar | 2.376 | 2.079 | 1.757 | 1.531 | 1.405 | 1.354 |
| T maximum, 9 microbar | 2751 | 3642 | 4967 | 5864 | 6328 | 6660 |
| T maximum, 1 microbar | 2810 | 3685 | 5064 | 6344 | 7070 | 7231 |
| v at the boundary, 9 microbar | 8.97 | 8.28 | 7.34 | 6.85 | 6.74 | 6.61 |
| v at the boundary, 1 microbar | 8.84 | 8.08 | 7.28 | 6.80 | 6.62 | 6.57 |

Read across, **three things change with the base level and one does not.**

* **The rate becomes flat.** 10.28 or 10.29 on every rung (10.28 to 10.30
  counting both cycles), against 10.32 to 10.42 over the same six rungs at
  9 microbar. Part of the old spread was the `He_diffusion` switch
  and part was the measure of section 127; what is left here is smaller than the
  printed precision. Whatever sets the rate on this planet, it is not the helium
  abundance, and at 1 microbar it is not the base composition either.
* **Everything molecular moves inward and thins.** The H2 front moves in by
  0.03 to 0.06 R_p on every rung that still has one, and n(H3+) at the peak
  falls by a factor 1.8 at He/H = 0.0793, 14 at He/H = 1 and 4e5 at He/H = 30.
  This is the section 136.3 direction, now measured on converged states rather
  than on 12,000-step snapshots.
* **The wind gets hotter.** T maximum is 1 to 12% higher on every rung, the
  increase growing with He/H (1.2% at He/H = 0.3, 11.7% at He/H = 10). The
  ionization front moves OUT on five rungs, by 0.04 to 0.16 R_p, and IN by
  0.20 R_p at He/H = 0.0793 -- but the He/H >= 3 entries of section 13.2 are the
  asterisked ones. Against the section 125 corrections instead, the He/H = 1
  front is where it was (1.758 -> 1.757) and the He/H = 3 and 30 fronts move IN
  by 0.09 R_p each.
* **The monotone ordering in He/H survives unchanged.** More helium still means
  a thinner molecular layer, a hotter wind, an ionization front closer to the
  planet and a slower outflow, and every one of those rows is monotone.

**Ledger.** No rung activated the energy floor, no rung fired the positivity
limiter of section 138.2 in any form, and no trial state was refused for an
uncertified composition. Root acceptance is clean in the sense that matters --
zero projected/handback roots and zero cells with no in-simplex starting point
on all six -- but the relaxation amnesty is used, 52 to 756 non-root
acceptances over the whole marching history of each run, with He/H = 10 the
worst. The largest element-budget violation is 1e-1 to 2e-1 on four of the six
rungs, larger than the 1e-3 to 1e-2 of the 9 microbar ladder's low rungs; that
is a diagnostic maximum over the whole run, not a property of the final state,
and it has not been traced to a radius here. **Every rung admitted a subsonic
base inflow during its cold start**, max Mach 3.5e-2 to 7.8e-2 rising
monotonically with He/H, where none of the six 9 microbar rungs did.

### 14.3 Reproducibility: a second chain cycle on every rung

Each converged state was restarted from its own output the same way
(`Load IC? True`, `Secondary_ionization: Immediate`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-3`, `Solver: Newton 100.0`) and re-converged.
**All six closed on the first restart, `info = 0` on both gates.**

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 |
|---|---|---|---|---|---|---|
| log10 Mdot, cycle 1 / 2 | 10.29 / 10.30 | 10.29 / 10.29 | 10.29 / 10.29 | 10.28 / 10.28 | 10.29 / 10.29 | 10.29 / 10.30 |
| `f(H2) = 0.5` [R_p] | 1.0431 / 1.0423 | 1.03302 / 1.03296 | 1.02025 / 1.02025 | 1.00995 / 1.01009 | -- | -- |
| `f(H2) = 1e-4` [R_p] | 1.1690 / 1.1614 | 1.11380 / 1.11348 | 1.08824 / 1.08824 | 1.0637 / 1.0617 | 1.05183 / 1.05156 | 1.0465 / 1.0428 |
| n(H3+) peak [cm^-3] | 1.760e5 / 1.767e5 | 4.053e5 / 4.037e5 | 2.2581e4 / 2.2572e4 | 5966.9 / 5969.5 | 874.41 / 874.21 | 0.36395 / 0.36393 |
| ionization front [R_p] | 2.376 / 2.393 | 2.0788 / 2.0788 | 1.75737 / 1.75737 | 1.531 / 1.538 | 1.405 / 1.416 | 1.354 / 1.372 |
| T maximum [K] | 2810 / 2814 | 3685 / 3684 | 5064 / 5063 | 6344 / 6332 | 7070 / 6994 | 7231 / 7104 |

**The reproducibility is not one number, and it splits by rung and by
quantity.** The He/H = 1 rung reproduces its H2 front to six significant digits
(1.02025 twice, and 1.08824 twice at `f = 1e-4`); He/H = 0.3 reproduces it to
five. From He/H = 3 upward the molecular layer still reproduces to three or four
digits but **the wind does not**: the ionization front moves by 0.43%
(He/H = 3), 0.72% (10) and 1.32% (30), and T maximum by 0.20%, 1.08% and 1.76%
on the same three. The gate rung sits between them, its `f = 0.5` front
reproducing to three digits (1.0431 against 1.0423) and its ionization front to
two (2.376 against 2.393). `log10 Mdot` reproduces to the printed digit
everywhere except the two rungs that read 10.29 and 10.30, i.e. one in the last
digit. The ionization-front entries that read identically (He/H = 0.3 and 1) are
exactly equal because the measure is the radius of the first cell with
`x_HII > 0.5` and both cycles land on the same cell -- that is agreement to a
cell width (9.2e-3 R_p at He/H = 1, 1.3e-2 at He/H = 0.3), not to six digits.

The measure this is against is section 10 of `docs/p44_base_sawtooth.md`, which
reproduced `r_front` to five digits at 1 microbar on the He/H = 0.0793
configuration of that section (1.03551 twice). That is met on the He/H = 1 and
0.3 rungs and not on the gate or the helium-rich ones. On the two most
helium-rich rungs the second cycle's flux spread is 4 to 5 times smaller than
the first's (4.27e-3 -> 8.50e-4 at He/H = 10, 3.59e-3 -> 8.55e-4 at 30), so the
two states are both inside the gate but not at the same distance from it, and
the wind differences track that.

### 14.4 What the base does on this ladder

Two features of the base are worth recording because both are new at 1 microbar
and both are reproducible across the chain cycle.

**(a) The helium-rich rungs lose their molecular layer inside the first cell.**
The handoff states `x2 = 0.985448` at the ghost on every rung, by construction.
What cell 1 carries after the run is 0.959 at He/H = 0.0793 and 0.988 at
He/H = 0.3, but 0.354 at He/H = 10 and **4.56e-4 at He/H = 30** -- three decades
below the ghost, across 2e-4 R_p. n(H3+) at He/H = 30 never exceeds
0.36 cm^-3 anywhere in the domain, so that rung has no molecular layer to
describe. The likely reading is that this is the interior chemistry rather than
the boundary: n_e at cell 1 rises monotonically along the ladder, 7.6e5 to
5.3e7 cm^-3, because each H nucleus is accompanied by up to 30 helium atoms
whose photoionization supplies the electrons, and every H3+ destruction branch
is proportional to n_e. That reading is not tested here -- no A/B was run on the
helium photoionization -- and the alternative, that the first cell cannot
resolve the composition gradient the handoff sets up, is not excluded.

**(b) The He/H = 0.3 rung carries a cold cell 1 and an alternating temperature
over the first six cells.** Its base cells read T = 312, 1296, 1427, 1433, 1144,
1034 K against a ghost at 1140 K, with `v` alternating +195, -28, +41, -18 cm/s,
and n(H2) at cell 1 is 3.6 times the ghost value. The second chain cycle
reproduces it (T at cell 1 = 316 K). **No other rung does this**: the five
others have monotone T through the same cells (for instance 1098, 1073, 1049,
1026, 1005, 985 K at He/H = 0.0793), with the same single-cell velocity
artifact `v(1)` = -217 to -354 cm/s that section 7.2 of
`docs/p44_base_sawtooth.md` identifies as a first-order-in-`dr` discretization
error. The He/H = 0.3 rung is therefore the only one on which that artifact has
reached the temperature, and its `T minimum` entry in the table above is that
cell, not a property of the molecular layer. Read the He/H = 0.3 column with
that in mind; the layer's own minimum on that rung, outside the first six cells,
is 471 K at 1.0590 R_p, and the second cycle puts it at 471 K and 1.0590.

> **Paragraph (b) is a property of this grid, and section 14.7 measures what
> happens on two more.** The cold cell 1 and the alternation are GONE on a
> 2x refined base grid (`Base grid [dr,cells]: 1.0e-4 100`, `Grid cells: 550`,
> cold start, four JFNK solves all `info = 0`, flux spread 4.667e-3): T over the
> first six cells is monotone -- 1104, 1095, 1086, 1078, 1070, 1061 K -- `v` is
> -125, +19.7, +6.2, +9.6, +8.1, +8.6 cm/s, and n(H2) at cell 1 is 0.96 of the
> ghost value rather than 3.6 times it. **The artifact was not cosmetic**: with
> it the rung's peak `n(H3+)` falls 4.053e5 -> 7.535e4 cm^-3, the layer minimum
> quoted just above moves 471 K at 1.0590 -> 341 K at 1.0418 R_p, and
> `log10 Mdot` moves 10.29 -> 10.26.
>
> **But the claim that it is "the only one" does not survive the refinement.**
> On the 2x grid the He/H = 10 rung acquires the same alternation (977, 939,
> 958, 1016, 1057, 1075 K over its first six cells, reproduced by its second
> restart cycle), so what the first cell does is a property of the grid and the
> composition together and not of this rung. `docs/p23_thermal_budget.md`
> section 8.2 has the He/H = 0.3 side-by-side table; section 14.7.2 has all six
> rungs on both grids, and section 14.7.4 shows that the layer these cells sit
> in is oscillating rather than steady on either grid.
>
> Paragraph (a) is NOT stale: the same refinement leaves the He/H = 30 first
> cell where it was (`x2` 4.563e-4 -> 4.593e-4, the He II jump across the ghost
> 842x -> 833x), so that drop is the interior chemistry and not the resolution.
> `docs/p23_thermal_budget.md` section 8.1 also identifies the reaction that
> carries it -- He+ + H2 -> H+ + H + He (R17), 93.5 percent of the H2
> destruction at cell 1, with R20 a further 6.4 percent -- which settles the
> "not tested here" caveat above in favor of the interior-chemistry reading.

### 14.5 The published profiles at 1 microbar

> **PRE-P53 (marked 2026-09-03).** Our column below is a state converged before
> the composition-dependent caloric EOS and the molecular reaction heat
> existed. Section 14.8.4 has the temperature comparison re-measured; the
> conclusion of this section is unchanged by it, and the published columns are
> the same digitization.

The material of `docs/p23_published_profiles.md` section 6, with our column
replaced by the converged 1 microbar He/H = 0.0793 state. The published columns
are unchanged -- they are the same digitization, re-read from the same arrays --
and are reproduced here so the three sit side by side. `f(H2)` is the fraction
of H nuclei bound in H2, `q(H2)` the volume mixing ratio against all heavy
particles. Entries at `r/r_base = 1.00` come from the near-vertical base columns
of the published figures and carry the caveat of that document's section 4.

| r/r_base | 1.00 | 1.05 | 1.10 | 1.15 | 1.20 | 1.30 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|---|---|---|
| **f(H2)**, ours at 1 microbar | 0.959 | 0.336 | 3.08e-3 | 1.59e-4 | 6.43e-5 | 3.52e-5 | 1.62e-5 | 4.28e-6 | 1.27e-6 |
| **f(H2)**, ours at 9 microbar (section 13.6) | 0.999 | 0.966 | 0.516 | 1.42e-2 | 1.10e-3 | 7.2e-5 | 2.8e-5 | 6.4e-6 | 1.5e-6 |
| **f(H2)**, Koskinen Fig. 8 | 0.984 | 0.982 | 0.976 | 0.963 | 0.936 | 0.843 | 0.688 | 0.508 | 0.403 |
| **f(H2)**, Frelikh Fig. 21 | 0.976 | 0.606 | 0.438 | 0.269 | 0.200 | 0.096 | 0.048 | 0.022 | 0.013 |
| **q(H2)**, ours at 1 microbar | 0.799 | 0.185 | 1.43e-3 | 7.36e-5 | 2.98e-5 | 1.63e-5 | 7.49e-6 | 1.98e-6 | 5.89e-7 |
| **q(H2)**, Koskinen Fig. 8 | 0.830 | 0.822 | 0.825 | 0.801 | 0.758 | 0.642 | 0.467 | 0.308 | 0.230 |
| **q(H2)**, Frelikh Fig. 21 | 0.817 | 0.387 | 0.254 | 0.142 | 0.102 | 0.046 | 0.023 | 0.010 | 6.1e-3 |
| **n(H3+)** [cm^-3], ours at 1 microbar | 1.45e5 | 743 | 1.11 | 1.88e-2 | 6.02e-3 | 1.69e-3 | 1.94e-4 | 6.74e-6 | 4.35e-7 |
| **n(H3+)**, ours at 9 microbar | 2.9e5 | 2.0e5 | 1.1e3 | 13.9 | 0.18 | 6.4e-3 | 7.9e-4 | 1.7e-5 | 6.9e-7 |
| **n(H3+)**, Koskinen Fig. 8 | 1.9e4 | 3.1e3 | 5.0e3 | 7.0e3 | 9.1e3 | 1.2e4 | 1.1e4 | 7.2e3 | 1.9e3 |
| **n(H3+)**, Frelikh Fig. 21 | 1.5e5 | 7.7e3 | 1.5e4 | 8.8e3 | 4.5e3 | 1.3e3 | 2.4e2 | 14.2 | 0.99 |
| **n_e** [cm^-3], ours at 1 microbar | 7.63e5 | 4.99e7 | 8.39e7 | 7.58e7 | 7.35e7 | 7.22e7 | 6.46e7 | 3.89e7 | 1.56e7 |
| **n_e**, Koskinen Fig. 8 (`em`) | 1.6e5 | 5.0e6 | 4.4e6 | 7.2e6 | 9.6e6 | 1.0e7 | 5.7e6 | 1.8e6 | 7.8e5 |
| **n_e**, Frelikh Fig. 21 (ion sum) | 1.6e7 | 2.8e8 | 2.0e8 | 1.9e8 | 1.6e8 | 1.2e8 | 6.6e7 | 2.8e7 | 1.1e7 |
| **T** [K], ours at 1 microbar | 1098 | 494 | 626 | 1167 | 1644 | 2273 | 2753 | 2671 | 2126 |
| **T**, ours at 9 microbar | 1100 | 795 | 775 | 502 | 759 | 1610 | 2470 | 2720 | 2250 |
| **T**, Koskinen Fig. 7 | 1082 | 1180 | 1566 | 1717 | 1901 | 2267 | 2871 | 4120 | 4742 |
| **T**, Frelikh Fig. 20(d) | 939 | 1137 | 2095 | 2470 | 2774 | 3253 | 3790 | 4151 | 4091 |

Landmarks, in the same normalization:

| | ours, 1 microbar | ours, 9 microbar | Koskinen Model A | Frelikh super-Earth |
|---|---|---|---|---|
| `f(H2)` down by 1 decade at | 1.062 r_base | 1.132 | not reached in the domain (7.16) | 1.296 |
| `f(H2)` down by 4 decades at | 1.172 | 1.273 | not reached | not reached (4.13) |
| `f(H2)` at `r/r_base` = 3 | 1.27e-6 | 1.5e-6 | 0.403 | 0.013 |
| n(H3+) maximum | 1.76e5 cm^-3 at 1.004 | 3.09e5 at 1.023 | 1.2e4, flat over 1.23-1.54 | 1.5e4 at 1.10 |
| n(H3+) at 2.00 r_base | 6.7e-6 | 1.7e-5 | 7.2e3 | 14.2 |
| T maximum, and where | 2810 K at 1.646 | 2751 K at 1.815 | 5527 K at 7.16 (still rising) | 4170 K at 2.10 |
| n_e maximum | 8.4e7 cm^-3 at 1.093 | 7.8e7 at 1.19 | 1.1e7 near 1.25 | 3.2e8 near 1.07 |
| outflow speed at the boundary | 8.84 km/s at 4.73 | 8.97 at 4.73 | 9.77 at 7.16 | 5.8 near 3.6 (their Figure 20 marker; the caption says 4.4 -- `p23_published_profiles.md` section 4) |

This section states the material and draws no verdict from it; the P23 gate is
section 13.5's and the transported-closure question is section 13.7's, and
neither is re-decided here. Two things about the material itself are worth
recording. First, **the direction section 13.6 measured is confirmed on
converged states**: the lighter base gives less H2 at altitude, not more, so the
base-depth caveat of section 13.5 is settled in the sense that raising our base
to Koskinen's own level widens the disagreement instead of closing it. Second,
**the temperature comparison changes shape**. At 9 microbar our profile had a
cold dip at the front (502 K at 1.15 r_base) and lay 401 to 2492 K below
Koskinen over 1.15-3.00 r_base; at 1 microbar the dip has moved inward to 1.05 r_base
(494 K), and through the front we have closed most of the gap: at 1.15, 1.20
and 1.50 r_base we are 550, 257 and 118 K below Koskinen where at 9 microbar we
were 1215, 1142 and 401 K below, and at 1.30 r_base the two profiles touch
(2273 K against 2267 K). Beyond the front they separate again, and by more than
before: at 3.0 r_base the gap is 2616 K against 2492 K at 9 microbar.

The figure is `docs/lower_atmosphere_figs/fig_ladder_1ubar.png`: `f(H2)`,
`n(H3+)`, `n_e` and `T` against `r/r_base`, one curve per rung of this ladder,
with the two digitized published profiles as points. It is drawn by
`make_fig.py` in the campaign scratch directory beside the runs; it is not
regenerated by `docs/lower_atmosphere_figs/make_figures.py`.

### 14.6 What was not measured

* **The carrier transport is off on every rung**, so section 13.7's question is
  untouched: this is the local-equilibrium closure re-measured at a new base
  level, not a transported ladder.
* **Nothing here re-runs the He/H = 30 base**, so whether its missing molecular
  layer is the interior chemistry or a resolution limit at the first cell is
  stated as a reading and not as a measurement. The A/B that would decide it --
  the same rung on a refined base grid, and the same rung with the helium
  photoionization suppressed -- was not run.
* **The element-budget violations of 1e-1 to 2e-1** on four rungs are reported
  as the run prints them and were not traced to a radius or a species. They are
  larger than the 9 microbar ladder's and that difference is not explained here.
* **The 200,000-step budget of the gate rung was not extended.** Whether it
  would have crossed `du = 1e-2` on its own at, say, 260,000 steps is not
  known; what is known is that its restart converged in 3.7 seconds from the
  state the budget ended on.
* **`log10 Mdot` is quoted to the two decimals the run prints.** The ladder's
  whole spread is 0.02 dex, which is at that precision, so the statement that
  the rate is flat is a statement that the differences are below what is
  printed, not a measurement of how much smaller they are.

### 14.7 Grid convergence of the ladder

**Judgment first, and it is not the one this section was set up to reach. The
six rungs were re-converged on a base grid refined by a factor 2, and the gate
rung and He/H = 0.3 also on a factor 4; all fourteen states return `info = 0` on
BOTH gates of section 133, and every one was reproduced by a second restart
cycle. But the quantities that move with the grid do not settle as a
discretization error settles. Only two do: the cell-1 velocity artifact, which
is first order in `dr` on every rung (ratios 1.80 to 2.84 from 1x to 2x, 2.19
and 3.27 from 2x to 4x), and the base composition `x2` at cell 1, which
reproduces to four digits from the 2x grid on. Everything else -- the rate, the
H2 fronts, the ionization front, the layer temperature minimum, and the
base-layer energy closure -- moves by more from 2x to 4x than it did from 1x to
2x on the gate rung, and changes sign there. The reason is measured, not
inferred: restarting each `info = 0` state under pure time marching and sampling
the mass flux at r = 1.03 every 500 steps shows that ratio swinging by 0.16 to
0.50 of the wind value over 5000 steps on all four states tested.** The base
layer is not steady, so the closure ratio, the sub-1.2 flux spread and the
front positions read off a converged state are phase samples of an oscillation.
Refining the grid does not remove that oscillation and does not systematically
reduce it: its swing grows from 0.16 to 0.20 at He/H = 0.3 and falls from 0.50
to 0.34 at He/H = 30. **The right reading of the tables below is therefore how
unsteady each grid's layer is, not how large its discretization error is**;
this is item (AA) of `TO_BE_DONE.md` and section 9.2 of
`docs/p23_thermal_budget.md`, and the P54 measurement of the same oscillation
(damped acoustic mode, period about 2100 steps, base to 1.10 R_p) is the
mechanism.

Measured 2026-09-03 with a copy of the tree's `EXHALE.x` taken at 10:58 KST
(md5 `9540ff5e9ee38422fc4ddca9febb5973`, tree file written 10:31 KST). That
binary is byte-identical to the section 14 binary on this configuration: an
800-step cold start of the He/H = 1 rung gives the same `Hydro_ioniz.txt` and
`Ion_species.txt` under both. The tree was rebuilt again later in the day; every
number here is from the 10:58 copy. Runs are in the campaign directory
`.../scratchpad/ladder2x/`, `OMP_NUM_THREADS=8`, all eight cold starts launched
at once.

**The three grids.** `Base grid [dr,cells]` sets the uniform base region at
fixed 0.01 R_p extent and `Grid cells` was raised in step, so the stretched
region keeps its 450 cells and the outer grid is unchanged -- the same ladder
`docs/p44_base_sawtooth.md` section 7.2 uses.

| label | keys | realized `dr` [R_p] | cells |
|---|---|---|---|
| 1x | (defaults) | 1.930e-4 | 500 |
| 2x | `Base grid [dr,cells]: 1.0e-4 100` + `Grid cells: 550` | 9.601e-5 | 550 |
| 4x | `Base grid [dr,cells]: 5.0e-5 200` + `Grid cells: 650` | 4.776e-5 | 650 |

#### 14.7.1 How each state was converged

Every rung was cold-started with the section 14.1 recipe (`PLM+WENO3`,
`du_th [PLM,WENO3]: 0.5 1.0e-3`, `Solver: Newton`, `Secondary_ionization` at its
STAGED default) and left to the automatic hand-off. **All eight converged that
way with no restart, the gate rung included** -- where at 1x that rung exhausted
its 200,000-step budget at `du = 1.5e-2` and needed the section 12.6 restart.
Its 2x cold start reached `info = 0` at step 161,323 in 88 m 51 s and its 4x
cold start at step 321,990 in 3 h 00 m.

`dt` is set by the base cells, and the ladder reproduces that to 1 percent on
every rung: the PLM to WENO3 switch, which is a physical state rather than a
step count, moves from 14,403-16,800 at 1x to 28,993-33,786 at 2x (ratios 2.005
to 2.013) and to 65,928 and 67,936 at 4x (2.009 and 2.011 against their own 2x).
The staged secondary-ionization switch gives the same ratios, 1.999 to 2.013.

**Cost, measured on this ladder.** Per step, 27.4-29.3 ms at 1x against
32.7-35.9 ms at 2x, a median +19 percent (the +12 percent of
`docs/p44_base_sawtooth.md` section 9.5 was measured with an alternating
protocol on an idle machine; this campaign shared 72 cores with other work).
Per unit physical time, exactly 2.0 times the steps. Wall clock to convergence
is 2.24 to 2.43 times the 1x value, median 2.40, against the 2.25 that section
predicts. **On the gate rung the refinement is cheaper than not refining**: 88 m
to an `info = 0` state at 2x against 84 m that reached no solution at 1x.

#### 14.7.2 The six rungs on two grids

Section 14.2's states are the 1x column (the gate rung's is its restart, the
only 1x state of that rung that converged). All fourteen states here returned
`info = 0` on both gates; `||R||` runs 1.9e-6 to 9.6e-6 and the gate's flux
spread 1.2e-3 to 4.9e-3.

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 |
|---|---:|---:|---:|---:|---:|---:|
| **log10 Mdot**, 1x | 10.29 | 10.29 | 10.29 | 10.28 | 10.29 | 10.29 |
| **log10 Mdot**, 2x | 10.25 | 10.26 | 10.27 | 10.28 | **10.34** | 10.30 |
| `f(H2) = 0.5` [R_p], 1x | 1.0431 | 1.0330 | 1.0202 | 1.0099 | below cell 1 | below cell 1 |
| `f(H2) = 0.5` [R_p], 2x | 1.0375 | 1.0308 | 1.0210 | 1.0091 | below cell 1 | below cell 1 |
| `f(H2) = 1e-4` [R_p], 1x | 1.1690 | 1.1138 | 1.0882 | 1.0637 | 1.0518 | 1.0465 |
| `f(H2) = 1e-4` [R_p], 2x | 1.1390 | 1.0922 | 1.0713 | 1.0684 | 1.0088 | 1.0093 |
| n(H3+) peak [cm^-3], 1x | 1.760e5 | 4.053e5 | 2.258e4 | 5.967e3 | 8.744e2 | 3.639e-1 |
| n(H3+) peak [cm^-3], 2x | 1.803e5 | **7.535e4** | 2.263e4 | 5.997e3 | **1.556e3** | 3.681e-1 |
| H3+ peak at [R_p], 1x | 1.0044 | 1.0002 | 1.0002 | 1.0002 | 1.0002 | 1.0002 |
| H3+ peak at [R_p], 2x | 1.0051 | 1.0022 | 1.0001 | 1.0001 | 1.0002 | 1.0001 |
| T minimum [K], 1x | 424 | **312** | 629 | 636 | 747 | 946 |
| T minimum [K], 2x | 318 | 341 | 470 | 639 | 939 | 1033 |
| T maximum [K], 1x | 2810 | 3685 | 5064 | 6344 | 7070 | 7231 |
| T maximum [K], 2x | 2847 | 3728 | 5110 | 6334 | 6794 | 7164 |
| ionization front [R_p], 1x | 2.3762 | 2.0788 | 1.7574 | 1.5312 | 1.4054 | 1.3541 |
| ionization front [R_p], 2x | 2.3047 | 2.0296 | 1.7269 | 1.5280 | **1.4725** | 1.3680 |
| `x2` at cell 1, 1x | 0.95876 | 0.98833 | 0.88685 | 0.74534 | 0.35378 | 4.5634e-4 |
| `x2` at cell 1, 2x | 0.95845 | **0.94034** | 0.88674 | 0.74603 | **0.43470** | 4.5925e-4 |
| T monotone over cells 1-6, 1x | yes | **NO** | yes | yes | yes | yes |
| T monotone over cells 1-6, 2x | yes | yes | yes | yes | **NO** | yes |
| `v(1)` [cm/s], 1x | -220.0 | -353.9 | -294.0 | -313.1 | -332.8 | -317.0 |
| `v(1)` [cm/s], 2x | -103.9 | -124.6 | -141.1 | -174.0 | -155.2 | -153.2 |
| `v(1)` ratio 1x/2x | 2.12 | 2.84 | 2.08 | 1.80 | 2.14 | 2.07 |
| closure at 1.03, 1x | 0.52 | 0.39 | 0.49 | 0.79 | 1.18 | 1.68 |
| closure at 1.03, 2x | 1.00 | 0.97 | 0.83 | 0.65 | **0.37** | 1.36 |
| flux spread `r >= 1.03`, 1x | 3.31e-1 | 4.60e-1 | 2.55e-1 | 1.14e-1 | 1.01e-1 | 1.53e-1 |
| flux spread `r >= 1.03`, 2x | 8.79e-2 | 1.77e-2 | 7.94e-2 | 2.68e-1 | 3.50e-1 | 2.15e-1 |
| marching steps, 1x / 2x | 200,000+2 / 161,323 | 60,552 / 119,856 | 56,226 / 111,058 | 50,426 / 99,156 | 55,335 / 108,818 | 44,878 / 87,697 |

The closure ratio is `heat / (cool + p div v + div[r^2 rho u v]/r^2)` evaluated
from the written profile with the solver's own `heat` and `cool` columns, the
form item (AA) tabulates; the total-energy form of section 9.2 of
`docs/p23_thermal_budget.md` gives the same number to 0.01-0.1 on every state
here and is not tabulated separately. The flux spread is `(max - min)/|mean|` of
`rho v r^2` over the radii shown, the gate's own measure moved inward from its
`r >= 1.2` window.

**Two rungs carry a first-cell artifact at 1x, and refining moves it rather
than removing it.** At 1x the He/H = 0.3 rung has the cold cell 1 and the
alternating temperature of section 14.4(b); at 2x that rung is clean (T
monotone, 1104 to 1061 K) and **He/H = 10 has become the alternating one**
(977, 939, 958, 1016, 1057, 1075 K, reproduced by its second cycle at
1006, 945, 951, 1001, 1047, 1071 K). The quantities the first cell sets follow
it: peak n(H3+) falls by 5.4 on He/H = 0.3 (4.05e5 to 7.54e4) and rises by 1.8
on He/H = 10 (874 to 1556), and `x2` at cell 1 moves by 4.9 percent and 23
percent on those two while the other four move by 0.01 to 0.6 percent. This
confirms section 14.4(b) and the He/H = 30 half of `p23_thermal_budget.md`
section 8.1 -- and shows the artifact is not a property of one rung.

**Two of the p23 refined runs are reproduced exactly.** The He/H = 30 2x run
here converged at the same step (87,697) with the same `||R||` (5.343e-6), flux
spread (4.936e-3), `log10 Mdot` (10.30), `x2` at cell 1 (4.593e-4), `v(1)`
(-153.2 cm/s), T maximum (7164 K) and ionization front (1.3680 R_p) that
`p23_thermal_budget.md` section 8.1 records; the He/H = 0.3 2x run reached WENO3
at step 33,772 and the JFNK at 119,856, both as section 8.2 records, with the
same first-six-cell temperatures and peak n(H3+).

#### 14.7.3 The third grid, and what it says about the second

He/H = 0.3 and the gate rung also carry a 4x state. The criterion set for this
campaign is that a quantity counts as converged from the 2x grid on when its
2x to 4x change is less than a third of its 1x to 2x change.

| quantity | He/H = 0.3: 1x / 2x / 4x | ratio | gate: 1x / 2x / 4x | ratio |
|---|---|---:|---|---:|
| log10 Mdot | 10.29 / 10.26 / 10.27 | 0.33 | 10.29 / 10.25 / 10.32 | 1.75 |
| `f(H2) = 0.5` [R_p] | 1.0330 / 1.0308 / 1.0318 | 0.44 | 1.0431 / 1.0375 / 1.0470 | 1.70 |
| `f(H2) = 1e-4` [R_p] | 1.1138 / 1.0922 / 1.0995 | 0.34 | 1.1690 / 1.1390 / 1.1759 | 1.23 |
| n(H3+) peak [cm^-3] | 4.053e5 / 7.535e4 / 7.549e4 | 0.00 | 1.760e5 / 1.803e5 / 1.680e5 | 2.86 |
| H3+ peak at [R_p] | 1.0002 / 1.0022 / 1.0024 | 0.09 | 1.0044 / 1.0051 / 1.0030 | 3.28 |
| T maximum [K] | 3685 / 3728 / 3712 | 0.35 | 2810 / 2847 / 2797 | 1.34 |
| T minimum [K] | 312 / 341 / 387 | 1.52 | 424 / 318 / 459 | 1.32 |
| ionization front [R_p] | 2.0788 / 2.0296 / 2.0361 | 0.13 | 2.3762 / 2.3047 / 2.4404 | 1.90 |
| `x_HII` at the boundary | 0.8582 / 0.8630 / 0.8609 | 0.44 | 0.8156 / 0.8239 / 0.8051 | 2.27 |
| `v` at the boundary [km/s] | 8.079 / 8.026 / 8.024 | 0.04 | 8.839 / 8.789 / 8.714 | 1.48 |
| `n_e` maximum [cm^-3] | 9.49e7 / 9.70e7 / 9.61e7 | 0.43 | 8.42e7 / 8.89e7 / 8.23e7 | 1.40 |
| `x2` at cell 1 | 0.98833 / 0.94034 / 0.94027 | 0.00 | 0.95876 / 0.95845 / 0.95837 | 0.29 |
| `\|v(1)\|` [cm/s] | 353.9 / 124.6 / 56.9 | 0.30 | 220.0 / 103.9 / 31.8 | 0.62 |
| closure at 1.03 | 0.39 / 0.97 / 0.75 | 0.38 | 0.52 / 1.00 / 0.39 | 1.25 |
| flux spread `r >= 1.03` | 4.60e-1 / 1.77e-2 / 1.23e-1 | 0.24 | 3.31e-1 / 8.79e-2 / 4.42e-1 | 1.45 |

**The two rungs do not agree, and that is the finding.** On He/H = 0.3 nine of
the fifteen quantities meet the criterion and the six that do not are mostly
quantities whose 2x to 4x move is tiny in absolute terms -- `log10 Mdot` moves
0.01 dex, which is the printed digit, and `f(H2) = 0.5` moves 0.001 R_p, which
is 0.1 percent. The genuine exception on that rung is the layer temperature
minimum, 312 to 341 to 387 K, which is still climbing. **On the gate rung only
`x2` at cell 1 meets it**: every other quantity moves further from 2x to 4x
than from 1x to 2x, and fifteen of the seventeen reverse direction (the
exceptions being the outer wind speed and `|v(1)|`, both of which keep falling),
so the sequence is not converging in `dr` at all.

**`|v(1)|` is the one quantity that behaves like a discretization error on both
rungs and on all six**: 353.9, 124.6, 56.9 cm/s on He/H = 0.3 (ratios 2.84 and
2.19) and 220.0, 103.9, 31.8 on the gate (2.12 and 3.27), i.e. first order in
`dr` or a little better, confirming section 7.2 of
`docs/p44_base_sawtooth.md` on the whole ladder rather than on one
configuration. It fails the campaign's own ratio test (0.30 and 0.62) only
because that test is written for a quantity that converges to a limit, and this
one converges to zero.

#### 14.7.4 The layer is unsteady, and that is why the sequence does not converge

Each `info = 0` state was restarted under **pure time marching** -- the JFNK
hand-off and the `du` stop both put out of reach (`Solver: Newton 1.0e-12`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-12`) -- and marched a fixed number of further
steps, one run per sample. `rho v r^2` at r = 1.03 R_p, normalized by its value
at r = 2 R_p (where the wind is flat), is the P54 deficit:

| He/H | grid | samples at 0, 500, ..., 5000 further steps | min | max | swing |
|---|---|---|---:|---:|---:|
| 0.3 | 1x | 0.562 0.723 0.588 0.597 0.642 0.655 0.621 0.611 0.634 0.640 0.625 | 0.562 | 0.723 | **0.161** |
| 0.3 | 2x | 1.006 1.132 0.957 0.937 0.948 0.991 1.026 1.032 1.015 0.998 0.993 | 0.937 | 1.132 | **0.195** |
| 30 | 1x | 0.850 0.940 1.353 1.114 1.009 1.146 1.156 1.141 1.153 1.154 1.146 | 0.850 | 1.353 | **0.503** |
| 30 | 2x | 0.790 0.870 0.819 1.029 1.129 1.085 0.973 0.882 0.861 0.898 0.951 | 0.790 | 1.129 | **0.339** |

**Every one of the four states leaves the value it was certified at.** The state
the gates accepted is one sample of a swing of 16 to 50 percent of the wind
flux, and 5000 steps is about two of the roughly 2100-step periods P54
measures. **Base refinement does not remove the oscillation and does not
systematically damp it**: the swing grows by 21 percent at He/H = 0.3 and falls
by 33 percent at He/H = 30. What the refinement does change on He/H = 0.3 is the
mean level, from about 0.63 to about 1.00 -- i.e. the 37 percent time-averaged
deficit at 1x is gone at 2x while the swing about it is not.

That is the reading for the whole of 14.7.2 and 14.7.3. **The closure ratio and
the sub-1.2 flux spread are not measures of discretization error on these
states; they are samples of an unsteady layer**, and their scatter across grids
(gate rung: closure at 1.03 of 0.52, 1.00, 0.39) is of the same size as their
scatter across restart cycles at one grid. The cycle scatter was measured on
every state:

| He/H | 1x cycle 1 | 1x cycle 2 | 2x cycle 1 | 2x cycle 2 |
|---|---|---|---|---|
| 0.0793 | 0.55 / 0.52 / 0.72 | 0.56 / 0.52 / 0.71 | 0.91 / 1.00 / 1.12 | 0.87 / 0.87 / 0.98 |
| 0.3 | 0.35 / 0.39 / 0.63 | 0.41 / 0.43 / 0.66 | 0.98 / 0.97 / 1.00 | **2.55 / 1.36 / 1.01** |
| 1 | 0.52 / 0.49 / 0.75 | **6.01 / 1.06 / 0.88** | 0.84 / 0.83 / 0.99 | **2.86 / 1.27 / 1.01** |
| 3 | 0.72 / 0.79 / 1.12 | 0.77 / 0.81 / 1.04 | 0.55 / 0.65 / 2.02 | 0.52 / 0.57 / 0.94 |
| 10 | 0.84 / 1.18 / 0.80 | 0.86 / 0.99 / 1.02 | 0.33 / 0.37 / 0.65 | 0.33 / 0.38 / 0.65 |
| 30 | **1.43 / 1.68 / 0.93** | **0.92 / 0.96 / 1.00** | 1.01 / 1.36 / 0.92 | 0.78 / 0.92 / 1.00 |

(closure at 1.02 / 1.03 / 1.05; all twenty-four states are `info = 0` on both
gates.) **This corrects one reading in item (AA).** That item compares
"0.92 / 0.96 / 1.00" at 1x against "1.01 / 1.36 / 0.92" at 2x on He/H = 30 and
concludes that refinement spoils that rung. The 1x figure it uses is the second
cycle; the first cycle of the same 1x state reads 1.43 / 1.68 / 0.93, and the
second cycle of the 2x state reads 0.78 / 0.92 / 1.00. The cycle-to-cycle spread
on that rung (0.51 at r = 1.02 on the 1x grid, 0.23 on the 2x) is comparable to
or larger than the grid-to-grid difference, so **"refinement spoils He/H = 30"
is not supported by these numbers.** The He/H = 0.3 half of the same comparison
does survive, but weakly: 0.39 at 1x against 0.97 and 1.36 on the two 2x cycles.
**The one rung whose improvement is reproducible across both cycles is the gate
rung**, 0.52 and 0.52 at 1x against 1.00 and 0.87 at 2x -- which is the
measurement item (AA) asked for as its cheapest next step, and it comes out on
the side of the discretization reading at that rung and that pair of grids. Its
4x state then reads 0.39, so even there the answer does not persist to a third
grid.

**What does reproduce across cycles.** `log10 Mdot` repeats to the printed digit
on six of the eight refined-grid pairs; the two that do not are the gate rung at
2x (10.25 then 10.27) and He/H = 3 at 2x (10.28 then 10.29), i.e. one in the
last digit. The He/H = 10 rate of 10.34 -- 0.05 dex above its 1x value, the
largest rate move on the ladder -- is identical in both cycles, so that move is
the grid and not the cycle. `f(H2) = 0.5` repeats to six digits at He/H = 0.3
and 1 on the 2x grid (1.03079 and 1.02099 twice) and on the 4x (1.03177 twice);
the ionization front repeats exactly at He/H = 0.3, 1 and 10, and T maximum to
the printed digit at 0.3 and 1 and to 1 K at 10.

#### 14.7.5 What this does and does not say about the input files

The candidate this campaign was set up to decide is whether the molecular
regression and diagnostic-arm inputs should carry
`Base grid [dr,cells]: 1.0e-4 100` and `Grid cells: 550` **as a property of the
input file** -- section 136 and `docs/p44_base_sawtooth.md` section 9.5 having
already decided that a grid must not be gated on a physics key, since that would
make the molecules-on and molecules-off runs of the same planet use different
grids and break Gate 0 of `docs/lower_atmosphere_coupling.md`.

**The case for writing it in is narrower than section 9.5 anticipated, and it
does not rest on grid convergence.** What the refinement demonstrably buys is:
the first-cell velocity artifact halves, first order in `dr`, on all six rungs;
the gate rung converges from cold in 88 minutes where at 1x it did not converge
in 84; and the 37 percent time-averaged base flux deficit on He/H = 0.3
disappears. What it does not buy is a grid-converged solution -- the gate rung's
own quantities move further from 2x to 4x than from 1x to 2x -- nor a steady
base layer, nor a smaller oscillation.

**It also relocates rather than removes the first-cell artifact**, from
He/H = 0.3 at 1x to He/H = 10 at 2x. A regression matrix pinned to 2x would
therefore be pinning a different artifact, not an artifact-free state.

The cost, measured here rather than quoted: +19 percent per step, exactly 2.0
times the steps per unit physical time, 2.24 to 2.43 times the wall clock to
convergence. For the regression matrix specifically, section 9.5's two
consequences stand unchanged and were not re-measured -- all six `mol_*` cases
are 12,000-step relaxation snapshots whose physical duration would halve, and
the interior cell count changes 500 to 550, so the goldens cannot be diffed
against the old ones and have to be re-snapshotted.

**No prescription is chosen here**, and this section does not recommend one. The
finding it hands to that decision is that the base layer is time dependent on
every grid tested, so a grid choice cannot be justified as "the converged one";
it can only be justified by what it demonstrably fixes, which on this evidence
is the cell-1 velocity artifact, the cold-start convergence of the gate rung,
and the mean base flux deficit of one rung.

The figure is `docs/lower_atmosphere_figs/fig_ladder_grid.png`: `f(H2)`,
`n(H3+)` and `T` against radius with one line style per grid, the cell-1
velocity against the realized base cell size, the normalized mass flux that a
steady state would make flat, and the r = 1.03 sample series of 14.7.4. It is
drawn by `make_fig_grid.py` in the campaign scratch directory beside the runs;
it is not regenerated by `docs/lower_atmosphere_figs/make_figures.py`.

#### 14.7.6 What was not measured

* **The 4x grid was run on two rungs only**, the gate and He/H = 0.3. Whether
  the He/H = 10 alternation that appears at 2x survives a further refinement, or
  moves to another rung again, is not known.
* **The oscillation was sampled on four states**, He/H = 0.3 and 30 at 1x and
  2x, at 500-step resolution over 5000 steps. No 4x state was sampled, no
  period was fitted, and the gate rung -- the one every published number rests
  on -- was not sampled.
* **The samples are single realizations.** Each is one run of a fixed step
  count; the swing quoted is max minus min over eleven samples, not an amplitude
  from a fitted mode.
* **Nothing here explains why the He/H = 10 rate moves 0.05 dex** with the grid
  while the other five move 0.00 to 0.04, nor why the alternation lands on that
  rung at 2x.
* **The carrier transport is off on every run here**, as in section 14, so this
  is the local-equilibrium closure measured on three grids and not a transported
  ladder.
* **The regression matrix was not re-run on any refined grid.** The cost figures
  for it in 14.7.5 are section 9.5's, not re-measured here.

### 14.8 The ladder with the caloric EOS and reaction heat

> **SUPERSEDED BY SECTION 14.9 (marked 2026-09-04).** These seven states were
> finished under a binary that predates block I (sections 148-155 of
> `docs/Update_EXHALE.md`). Section 14.9 re-finishes the same seven states
> under the current tree binary (md5 `063b2673`) and quotes what changed. The
> physical rows below reproduce there to the digit they are printed to; the
> rows that do NOT carry over are the reported `||R||` (block I evaluates it at
> the state's own composition, and it rises by 3 to 16 times), the base
> velocity (the sign of the cell-1 velocity flips on every rung), and the three
> `rho v r^2` spread rows (block I's gate reads the Riemann face flux, and the
> deep flux hole these states carry in cells 2 to 4 is gone). Read this section
> for the P53 physics and section 14.9 for the numbers.

> **PROVISIONAL, under the current residual norm (marked 2026-09-03).** Every
> state below was accepted by the two gates of section 133 as the norm stands
> in the tree at 17:38 KST. Whether that norm should be volume weighted is
> under review: on this planet a G23 root reads a `rho v r^2` spread of
> 3.4e-3 over `r >= 1.03`, which is not flat, and if the review changes the
> weighting these states will be re-finished once more. The runs and the
> scripts are kept for that. Read the numbers as measured under the present
> norm, not as final.

**Judgment first. The composition-dependent caloric EOS and the molecular
reaction heat change the molecular layer and leave the wind almost alone, and
the reaction heat is the larger of the two by a wide margin -- it is 71 percent
of the total heating at its peak, not a correction. The layer gets hotter and
its temperature minimum moves out: on the gate rung the minimum rises from
424 K at 1.066 to 541 K at 1.076 `r_base`, and it rises on five of the six
rungs. `gamma_eff` at the base runs from 1.409 at He/H = 0.0793 to 1.667 at
He/H = 30, i.e. the correction is proportional to how much H2 the rung still
carries, and it is back to 5/3 by 1.10 `r_base` on every rung. `log10 Mdot`
moves UP on every rung, by 0.01 to 0.06 dex, and the ladder stays flat:
10.30 to 10.34 across a factor 378 in He/H, a spread of 0.04 dex against
0.01 dex before. On the P23 comparison the direction does not change -- our H2
is still gone by 1.05 `r_base` where Koskinen's is not -- but the temperature
deficit at the front narrows: at 1.15, 1.20 and 1.30 `r_base` we are now 734,
453 and 152 K below Koskinen where before P53 we were 550, 257 and -6 K, and
at 1.05 the deficit falls from 686 K to 547 K. Beyond the front nothing moves:
the 3.00 `r_base` gap is 2578 K against 2616 K.**

> **THE FLUX GATE HAS SINCE BEEN REDEFINED (marked 2026-09-03, P55).** Every
> `rho v r^2` spread in this section is the CELL-CENTRE product, which is what
> the gate measured when these runs were made. `flux_spread_of_state` now
> measures the **Riemann FACE mass flux** instead, because that is what the
> finite-volume scheme conserves; section 11.3 of `docs/p55_base_mode.md`
> measures the two functionals differing by six orders on an accepted state
> (face flux flat to 3.2e-10 above 1.10 `r_base` where the centred product read
> 2.5e-4). The numbers below are therefore still correct AS A CELL-CENTRE
> DIAGNOSTIC -- and that is the honest reading of the paragraph that follows,
> which is about the layer and not about the gate -- but the row labelled "the
> gate's own" no longer names the acceptance test. Re-finishing the gate rung
> under the P55 binary returns a face-flux spread of 9.7e-13 on a state whose
> centred spread is 2.5e-4.

**What "converged" covers, and what it does not.** Every rung below is
`info = 0` by both gates of section 133, and both gates are gates on the WIND.
P54 measured that the layer below 1.15 `r_base` is not steady -- a 30 percent
mass-flux deficit at 1.03 `r_base`, carried by a damped base sound wave -- and
neither gate sees it; section 14.7 reaches the same conclusion from the other
side, by sampling that flux under time marching. The three `rho v r^2` spread
rows below make it explicit: the gate's own window (`r >= 1.2`) reads 2.5e-4 to
3.4e-4 on every rung, and the same measure over `r >= 1.03` reads 3.4e-3 to
1.4e-2, one to two orders larger. **Read the wind values as steady to the gate
and the layer values as not.** Cell 1 is excluded from all three windows: P54
measured its cell-center flux at -21 F_0, an artifact, while from
`r >= 1.002` the cell-center and face fluxes agree to 1 to 4 percent.

Measured 2026-09-03 in `.../scratchpad/ladder3/`. Two binaries are involved and
the distinction matters for every number here.

* The **cold starts** used the tree's `EXHALE.x` as built at 14:20 KST
  (md5 `ec1f367c` is NOT this one; this one is `ea9e7492`), the binary that
  carries P53 and passes the block-G goldens 10/10.
* Each of those seven `info = 0` states was then **re-finished** under the
  17:38 KST binary (md5 `ec1f367c`, P54 + G23 + P55, 10/10 identical), handed
  straight to the JFNK with `Load IC? True`, `Solver: Newton 100.0`,
  `du_th [PLM,WENO3]: 1.0e9 1.0e-3` and `Secondary_ionization: Immediate`.
  **All seven closed on the first restart**, in 6 seconds to 3 m 53 s.
  The table quotes the re-finished states.

Six rungs on the default grid plus the He/H = 0.3 rung on a base grid refined
by a factor 2 (`Base grid [dr,cells]: 1.0e-4 100`, `Grid cells: 550`), the
grid section 8.2 of `docs/p23_thermal_budget.md` introduced. `OMP_NUM_THREADS=4`,
all seven launched at once. No metals, no `opacity.inp`, no carrier transport;
`Secondary_ionization` at its STAGED default for the cold start. The inputs are
those of section 14.1 taken as they stand.

#### 14.8.1 The rungs

Rows carry the meaning section 14.2 gives them. `log10 Mdot` is the number the
run prints. The three spread rows are cell-center `rho v r^2` measured over the
window named, ghosts dropped and cell 1 dropped, as the paragraph above sets
out; the first of them is the gate's own window and reproduces the gate's
printed value.

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 | 0.3 (2x base grid) |
|---|---:|---:|---:|---:|---:|---:|---:|
| state | `heh0p0793_p54_r1` | `heh0p3_p54_r1` | `heh1_p54_r1` | `heh3_p54_r1` | `heh10_p54_r1` | `heh30_p54_r1` | `heh0p3_x2grid_p54_r1` |
| n(JFNK), info | 2, 0 | 0 | 0, 0 | 2, 2, 0 | 0 | 0 | 0, 0 |
| final `\|\|R\|\|` (section 133 scale) | 9.20e-06 | 5.92e-06 | 8.19e-06 | 4.63e-06 | 7.38e-06 | 5.60e-06 | 3.80e-06 |
| flux spread, the gate's own | 2.50e-04 | 2.55e-04 | 2.72e-04 | 2.80e-04 | 2.87e-04 | 2.90e-04 | 3.36e-04 |
| du of the solution | 7.55e-04 | 1.61e-04 | 7.04e-04 | 1.92e-03 | 1.22e-03 | 2.05e-04 | 8.96e-04 |
| rho v r^2 spread, r >= 1.20 (the gate window) | 2.50e-04 | 2.55e-04 | 2.72e-04 | 2.80e-04 | 2.87e-04 | 2.90e-04 | 3.36e-04 |
| rho v r^2 spread, r >= 1.10 | 2.37e-03 | 2.79e-03 | 1.34e-03 | 3.87e-04 | 2.87e-04 | 4.94e-04 | 3.21e-03 |
| **rho v r^2 spread, r >= 1.03** | **3.43e-03** | **4.42e-03** | **8.35e-03** | **1.27e-02** | **1.43e-02** | **1.06e-02** | **4.78e-03** |
| **log10 Mdot [g/s]** | **10.34** | **10.32** | **10.31** | **10.30** | **10.30** | **10.30** | **10.32** |
| cold-start marching steps (P53 stage) | 209913 | 52275 | 186856 | 43284 | 52040 | 44648 | 103185 |
| cold-start wall clock, 4 threads | 1 h 58 m 24 s | 0 h 26 m 00 s | 1 h 49 m 41 s | 0 h 20 m 44 s | 0 h 24 m 59 s | 0 h 20 m 25 s | 0 h 53 m 14 s |
| re-finish marching steps | 2002 | 2 | 2 | 4002 | 2 | 2 | 2 |
| re-finish wall clock, 4 threads | 0 h 01 m 23 s | 0 h 00 m 07 s | 0 h 00 m 18 s | 0 h 03 m 52 s | 0 h 00 m 06 s | 0 h 00 m 06 s | 0 h 00 m 14 s |
| x2 at the ghost (the handoff) | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 |
| x2 at cell 1 | 0.9579 | 0.9393 | 0.8847 | 0.7414 | 0.3477 | 0.0004463 | 0.9392 |
| H2 front `f(H2) = 0.5` [R_p] | 1.0449 | 1.0352 | 1.0212 | 1.0086 | below cell 1 | below cell 1 | 1.0352 |
| `f(H2) = 1e-2` [R_p] | 1.0838 | 1.0662 | 1.0394 | 1.0204 | 1.0064 | below cell 1 | 1.0662 |
| `f(H2) = 1e-4` [R_p] | 1.1851 | 1.1298 | 1.0928 | 1.0700 | 1.0539 | 1.0423 | 1.1298 |
| T at the H2 front [K] | 665 | 797 | 995 | 1104 | -- | -- | 797 |
| H3+ peak [R_p] | 1.0021 | 1.0002 | 1.0002 | 1.0002 | 1.0002 | 1.0002 | 1.0001 |
| n(H3+) at the peak [cm^-3] | 1.46e+05 | 7.08e+04 | 2.20e+04 | 5.81e+03 | 8.46e+02 | 3.54e-01 | 7.09e+04 |
| n(H3+) > 1 cm^-3 out to [R_p] | 1.1023 | 1.0821 | 1.0481 | 1.0237 | 1.0095 | never > 1 | 1.0824 |
| ionization front `x_HII = 0.5` [R_p] | 2.4804 | 2.1189 | 1.7761 | 1.5579 | 1.4259 | 1.3720 | 2.1193 |
| `x_HII` at the outer boundary | 0.803 | 0.852 | 0.895 | 0.917 | 0.926 | 0.928 | 0.852 |
| n_e at cell 1 [cm^-3] | 7.73e+05 | 1.68e+06 | 5.61e+06 | 1.85e+07 | 4.32e+07 | 5.33e+07 | 1.68e+06 |
| n_e maximum [cm^-3] (at [R_p]) | 8.04e+07 (1.113) | 9.26e+07 (1.232) | 1.21e+08 (1.154) | 1.50e+08 (1.093) | 1.71e+08 (1.066) | 1.72e+08 (1.054) | 9.24e+07 (1.230) |
| T minimum [K] (at [R_p]) | 541 (1.0759) | 584 (1.0710) | 655 (1.0606) | 705 (1.0481) | 762 (1.0320) | 860 (1.0190) | 584 (1.0711) |
| T minimum beyond cell 6 [K] (at [R_p]) | 541 (1.0759) | 584 (1.0710) | 655 (1.0606) | 705 (1.0481) | 762 (1.0320) | 860 (1.0190) | 584 (1.0711) |
| T at 1.15 R_p [K] | 983 | 1643 | 3119 | 4426 | 4714 | 4129 | 1643 |
| T at 1.50 R_p [K] | 2683 | 3648 | 5021 | 6220 | 6915 | 7095 | 3648 |
| T at 2.00 R_p [K] | 2687 | 3330 | 4335 | 5237 | 5785 | 5971 | 3330 |
| T at 3.00 R_p [K] | 2164 | 2570 | 3215 | 3781 | 4123 | 4242 | 2570 |
| T maximum [K] (at [R_p]) | 2787 (1.695) | 3655 (1.538) | 5042 (1.442) | 6278 (1.416) | 6961 (1.431) | 7114 (1.458) | 3655 (1.543) |
| outflow speed at the boundary [km/s] | 8.83 | 8.12 | 7.30 | 6.79 | 6.54 | 6.46 | 8.10 |
| **gamma_eff minimum** (at [R_p]) | **1.4088** (1.0002) | **1.4605** (1.0002) | **1.5444** (1.0002) | **1.6168** (1.0002) | **1.6583** (1.0002) | **1.6667** (1.0002) | **1.4604** (1.0001) |
| gamma_eff at cell 1 | 1.4088 | 1.4605 | 1.5444 | 1.6168 | 1.6583 | 1.6667 | 1.4604 |
| heat_mol_chem / heat_total at 1.02 | 0.706 | 0.705 | 0.647 | 0.564 | 0.122 | 0.032 | 0.705 |
|   the same at 1.05 | 0.548 | 0.488 | 0.177 | 0.064 | 0.015 | 0.013 | 0.487 |
|   the same at 1.10 | 0.068 | 0.026 | 0.014 | 0.010 | 0.008 | 0.011 | 0.026 |
|   its maximum (at [R_p]) | 0.707 (1.0241) | 0.708 (1.0164) | 0.705 (1.0060) | 0.703 (1.0002) | 0.663 (1.0002) | 0.359 (1.0002) | 0.708 (1.0164) |
| energy-floor hits | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| positivity limiter: face states | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| positivity: ghost cells / face fluxes | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| non-root acceptances (all ledgers) | 0 | 38 | 5 | 2 | 1 | 0 | 66 |
| projected/handback roots | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| cells with no in-simplex start | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| largest element-budget violation | 2.0e-01 | 0.0e+00 | 5.6e-03 | 2.6e-04 | 0.0e+00 | 0.0e+00 | 0.0e+00 |
| trial states refused (section 138.3) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Read down the `gamma_eff` row and the three reaction-heat rows together: both
are functions of how much H2 the rung still has at that radius, and both fall
along the ladder in step with `x2 at cell 1`. At He/H = 30, whose first cell
holds 4.5e-4 of its H nuclei in H2, `gamma_eff` is 1.6667 to four decimals and
the reaction heat is 3 percent of the total at 1.02 `r_base`; at the gate rung
it is 1.409 and 71 percent.

#### 14.8.2 What P53 changed, rung by rung

The section 14.2 states re-measured with the same code beside the new ones, so
the two entries of every pair are like for like. The pre-P53 column is the
section 14.2 cycle-1 state (the cold run on five rungs, the restart on the gate
rung); the 2x-grid pair takes its pre-P53 side from the section 8.2 run of
`docs/p23_thermal_budget.md`.

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 | 0.3 (2x base grid) |
|---|---:|---:|---:|---:|---:|---:|---:|
| (pre -> p54) |  |  |  |  |  |  |  |
| log10 Mdot [g/s] | 10.29 -> 10.34 | 10.29 -> 10.32 | 10.29 -> 10.31 | 10.28 -> 10.30 | 10.29 -> 10.30 | 10.29 -> 10.30 | 10.26 -> 10.32 |
| rho v r^2 spread, r >= 1.20 | 4.98e-03 -> 2.50e-04 | 2.58e-03 -> 2.55e-04 | 1.78e-03 -> 2.72e-04 | 2.43e-03 -> 2.80e-04 | 4.27e-03 -> 2.87e-04 | 3.59e-03 -> 2.90e-04 | 4.67e-03 -> 3.36e-04 |
| rho v r^2 spread, r >= 1.10 | 4.15e-02 -> 2.37e-03 | 5.34e-03 -> 2.79e-03 | 5.12e-03 -> 1.34e-03 | 4.74e-03 -> 3.87e-04 | 1.03e-02 -> 2.87e-04 | 1.56e-02 -> 4.94e-04 | 1.22e-02 -> 3.21e-03 |
| rho v r^2 spread, r >= 1.03 | 3.31e-01 -> 3.43e-03 | 4.60e-01 -> 4.42e-03 | 2.55e-01 -> 8.35e-03 | 1.14e-01 -> 1.27e-02 | 1.01e-01 -> 1.43e-02 | 1.53e-01 -> 1.06e-02 | 1.77e-02 -> 4.78e-03 |
| `f(H2) = 0.5` front [R_p] | 1.0431 -> 1.0449 | 1.0330 -> 1.0352 | 1.0202 -> 1.0212 | 1.0099 -> 1.0086 | below cell 1 -> below cell 1 | below cell 1 -> below cell 1 | 1.0308 -> 1.0352 |
| `f(H2) = 1e-2` [R_p] | 1.0799 -> 1.0838 | 1.0643 -> 1.0662 | 1.0397 -> 1.0394 | 1.0224 -> 1.0204 | 1.0080 -> 1.0064 | below cell 1 -> below cell 1 | 1.0531 -> 1.0662 |
| `f(H2) = 1e-4` [R_p] | 1.1690 -> 1.1851 | 1.1138 -> 1.1298 | 1.0882 -> 1.0928 | 1.0637 -> 1.0700 | 1.0518 -> 1.0539 | 1.0465 -> 1.0423 | 1.0922 -> 1.1298 |
| x2 at cell 1 | 0.9588 -> 0.9579 | 0.9883 -> 0.9393 | 0.8869 -> 0.8847 | 0.7453 -> 0.7414 | 0.3538 -> 0.3477 | 0.0004563 -> 0.0004463 | 0.9403 -> 0.9392 |
| n(H3+) peak [cm^-3] | 1.760e+05 -> 1.456e+05 | 4.053e+05 -> 7.081e+04 | 2.258e+04 -> 2.199e+04 | 5.967e+03 -> 5.813e+03 | 8.744e+02 -> 8.459e+02 | 3.639e-01 -> 3.540e-01 | 7.535e+04 -> 7.086e+04 |
| H3+ peak [R_p] | 1.0044 -> 1.0021 | 1.0002 -> 1.0002 | 1.0002 -> 1.0002 | 1.0002 -> 1.0002 | 1.0002 -> 1.0002 | 1.0002 -> 1.0002 | 1.0022 -> 1.0001 |
| T minimum of the layer [K] | 424 -> 541 | 471 -> 584 | 629 -> 655 | 636 -> 705 | 747 -> 762 | 946 -> 860 | 341 -> 584 |
|   its radius [R_p] | 1.0656 -> 1.0759 | 1.0590 -> 1.0710 | 1.0566 -> 1.0606 | 1.0419 -> 1.0481 | 1.0320 -> 1.0320 | 1.0248 -> 1.0190 | 1.0418 -> 1.0711 |
| T at 1.15 R_p [K] | 1167 -> 983 | 1995 -> 1643 | 3304 -> 3119 | 4701 -> 4426 | 4730 -> 4714 | 4118 -> 4129 | 2397 -> 1643 |
| T at 1.50 R_p [K] | 2753 -> 2683 | 3685 -> 3648 | 5033 -> 5021 | 6258 -> 6220 | 7011 -> 6915 | 7206 -> 7095 | 3722 -> 3648 |
| T at 2.00 R_p [K] | 2671 -> 2687 | 3303 -> 3330 | 4322 -> 4335 | 5227 -> 5237 | 5837 -> 5785 | 6045 -> 5971 | 3283 -> 3330 |
| T at 3.00 R_p [K] | 2126 -> 2164 | 2529 -> 2570 | 3195 -> 3215 | 3761 -> 3781 | 4159 -> 4123 | 4296 -> 4242 | 2497 -> 2570 |
| T maximum [K] | 2810 -> 2787 | 3685 -> 3655 | 5064 -> 5042 | 6344 -> 6278 | 7070 -> 6961 | 7231 -> 7114 | 3728 -> 3655 |
| ionization front [R_p] | 2.3762 -> 2.4804 | 2.0788 -> 2.1189 | 1.7574 -> 1.7761 | 1.5312 -> 1.5579 | 1.4054 -> 1.4259 | 1.3541 -> 1.3720 | 2.0296 -> 2.1193 |
| n_e maximum [cm^-3] | 8.420e+07 -> 8.040e+07 | 9.492e+07 -> 9.255e+07 | 1.231e+08 -> 1.214e+08 | 1.542e+08 -> 1.499e+08 | 1.779e+08 -> 1.710e+08 | 1.732e+08 -> 1.720e+08 | 9.697e+07 -> 9.244e+07 |
| v at the boundary [km/s] | 8.84 -> 8.83 | 8.08 -> 8.12 | 7.28 -> 7.30 | 6.80 -> 6.79 | 6.62 -> 6.54 | 6.57 -> 6.46 | 8.03 -> 8.10 |

**Four things move and two do not.**

* **The rate moves up and stays flat.** Every rung gains 0.01 to 0.06 dex, and
  the 2x-grid rung gains 0.06. The ladder's own spread goes 0.01 -> 0.04 dex,
  which is still at the printed precision, so "the rate does not depend on the
  helium abundance" survives P53.
* **The layer gets hotter and its minimum moves out.** The layer temperature
  minimum rises on five rungs (424 -> 541 K on the gate rung, 471 -> 584 at
  He/H = 0.3, 629 -> 655 at 1, 636 -> 705 at 3, 747 -> 762 at 10) and falls on
  one (946 -> 860 at He/H = 30, the rung with no molecular layer). Its radius
  moves out by 0.004 to 0.012 `R_p` on the same five. This is the direction the
  reaction heat implies and the sign is the same on every rung that has H2.
* **The H2 layer gets slightly deeper and the fronts move out.** `f(H2) = 1e-4`
  moves out on five rungs, by 0.002 to 0.016 `R_p`, and the `f = 0.5` front by
  up to 0.002. The molecular layer is not removed or created by P53; it is
  displaced.
* **The ionization front moves out on every rung**, by 0.02 to 0.10 `R_p`.
* **The wind does not move.** T maximum falls by 0.4 to 1.6 percent, the
  outflow speed at the boundary changes by at most 1.7 percent, and T at 2.00
  and 3.00 `r_base` moves by less than 2 percent on every rung. Whatever P53
  does, it does it below 1.2 `r_base`.
* **The base composition does not move.** `x2` at cell 1 changes in the fourth
  digit on five rungs. The one exception is He/H = 0.3 on the default grid
  (0.9883 -> 0.9393), and that is the section 14.4(b) grid artifact
  disappearing, not P53: the 2x-grid rung reads 0.9403 -> 0.9392, i.e. it was
  already at the artifact-free value before P53 and stays there.

**The grid sensitivity that section 8.2 of `p23_thermal_budget.md` recorded is
gone.** There the default and 2x base grids gave `log10 Mdot` 10.29 and 10.26,
a 0.03 dex difference, and the He/H = 0.3 rung's peak `n(H3+)` differed by a
factor 5.4. On the re-finished states the two grids give 10.32 and 10.32, peak
`n(H3+)` 7.08e4 against 7.09e4 cm^-3, layer temperature minimum 584 K against
584 K, and ionization front 2.1189 against 2.1193 `R_p`. **The likely reading
is that most of what section 8.2 measured as a grid effect was the convergence
measure rather than the discretization**, since the two states agree once both
are finished on the face-flux measure. That reading is not isolated here -- the
two changes (P53 and the P54 measure) arrived together and no run separates
them -- and section 14.7 measured the layer to be unsteady on every grid, which
is the competing explanation.

#### 14.8.3 What the re-finish changed

The same seven states before and after the P54 re-finish, i.e. the same physics
under two convergence measures. This is the column that says how much of
section 14.8.2 is P53 and how much is the measure.

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 | 0.3 (2x base grid) |
|---|---:|---:|---:|---:|---:|---:|---:|
| (p53 -> p54) |  |  |  |  |  |  |  |
| log10 Mdot [g/s] | 10.36 -> 10.34 | 10.33 -> 10.32 | 10.34 -> 10.31 | 10.37 -> 10.30 | 10.33 -> 10.30 | 10.29 -> 10.30 | 10.36 -> 10.32 |
| rho v r^2 spread, r >= 1.20 | 2.96e-03 -> 2.50e-04 | 1.33e-03 -> 2.55e-04 | 2.20e-03 -> 2.72e-04 | 3.83e-03 -> 2.80e-04 | 2.64e-03 -> 2.87e-04 | 2.98e-04 -> 2.90e-04 | 2.41e-03 -> 3.36e-04 |
| rho v r^2 spread, r >= 1.10 | 5.23e-02 -> 2.37e-03 | 2.72e-02 -> 2.79e-03 | 5.45e-02 -> 1.34e-03 | 3.51e-02 -> 3.87e-04 | 1.43e-02 -> 2.87e-04 | 1.41e-03 -> 4.94e-04 | 1.13e-01 -> 3.21e-03 |
| rho v r^2 spread, r >= 1.03 | 2.57e-01 -> 3.43e-03 | 6.27e-01 -> 4.42e-03 | 2.11e-01 -> 8.35e-03 | 1.25e-01 -> 1.27e-02 | 2.63e-01 -> 1.43e-02 | 3.40e-02 -> 1.06e-02 | 4.27e-01 -> 4.78e-03 |
| `f(H2) = 0.5` front [R_p] | 1.0488 -> 1.0449 | 1.0361 -> 1.0352 | 1.0189 -> 1.0212 | 1.0014 -> 1.0086 | below cell 1 -> below cell 1 | below cell 1 -> below cell 1 | 1.0386 -> 1.0352 |
| `f(H2) = 1e-2` [R_p] | 1.0929 -> 1.0838 | 1.0679 -> 1.0662 | 1.0389 -> 1.0394 | 1.0180 -> 1.0204 | 1.0054 -> 1.0064 | below cell 1 -> below cell 1 | 1.0696 -> 1.0662 |
| `f(H2) = 1e-4` [R_p] | 1.2106 -> 1.1851 | 1.1451 -> 1.1298 | 1.0654 -> 1.0928 | 1.0280 -> 1.0700 | 1.0105 -> 1.0539 | 1.0388 -> 1.0423 | 1.1638 -> 1.1298 |
| x2 at cell 1 | 0.9579 -> 0.9579 | 0.9392 -> 0.9393 | 0.8852 -> 0.8847 | 0.7444 -> 0.7414 | 0.3497 -> 0.3477 | 0.0004459 -> 0.0004463 | 0.9396 -> 0.9392 |
| n(H3+) peak [cm^-3] | 1.459e+05 -> 1.456e+05 | 7.067e+04 -> 7.081e+04 | 2.201e+04 -> 2.199e+04 | 1.496e+05 -> 5.813e+03 | 8.501e+02 -> 8.459e+02 | 3.536e-01 -> 3.540e-01 | 7.112e+04 -> 7.086e+04 |
| H3+ peak [R_p] | 1.0021 -> 1.0021 | 1.0002 -> 1.0002 | 1.0002 -> 1.0002 | 1.0025 -> 1.0002 | 1.0002 -> 1.0002 | 1.0002 -> 1.0002 | 1.0001 -> 1.0001 |
| T minimum of the layer [K] | 589 -> 541 | 621 -> 584 | 860 -> 655 | 461 -> 705 | 1084 -> 762 | 864 -> 860 | 708 -> 584 |
|   its radius [R_p] | 1.0997 -> 1.0759 | 1.0853 -> 1.0710 | 1.0876 -> 1.0606 | 1.0025 -> 1.0481 | 1.0044 -> 1.0320 | 1.0213 -> 1.0190 | 1.1038 -> 1.0711 |
| T at 1.15 R_p [K] | 803 -> 983 | 1401 -> 1643 | 2268 -> 3119 | 2839 -> 4426 | 4022 -> 4714 | 4106 -> 4129 | 1080 -> 1643 |
| T at 1.50 R_p [K] | 2621 -> 2683 | 3626 -> 3648 | 4950 -> 5021 | 6047 -> 6220 | 6816 -> 6915 | 7112 -> 7095 | 3555 -> 3648 |
| T at 2.00 R_p [K] | 2688 -> 2687 | 3333 -> 3330 | 4360 -> 4335 | 5300 -> 5237 | 5814 -> 5785 | 5977 -> 5971 | 3348 -> 3330 |
| T at 3.00 R_p [K] | 2179 -> 2164 | 2575 -> 2570 | 3261 -> 3215 | 3893 -> 3781 | 4174 -> 4123 | 4243 -> 4242 | 2612 -> 2570 |
| T maximum [K] | 2767 -> 2787 | 3639 -> 3655 | 4951 -> 5042 | 6047 -> 6278 | 6824 -> 6961 | 7131 -> 7114 | 3595 -> 3655 |
| ionization front [R_p] | 2.5168 -> 2.4804 | 2.1326 -> 2.1189 | 1.8249 -> 1.7761 | 1.6461 -> 1.5579 | 1.4641 -> 1.4259 | 1.3674 -> 1.3720 | 2.1834 -> 2.1193 |
| n_e maximum [cm^-3] | 7.881e+07 -> 8.040e+07 | 9.246e+07 -> 9.255e+07 | 1.181e+08 -> 1.214e+08 | 1.399e+08 -> 1.499e+08 | 1.651e+08 -> 1.710e+08 | 1.722e+08 -> 1.720e+08 | 8.978e+07 -> 9.244e+07 |
| v at the boundary [km/s] | 8.84 -> 8.83 | 8.13 -> 8.12 | 7.33 -> 7.30 | 6.85 -> 6.79 | 6.57 -> 6.54 | 6.47 -> 6.46 | 8.12 -> 8.10 |

**The re-finish is not cosmetic below the front.** It tightens the layer by one
to two orders on the `r >= 1.03` window (2.6e-1 -> 3.4e-3 on the gate rung,
6.3e-1 -> 4.4e-3 at He/H = 0.3, 4.3e-1 -> 4.8e-3 on the 2x grid) and it moves
quantities that the layer sets: T at 1.15 `r_base` by +180 to +1587 K, the
layer temperature minimum by -322 to +244 K, and on He/H = 3 the peak `n(H3+)`
by a factor 26 (1.50e5 -> 5.81e3 cm^-3, the peak moving 1.0025 -> 1.0002 `R_p`).
`log10 Mdot` falls by 0.01 to 0.07 dex, which is where the 0.08 dex spread of
the P53-only states goes: **the ladder is flat on the re-finished states and
was not on the states before them.** Above the front the re-finish does almost
nothing -- T at 2.00 and 3.00 `r_base` move by under 3 percent, the outflow
speed by under 1 percent.

The reading this supports is that **the P53-only states were phase samples of
the section 14.7 oscillation and the re-finished ones are much closer to its
mean**, not that P53 itself is uncertain: the pre-P53 states were finished on
the old measure too, so the pre/post pairing of 14.8.2 compares like with like
only in the sense that both sides carry the same kind of layer error. The
cleanest statement available from these runs is the one in the judgment: the
reaction heat is 71 percent of the heating at its peak and the layer responds;
how much of the 0.01-0.06 dex rate change survives a further norm change is
not settled here.

#### 14.8.4 Against Koskinen et al. (2022) Model A

The temperature comparison of section 14.5, with the P54 gate-rung state and
the P53-only state beside it. The published column is the same digitization.

| r/r_base | 1 | 1.05 | 1.1 | 1.15 | 1.2 | 1.3 | 1.5 | 2 | 3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **T [K], ours (P54)** | 1117 | 633 | 592 | 983 | 1448 | 2115 | 2683 | 2687 | 2164 |
| T [K], ours (P53 only) | 1117 | 811 | 589 | 803 | 1246 | 1968 | 2621 | 2688 | 2179 |
| T [K], ours before P53 (section 14.5) | 1098 | 494 | 626 | 1167 | 1644 | 2273 | 2753 | 2671 | 2126 |
| T [K], Koskinen Model A | 1082 | 1180 | 1566 | 1717 | 1901 | 2267 | 2871 | 4120 | 4742 |
| deficit T(K22) - ours, before P53 | -16 | 686 | 940 | 550 | 257 | -6 | 118 | 1449 | 2616 |
| deficit, P53 only | -35 | 369 | 977 | 914 | 655 | 299 | 250 | 1432 | 2563 |
| **deficit, P54** | -35 | 547 | 974 | 734 | 453 | 152 | 188 | 1433 | 2578 |
| ours/Koskinen, P54 | 1.03 | 0.54 | 0.38 | 0.57 | 0.76 | 0.93 | 0.93 | 0.65 | 0.46 |

**The deficit narrows through the front and nowhere else.** At 1.15, 1.20 and
1.30 `r_base` it falls from 550, 257 and -6 K to 734, 453 and 152 K -- which is
a WIDENING at those three radii, not a narrowing, and the narrowing is at 1.05
(686 -> 547 K). The honest summary is that P53 moved our profile's shape in the
layer rather than lifting it: our 1.05 `r_base` temperature rises 494 -> 633 K
while our 1.15 falls 1167 -> 983 K, because the reaction heat is deposited
below the front and the front itself moved out. Beyond 1.5 `r_base` nothing
changes: 1449 -> 1433 K at 2.00 and 2616 -> 2578 K at 3.00.

So the section 14.5 conclusion stands unchanged: **raising our base to
Koskinen's own level and then adding the two P53 channels still leaves our H2
gone by 1.05 `r_base` where theirs survives to 3, and the outer temperature gap
is untouched.** The P23 attribution -- that the gap is a composition effect
through 2 `r_base` and both a composition and an energy effect beyond -- is not
re-decided by these runs; section 10 of `docs/p23_thermal_budget.md` carries
the budget re-measured on this state.

#### 14.8.5 The two channels separated

**The two P53 channels are now separated on the gate rung, and the reaction
heat carries essentially all of the effect.** The re-finished He/H = 0.0793
state was re-finished once more with `Molecular reaction heat: False` and
everything else identical (`Load IC? True`, `Solver: Newton 100.0`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-3`, `Secondary_ionization: Immediate`, the same
P54 binary `ec1f367c`). It closed in 7.9 s, `info = 0` on both gates,
`||R|| = 5.087e-6`, gate flux spread 2.468e-4.

The caloric EOS has no input key -- `caloric_state_from_composition` takes the
legacy constant-gamma branch only when the run has no molecular chemistry -- so
the pair below is **caloric EOS alone** against **caloric EOS plus reaction
heat**, not "P53 off" against "P53 on". `gamma_eff` is the same 1.409 in both
columns, which is the check that the EOS is indeed on in B.

| | B: caloric EOS only | C: EOS + reaction heat | C - B |
|---|---:|---:|---:|
| JFNK `info`, both gates | 0 | 0 | -- |
| final `\|\|R\|\|` | 5.087e-06 | 9.198e-06 | 4.111e-06 |
| flux spread, gate window r >= 1.20 | 2.468e-04 | 2.501e-04 | 3.360e-06 |
| rho v r^2 spread, r >= 1.10 | 9.48e-04 | 2.37e-03 | 1.43e-03 |
| rho v r^2 spread, r >= 1.03 | 8.21e-03 | 3.43e-03 | -4.78e-03 |
| **log10 Mdot [g/s]** | 10.29 | 10.34 | 0.05 |
| heat_mol_chem/heat_total at 1.02 | 0.000 | 0.706 | 0.706 |
| gamma_eff minimum | 1.4090 | 1.4088 | -0.0002 |
| layer T minimum [K] | 372.1 | 541.4 | 169.3 |
|   its radius [r_base] | 1.0522 | 1.0759 | 0.0237 |
| H2 front `f = 0.5` [r_base] | 1.0394 | 1.0449 | 0.0055 |
| H2 front `f = 0.1` [r_base] | 1.0527 | 1.0664 | 0.0137 |
| `f = 1e-2` [r_base] | 1.0698 | 1.0838 | 0.0140 |
| `f = 1e-4` [r_base] | 1.1447 | 1.1851 | 0.0404 |
| T at the `f = 0.5` front [K] | 428.1 | 665.1 | 237.0 |
| H3+ peak radius [r_base] | 1.0037 | 1.0021 | -0.0015 |
| n(H3+) at the peak [cm^-3] | 1.539e+05 | 1.456e+05 | -8.305e+03 |
| ionization front [r_base] | 2.3762 | 2.4804 | 0.1042 |
| T at 1.05 r_base [K] | 374 | 633 | 259 |
| T at 1.10 r_base [K] | 788 | 592 | -196 |
| T at 1.15 r_base [K] | 1384 | 983 | -401 |
| T at 1.50 r_base [K] | 2793 | 2683 | -110 |
| T at 2.00 r_base [K] | 2663 | 2687 | 24 |
| T at 3.00 r_base [K] | 2108 | 2164 | 56 |
| T maximum [K] | 2827 | 2787 | -41 |
| outflow speed at the boundary [km/s] | 8.79 | 8.83 | 0.04 |

| Koskinen deficit T(K22) - ours [K] | B | C | C - B |
|---|---:|---:|---:|
| 1.05 r_base | 806 | 547 | -259 |
| 1.10 r_base | 778 | 974 | +196 |
| 1.15 r_base | 333 | 734 | +401 |
| 1.20 r_base | 68 | 453 | +385 |
| 1.30 r_base | -125 | 152 | +277 |
| 1.50 r_base | 78 | 188 | +110 |
| 2.00 r_base | 1457 | 1433 | -24 |
| 3.00 r_base | 2634 | 2578 | -56 |

**Reading.** `log10 Mdot` moves 10.29 -> 10.34, so **the whole +0.05 dex is the
reaction heat**; the caloric EOS on its own leaves the rate at the pre-P53
value of 10.29. The layer temperature minimum moves +169 K (372 -> 541 K) and
out by 0.024 `r_base`, the H2 fronts move out by 0.006 to 0.040 `r_base`, and
the ionization front by 0.104. Against Koskinen the two channels pull in
opposite directions over the front: the reaction heat closes the 1.05 `r_base`
deficit by 259 K and **opens** the 1.15, 1.20 and 1.30 deficits by 401, 385 and
277 K. That is the same redistribution section 141.4 of `docs/Update_EXHALE.md`
records, measured here on states finished under the face-flux measure.

**Against the section 141.4 table (pre-G23).** That A/B/C was measured on the
pre-G23 residual, restarting a state from the earlier campaign; the C column
here reproduces its C closely -- `log10 Mdot` 10.34 in both, layer minimum
541 K against 525 K, `f = 0.5` front 1.0449 against 1.0429, peak `n(H3+)`
1.456e5 against 1.469e5 cm^-3 -- while the B column differs by more (minimum
372 K against 417 K, `f = 0.1` front 1.0527 against 1.0598), so the reaction
heat's own increment on the rate agrees to 0.01 dex (+0.05 here, +0.046 there)
and its increment on the layer minimum does not (+169 K here, +108 K there).

#### 14.8.6 What was not measured

* **The channel split was measured on the gate rung only** (section 14.8.5).
  The other five rungs and the 2x-grid rung were not run with
  `Molecular reaction heat: False`, so their reaction-heat fraction rows say
  how large that channel is, not what each run would do without it.
* **The two changes arrived together.** Every P53 state was cold-started under
  one binary and re-finished under another, so "P53" and "the face-flux measure"
  are not separated except by section 14.8.3's before/after, which is itself a
  comparison of two finishes of the same physics. The two PHYSICS channels are
  separated, on one rung, in section 14.8.5.
* **The residual norm is under review** (the box at the head of this section).
  If it is volume weighted these seven states are superseded.
* **The carrier transport is off on every rung**, as in sections 14 and 14.7.
* **Only the He/H = 0.3 rung was run on a refined base grid.** The grid
  statement of 14.8.2 rests on that one pair; section 14.7 has the wider grid
  material, measured before the re-finish.
* **`log10 Mdot` is quoted to the two decimals the run prints**, and the whole
  ladder spread is 0.04 dex, so the flatness claim is again a statement that
  the differences are near what is printed.

The figure is `docs/lower_atmosphere_figs/fig_ladder_p53.png`: `f(H2)` and `T`
against `r/r_base` with the pre-P53 states dashed and the two digitized
published profiles as points, `gamma_eff` of the caloric EOS, and the share of
the total heating the molecular reaction heat carries, one curve per rung. It
is drawn by `make_fig.py` in `.../scratchpad/ladder3/`; it is not regenerated
by `docs/lower_atmosphere_figs/make_figures.py`.

### 14.9 The ladder under the block-I physics

**Judgment first. Block I (sections 148-155 of `docs/Update_EXHALE.md`) leaves
the physics of this ladder exactly where section 14.8 left it and fixes the
base. All seven rungs re-finish to `info = 0`, and on every one of them the
wind and the layer reproduce section 14.8 to the digit the table prints:
`log10 Mdot` is identical on all seven (10.34, 10.32, 10.31, 10.30, 10.30,
10.30, 10.32), the H2 fronts agree to 0.0004 `r_base`, the layer temperature
minimum to 1 K, the ionization front to 1e-4 `r_base`, `gamma_eff` to four
decimals, and `T` at 1.15, 1.50 and 3.00 `r_base` to at most 6 K. What DID
change is the base. The characteristic condition of section 152 removes the
spurious inflow at the first cell -- the cell-1 velocity flips from -15, -37,
-104, -298, -801 and -1089 cm/s to +14, +13, +10, +8, +8 and +7 cm/s on the
six default-grid rungs, and from -1.0 to +11.3 on the 2x-grid rung -- and with
it the deep mass-flux hole those states carried two to four cells above the
base: the smallest `rho v r^2` anywhere outside cell 1, in units of the wind
value, rises from 0.73, 0.50, 0.13, -0.10, -0.02, 0.66 and 0.86 to 0.99, 0.95,
0.99, 0.99, 0.99, 0.99 and 0.98. At He/H = 3 and He/H = 10 the old states had
a REVERSED cell-centre mass flux at 1.0006 `r_base`; they do not any more. The
face mass flux the gate now reads is flat to 1e-14 to 1e-11 in its own window
and to 1e-10 to 1e-9 over `r >= 1.03`, against the 3e-3 to 1.4e-2 the
cell-centre product still reads over the same window -- and that row is
unchanged to three digits, so the cell-centre non-flatness of a few parts in a
thousand ABOVE 1.03 `r_base` is untouched. The reported `||R||` rises by 3 to
16 times (9.2e-6 to 1.5e-4 on the gate rung), which is section 155's
own-composition norm doing what section 155 said it would; it changes no
acceptance here, because none of these inputs sets `Resid tol` and the
residual gate is off in `steady_gates_met`. What it DOES change is how hard
the JFNK works: five rungs close on the first restart as in section 14.8, but
He/H = 3 needs four rounds and He/H = 10 needs nine, and on the way the JFNK's
own stagnation, not the flux gate, is what refuses the state (section
14.9.3).**

Measured 2026-09-04 in `.../scratchpad/ladder4/`. One binary: the tree's
`EXHALE.x` as built 2026-09-04 02:23 KST, md5 `063b2673`, copied to the
scratch directory before the campaign and its md5 checked there. The seven
`info = 0` states of section 14.8 were handed straight to the JFNK with `Load
IC? True`, `Solver: Newton 100.0`, `du_th [PLM,WENO3]: 1.0e9 1.0e-3` and
`Secondary_ionization: Immediate`, restarted from their own output until both
gates were met, `OMP_NUM_THREADS=4`, all seven launched at once. Nothing else
in the inputs was touched.

**What acceptance means in these runs, and which test is which.** Three
different thresholds appear in these logs and confusing them makes the
convergence story unreadable:

* **The JFNK's own convergence test**, `||R|| < 1e-5 AND flux spread < 2e-5`.
  Failing it is what prints `(JFNK) gate NOT met: ...`. A JFNK that stagnates
  before meeting it prints `returning best iterate` and hands back `info = 2`.
* **`steady_gates_met` in `EXHALE_main`**, which is what the acceptance line
  reports: `||R|| < -1.000E+00 AND flux spread < 2.000E-05`. The negative
  residual threshold is `resid_th = -1`, i.e. **the residual gate is OFF**
  because no input here sets `Resid tol`; the binding test is the flux gate,
  and since P55 it reads the Riemann FACE mass flux over `r >= r_flux = 1.2`.
* **`du`**, which is not a gate here at all: `du_th` is 1.0e9 for the PLM
  stage, so the marching never stops on it and the Newton is armed at once.

So section 155's own-composition residual inflates the number these runs
REPORT and decides no acceptance, exactly as section 155.2 predicted for a
case with `Resid tol` unset -- but it is compared against the JFNK's own 1e-5,
and there it does bite. Cell 1 is excluded from every cell-centre spread
below, as in section 14.8.

#### 14.9.1 The rungs

Rows carry the meaning section 14.2 and section 14.8.1 give them, with three
additions: the three face-flux rows are the run's own printed windows (the
conserved functional the gate reads), the three cell-centre rows are the same
product section 14.8 tabulated, so the two sections can be compared directly,
and the cell-1 velocity row is section 152's signature.

| He/H | 0.0793 | 0.3 | 1 | 3 | 10 | 30 | 0.3 (2x base grid) |
|---|---:|---:|---:|---:|---:|---:|---:|
| state | `heh0p0793_r1` | `heh0p3_r1` | `heh1_r1` | `heh3_r4` | `heh10_r9` | `heh30_r1` | `heh0p3_x2grid_r1` |
| n(JFNK), info | 2, 2, 0 | 0 | 0, 0 | 0 | 0 | 0 | 0, 0 |
| accepted by the two gates | yes | yes | yes | yes | yes | yes | yes |
| final `\|\|R\|\|` (own-composition norm) | 1.49e-04 | 9.61e-05 | 5.25e-05 | 4.90e-05 | 3.47e-05 | 1.53e-05 | 1.48e-05 |
| **face-flux spread, the gate window `r >= 1.2`** | **1.13e-13** | **4.22e-14** | **2.38e-14** | **9.35e-15** | **4.57e-12** | **1.03e-10** | **1.37e-13** |
| face-flux spread, `r >= 1.10` | 5.13e-13 | 3.10e-13 | 1.48e-13 | 6.96e-14 | 3.32e-11 | 7.20e-10 | 1.96e-12 |
| face-flux spread, `r >= 1.03` | 1.54e-10 | 1.26e-10 | 1.11e-09 | 7.34e-12 | 4.91e-10 | 2.38e-09 | 2.74e-10 |
| cell-centre `rho v r^2` spread, `r >= 1.20` | 2.49e-04 | 2.55e-04 | 2.72e-04 | 2.80e-04 | 2.87e-04 | 2.90e-04 | 3.35e-04 |
| cell-centre `rho v r^2` spread, `r >= 1.10` | 2.37e-03 | 2.77e-03 | 1.32e-03 | 3.82e-04 | 2.87e-04 | 4.90e-04 | 3.20e-03 |
| cell-centre `rho v r^2` spread, `r >= 1.03` | 3.42e-03 | 4.41e-03 | 8.34e-03 | 1.26e-02 | 1.42e-02 | 1.05e-02 | 4.80e-03 |
| du of the solution | 2.72e-04 | 1.61e-04 | 7.04e-04 | 3.89e-04 | 2.01e-04 | 2.05e-04 | 8.96e-04 |
| **log10 Mdot [g/s]** | **10.34** | **10.32** | **10.31** | **10.30** | **10.30** | **10.30** | **10.32** |
| re-finish marching steps | 4002 | 2 | 2 | 2 | 2 | 2 | 2 |
| re-finish wall clock, 4 threads | 0 h 03 m 06 s | 0 h 00 m 14 s | 0 h 00 m 19 s | 0 h 00 m 10 s | 0 h 00 m 05 s | 0 h 00 m 04 s | 0 h 00 m 29 s |
| **velocity at cell 1 [cm/s]** | **+14.39** | **+12.95** | **+9.539** | **+8.292** | **+7.704** | **+7.065** | **+11.3** |
| x2 at the ghost (the handoff) | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 | 0.985448 |
| x2 at cell 1 | 0.9577 | 0.939 | 0.8841 | 0.7399 | 0.343 | 0.000431 | 0.9391 |
| H2 front `f(H2) = 0.5` [R_p] | 1.0450 | 1.0352 | 1.0212 | 1.0085 | below cell 1 | below cell 1 | 1.0352 |
| `f(H2) = 1e-2` [R_p] | 1.0837 | 1.0661 | 1.0393 | 1.0203 | 1.0063 | below cell 1 | 1.0661 |
| `f(H2) = 1e-4` [R_p] | 1.1847 | 1.1296 | 1.0927 | 1.0698 | 1.0538 | 1.0422 | 1.1297 |
| T at the H2 front [K] | 663 | 796 | 994 | 1102 | -- | -- | 796 |
| H3+ peak [R_p] | 1.0021 | 1.0002 | 1.0002 | 1.0002 | 1.0002 | 1.0002 | 1.0001 |
| n(H3+) at the peak [cm^-3] | 1.43e+05 | 7.00e+04 | 2.17e+04 | 5.73e+03 | 8.23e+02 | 3.36e-01 | 7.04e+04 |
| n(H3+) > 1 cm^-3 out to [R_p] | 1.1023 | 1.0821 | 1.0481 | 1.0237 | 1.0095 | never > 1 | 1.0824 |
| ionization front `x_HII = 0.5` [R_p] | 2.4804 | 2.1189 | 1.7761 | 1.5579 | 1.4259 | 1.3720 | 2.1193 |
| `x_HII` at the outer boundary | 0.804 | 0.852 | 0.895 | 0.917 | 0.926 | 0.928 | 0.852 |
| n_e at cell 1 [cm^-3] | 7.80e+05 | 1.70e+06 | 5.67e+06 | 1.87e+07 | 4.34e+07 | 5.33e+07 | 1.69e+06 |
| n_e maximum [cm^-3] (at [R_p]) | 8.04e+07 (1.112) | 9.26e+07 (1.232) | 1.21e+08 (1.154) | 1.50e+08 (1.093) | 1.71e+08 (1.066) | 1.72e+08 (1.054) | 9.25e+07 (1.230) |
| T minimum [K] (at [R_p]) | 540 (1.0759) | 584 (1.0710) | 655 (1.0606) | 704 (1.0481) | 762 (1.0320) | 859 (1.0187) | 584 (1.0711) |
| T minimum beyond cell 6 [K] (at [R_p]) | 540 (1.0759) | 584 (1.0710) | 655 (1.0606) | 704 (1.0481) | 762 (1.0320) | 859 (1.0187) | 584 (1.0711) |
| T at 1.15 R_p [K] | 986 | 1646 | 3124 | 4428 | 4720 | 4133 | 1645 |
| T at 1.50 R_p [K] | 2684 | 3649 | 5022 | 6220 | 6916 | 7096 | 3648 |
| T at 2.00 R_p [K] | 2686 | 3330 | 4335 | 5236 | 5785 | 5971 | 3330 |
| T at 3.00 R_p [K] | 2164 | 2570 | 3215 | 3781 | 4123 | 4242 | 2570 |
| T maximum [K] (at [R_p]) | 2787 (1.695) | 3655 (1.538) | 5043 (1.442) | 6277 (1.416) | 6963 (1.431) | 7114 (1.458) | 3655 (1.543) |
| outflow speed at the boundary [km/s] | 8.83 | 8.12 | 7.30 | 6.79 | 6.54 | 6.46 | 8.10 |
| **gamma_eff minimum** (at [R_p]) | **1.4089** (1.0002) | **1.4607** (1.0002) | **1.5446** (1.0002) | **1.6170** (1.0002) | **1.6584** (1.0002) | **1.6667** (1.0002) | **1.4604** (1.0001) |
| gamma_eff at cell 1 | 1.4089 | 1.4607 | 1.5446 | 1.6170 | 1.6584 | 1.6667 | 1.4604 |
| heat_mol_chem / heat_total at 1.02 | 0.707 | 0.705 | 0.646 | 0.560 | 0.121 | 0.032 | 0.705 |
|   the same at 1.05 | 0.549 | 0.488 | 0.175 | 0.063 | 0.015 | 0.013 | 0.487 |
|   the same at 1.10 | 0.067 | 0.025 | 0.014 | 0.010 | 0.008 | 0.011 | 0.025 |
|   its maximum (at [R_p]) | 0.707 (1.0241) | 0.708 (1.0164) | 0.705 (1.0058) | 0.703 (1.0002) | 0.663 (1.0002) | 0.353 (1.0002) | 0.708 (1.0164) |
| energy-floor hits | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| positivity limiter: face states | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| positivity: ghost cells / face fluxes | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| non-root acceptances (all ledgers) | 0 | 1 | 0 | 1 | 0 | 0 | 10 |
| projected/handback roots | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| cells with no in-simplex start | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| largest element-budget violation | 1.5e-02 | 2.8e-03 | 1.7e-09 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 1.8e-06 |
| trial states refused | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

#### 14.9.2 What block I changed, rung by rung

The section 14.8 state and the block-I re-finish of the same rung, so the two
entries of every pair are like for like. `min rho v r^2 over cells 2..N` is
the smallest cell-centre mass flux anywhere above cell 1, in units of the mean
over `r >= 1.2`; it is the measure on which the base changed.

**He/H = 0.0793**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.34 | 10.34 |
| reported `\|\|R\|\|` | 9.20e-06 | 1.49e-04 |
| velocity at cell 1 [cm/s] | -15.0 | +14.4 |
| min `rho v r^2` over cells 2..N, / wind | 0.7279 | 0.9922 |
|   its radius [r_base] | 1.0006 | 1.0008 |
| cell-centre spread, r >= 1.03 | 3.43e-03 | 3.42e-03 |
| H2 front f = 0.5 [r_base] | 1.0449 | 1.0450 |
| f = 1e-4 [r_base] | 1.1851 | 1.1847 |
| layer T minimum [K] | 541 | 540 |
|   its radius [r_base] | 1.0759 | 1.0759 |
| T at 1.15 r_base [K] | 983 | 986 |
| T at 1.50 r_base [K] | 2683 | 2684 |
| T at 3.00 r_base [K] | 2164 | 2164 |
| n(H3+) at the peak [cm^-3] | 1.46e+05 | 1.43e+05 |
| ionization front [r_base] | 2.4804 | 2.4804 |
| gamma_eff minimum | 1.4088 | 1.4089 |
| heat_mol_chem/heat_total at 1.02 | 0.706 | 0.707 |

**He/H = 0.3**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.32 | 10.32 |
| reported `\|\|R\|\|` | 5.92e-06 | 9.61e-05 |
| velocity at cell 1 [cm/s] | -37.3 | +12.9 |
| min `rho v r^2` over cells 2..N, / wind | 0.4971 | 0.9491 |
|   its radius [r_base] | 1.0006 | 1.0004 |
| cell-centre spread, r >= 1.03 | 4.42e-03 | 4.41e-03 |
| H2 front f = 0.5 [r_base] | 1.0352 | 1.0352 |
| f = 1e-4 [r_base] | 1.1298 | 1.1296 |
| layer T minimum [K] | 584 | 584 |
|   its radius [r_base] | 1.0710 | 1.0710 |
| T at 1.15 r_base [K] | 1643 | 1646 |
| T at 1.50 r_base [K] | 3648 | 3649 |
| T at 3.00 r_base [K] | 2570 | 2570 |
| n(H3+) at the peak [cm^-3] | 7.08e+04 | 7.00e+04 |
| ionization front [r_base] | 2.1189 | 2.1189 |
| gamma_eff minimum | 1.4605 | 1.4607 |
| heat_mol_chem/heat_total at 1.02 | 0.705 | 0.705 |

**He/H = 1**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.31 | 10.31 |
| reported `\|\|R\|\|` | 8.19e-06 | 5.25e-05 |
| velocity at cell 1 [cm/s] | -104.4 | +9.5 |
| min `rho v r^2` over cells 2..N, / wind | 0.1282 | 0.9948 |
|   its radius [r_base] | 1.0006 | 1.0100 |
| cell-centre spread, r >= 1.03 | 8.35e-03 | 8.34e-03 |
| H2 front f = 0.5 [r_base] | 1.0212 | 1.0212 |
| f = 1e-4 [r_base] | 1.0928 | 1.0927 |
| layer T minimum [K] | 655 | 655 |
|   its radius [r_base] | 1.0606 | 1.0606 |
| T at 1.15 r_base [K] | 3119 | 3124 |
| T at 1.50 r_base [K] | 5021 | 5022 |
| T at 3.00 r_base [K] | 3215 | 3215 |
| n(H3+) at the peak [cm^-3] | 2.20e+04 | 2.17e+04 |
| ionization front [r_base] | 1.7761 | 1.7761 |
| gamma_eff minimum | 1.5444 | 1.5446 |
| heat_mol_chem/heat_total at 1.02 | 0.647 | 0.646 |

**He/H = 3**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.30 | 10.30 |
| reported `\|\|R\|\|` | 4.63e-06 | 4.90e-05 |
| velocity at cell 1 [cm/s] | -298.0 | +8.3 |
| min `rho v r^2` over cells 2..N, / wind | -0.1035 | 0.9949 |
|   its radius [r_base] | 1.0006 | 1.0100 |
| cell-centre spread, r >= 1.03 | 1.27e-02 | 1.26e-02 |
| H2 front f = 0.5 [r_base] | 1.0086 | 1.0085 |
| f = 1e-4 [r_base] | 1.0700 | 1.0698 |
| layer T minimum [K] | 705 | 704 |
|   its radius [r_base] | 1.0481 | 1.0481 |
| T at 1.15 r_base [K] | 4426 | 4428 |
| T at 1.50 r_base [K] | 6220 | 6220 |
| T at 3.00 r_base [K] | 3781 | 3781 |
| n(H3+) at the peak [cm^-3] | 5.81e+03 | 5.73e+03 |
| ionization front [r_base] | 1.5579 | 1.5579 |
| gamma_eff minimum | 1.6168 | 1.6170 |
| heat_mol_chem/heat_total at 1.02 | 0.564 | 0.560 |

**He/H = 10**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.30 | 10.30 |
| reported `\|\|R\|\|` | 7.38e-06 | 3.47e-05 |
| velocity at cell 1 [cm/s] | -801.5 | +7.7 |
| min `rho v r^2` over cells 2..N, / wind | -0.0245 | 0.9937 |
|   its radius [r_base] | 1.0006 | 1.0100 |
| cell-centre spread, r >= 1.03 | 1.43e-02 | 1.42e-02 |
| H2 front f = 0.5 [r_base] | -- | -- |
| f = 1e-4 [r_base] | 1.0539 | 1.0538 |
| layer T minimum [K] | 762 | 762 |
|   its radius [r_base] | 1.0320 | 1.0320 |
| T at 1.15 r_base [K] | 4714 | 4720 |
| T at 1.50 r_base [K] | 6915 | 6916 |
| T at 3.00 r_base [K] | 4123 | 4123 |
| n(H3+) at the peak [cm^-3] | 8.46e+02 | 8.23e+02 |
| ionization front [r_base] | 1.4259 | 1.4259 |
| gamma_eff minimum | 1.6583 | 1.6584 |
| heat_mol_chem/heat_total at 1.02 | 0.122 | 0.121 |

**He/H = 30**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.30 | 10.30 |
| reported `\|\|R\|\|` | 5.60e-06 | 1.53e-05 |
| velocity at cell 1 [cm/s] | -1089.4 | +7.1 |
| min `rho v r^2` over cells 2..N, / wind | 0.6591 | 0.9923 |
|   its radius [r_base] | 1.0006 | 1.0100 |
| cell-centre spread, r >= 1.03 | 1.06e-02 | 1.05e-02 |
| H2 front f = 0.5 [r_base] | -- | -- |
| f = 1e-4 [r_base] | 1.0423 | 1.0422 |
| layer T minimum [K] | 860 | 859 |
|   its radius [r_base] | 1.0190 | 1.0187 |
| T at 1.15 r_base [K] | 4129 | 4133 |
| T at 1.50 r_base [K] | 7095 | 7096 |
| T at 3.00 r_base [K] | 4242 | 4242 |
| n(H3+) at the peak [cm^-3] | 3.54e-01 | 3.36e-01 |
| ionization front [r_base] | 1.3720 | 1.3720 |
| gamma_eff minimum | 1.6667 | 1.6667 |
| heat_mol_chem/heat_total at 1.02 | 0.032 | 0.032 |

**He/H = 0.3 (2x base grid)**

| | section 14.8 | block I |
|---|---:|---:|
| log10 Mdot [g/s] | 10.32 | 10.32 |
| reported `\|\|R\|\|` | 3.80e-06 | 1.48e-05 |
| velocity at cell 1 [cm/s] | -1.0 | +11.3 |
| min `rho v r^2` over cells 2..N, / wind | 0.8587 | 0.9802 |
|   its radius [r_base] | 1.0003 | 1.0002 |
| cell-centre spread, r >= 1.03 | 4.78e-03 | 4.80e-03 |
| H2 front f = 0.5 [r_base] | 1.0352 | 1.0352 |
| f = 1e-4 [r_base] | 1.1298 | 1.1297 |
| layer T minimum [K] | 584 | 584 |
|   its radius [r_base] | 1.0711 | 1.0711 |
| T at 1.15 r_base [K] | 1643 | 1645 |
| T at 1.50 r_base [K] | 3648 | 3648 |
| T at 3.00 r_base [K] | 2570 | 2570 |
| n(H3+) at the peak [cm^-3] | 7.09e+04 | 7.04e+04 |
| ionization front [r_base] | 2.1193 | 2.1193 |
| gamma_eff minimum | 1.4604 | 1.4604 |
| heat_mol_chem/heat_total at 1.02 | 0.705 | 0.705 |

**The base, and nothing else.** Read the pairs down: `log10 Mdot` is identical
on all seven rungs, `gamma_eff` to four decimals, the reaction-heat fraction
to 0.004, the ionization front to 1e-4 `r_base`, the layer temperature minimum
to 1 K, `T` at three radii to 6 K. The three rows that move are the reported
`||R||` (up 3 to 16 times), the cell-1 velocity (sign flipped on all seven)
and the layer flux minimum (from -0.10 ... 0.86 to 0.95 ... 0.99). This is
what section 152 was for, and the reading it supports is that the P44 base
sawtooth and the base-adjacent part of the P54 flux deficit were both
artifacts of the old ghost prescription rather than features of the solution:
removing them moved no physical number on any rung.

**What did NOT move, and is worth naming.** The cell-centre `rho v r^2` spread
over `r >= 1.03` is 3.42e-3, 4.41e-3, 8.34e-3, 1.26e-2, 1.42e-2, 1.05e-2 and
4.80e-3 -- the same values to three digits as section 14.8. The base fix
therefore closed the hole in cells 2 to 4 and left the non-flatness of a few
parts in a thousand in the cell-centre product from 1.03 to 1.2 `r_base`
exactly as it was. The face flux over the same window reads 1e-12 to 1e-9, so
the two functionals disagree by six to nine orders there, which is the same
disagreement section 11.3 of `docs/p55_base_mode.md` measured; nothing here
decides which one a physical steadiness claim should rest on.

**A base normalization that moved and is not traced.** `ntot_bc_per_H` in
`EXHALE_resolved.out` rises from 0.5434782, 0.6209816, 0.7536380, 0.8768190,
0.9552069 and 0.9841057 to 0.5464827, 0.6250068, 0.7592713, 0.8838170,
0.9630203 and 0.9922105, i.e. by 0.55 to 0.82 percent, monotonically in He/H,
on a `q_H2_base` and an `H_nuclei_in_H2_fraction` that are unchanged
(0.9854478 on every rung, both binaries). Since `n0` is derived from `p_base =
n0 k_B T0 ntot_bc`, the base number density falls by the same fraction, which
is 0.002 to 0.004 dex on the rate and below the two decimals it is printed to.
**Which block-I section changes it was not traced**, and it is the one
input-side number in this campaign that section 14.8 does not reproduce.

#### 14.9.3 What the two slow rungs cost, and what refused them

Five rungs close on the first restart, in 4 to 29 seconds, as in section 14.8.
**He/H = 3 needs four rounds and He/H = 10 needs nine**, and the two are the
rungs whose H2 dissociation front sits in the base cells (`x2 at cell 1` is
0.74 and 0.34, against 0.96 on the gate rung). Every failing attempt on both
puts its worst cell at 1.002 to 1.006 `r_base` in the ENERGY row, which is
that front.

He/H = 10 is worth reading round by round, because it separates the two
thresholds. Twenty-seven JFNK attempts over nine rounds, three per round:

| round | `\|\|R\|\|` of the three attempts | face-flux spread of the three | `log10 Mdot` |
|---|---|---|---:|
| 1 | 1.06, 0.77, 0.80 | 2.6e-3, 8.9e-3, 3.2e-3 | 10.32 |
| 2 | 1.08, 1.20, 1.23 | 2.5e-3, 2.1e-3, 2.1e-3 | 10.32 |
| 3 | 1.80, 1.02, 0.70 | 2.2e-3, 9.6e-3, 3.3e-3 | 10.32 |
| 4 | 1.52, 0.54, 0.59 | 2.3e-3, 2.3e-3, 2.1e-3 | 10.32 |
| 5 | 0.51, 0.84, 0.90 | 2.3e-3, 2.1e-3, 2.0e-3 | 10.32 |
| 6 | 0.49, 0.50, 0.49 | 2.0e-3, 1.7e-3, 1.6e-3 | 10.32 |
| 7 | 0.49, 0.49, 0.69 | 2.0e-3, 1.7e-3, 2.0e-3 | 10.32 |
| 8 | 0.21, 0.17, **0.083** | 1.2e-3, 1.2e-3, **1.51e-5** | 10.30 |
| **9** | **3.47e-5** | **4.57e-12** | **10.30** |

**It is converging, and slowly: rounds 1 to 4 hover between 0.5 and 1.8, and
from round 5 the residual walks down from 0.90 to 3.5e-5 over four rounds.**
Two things are worth taking from the table. First, the marching between
restarts is doing the work -- each of rounds 1 to 8 marches 4003 steps and
spends 3 m 39 s to 4 m 21 s on four threads before its three JFNK attempts,
and it is the marching that walks the layer down, not the Newton; round 9
needs 2 steps and 5.5 s. Second, **on the last attempt of
round 8 the FLUX gate was already met** (1.51e-5 against the 2e-5 threshold)
and the state was still refused, because the JFNK had stagnated (`returning
best iterate`) with an own-composition residual of 8.3e-2 and therefore
returned `info = 2`, which `EXHALE_main` never accepts. So the statement "the
flux gate is what refuses these states" is true for rounds 1 to 7 and FALSE
for round 8: there the binding condition was the JFNK's own `||R|| < 1e-5`,
evaluated at the state's own composition. That is the one place in this
campaign where section 155's norm changes an outcome rather than a report, and
it is a slower convergence, not a rejected solution: round 9 closes at `||R||
= 3.47e-5` with the flux spread four orders below tolerance.

The rate is insensitive to all of this: `log10 Mdot` prints 10.32 on rounds 1
through 7 and 10.30 on rounds 8 and 9, and the accepted value 10.30 is exactly
the section 14.8 value for this rung.

#### 14.9.4 What was not measured

* **One binary, one restart recipe.** Block I arrived as one block; nothing
  here separates section 150 from 151 from 152 from 155 on this ladder. Two of
  them can still be spoken about separately, and the two statements are
  different in kind:
  * **Section 150 cannot act here.** The Lyman-Werner channel is OFF on every
    rung, and the run says so itself: `WARNING H2 Lyman-Werner
    photodissociation OFF: no "Stellar LW flux" key and no spectrum to
    integrate one from`, followed by `The band exists whenever the star does,
    so this is a MISSING CHANNEL, not a modelling choice: H2 keeps a sink it
    should not.` `heat_H2_LW` is 0.00e+00 at every radius, before and after,
    so halving the dissociation rate multiplies zero. **This ladder cannot
    test section 150 at all**; doing so needs a rung that supplies
    `Stellar LW flux` or reads a spectrum file, and that is a separate
    campaign.
  * **Section 151 IS resolved here, and its effect is genuinely small --
    which makes this a measurement rather than a gap.** The power-law branch
    of `src/modules/init/set_energy_vectors.f90` lays `num_HI = 50`
    logarithmic points between `e_th_HI` = 13.6 eV and `e_th_HeI` = 24.6 eV,
    of which **13 fall inside the 15.4-18.0 eV window Backx et al. (1976)
    measured, at a spacing of 0.185 to 0.211 eV** (15.494, 15.679, 15.866,
    ... 17.862 eV). The threshold structure is carried on the grid at that
    resolution. With it carried, `heat_H2` at 1.05 `r_base` moves 7.84e-9 to
    7.29e-9, **-7 percent**, on a channel that is 7.6 percent of the
    photoheating there, i.e. **-0.6 percent of the local photoheating**. The
    85-times factor at threshold does not become a large change in the rate
    because the channel it corrects is a small part of the budget on this
    planet, and that is a statement about this planet, not about the energy
    grid.
* **The cell-centre / face-flux disagreement is recorded, not resolved**, and
  the layer's steadiness therefore still depends on which functional is
  believed.
* **`ntot_bc_per_H` moved 0.55 to 0.82 percent and was not traced.**
* **The nine-round path of He/H = 10 was not repeated**, so whether nine is
  the number or this instance was unlucky is not established; the two rungs
  that needed more than one round are the two whose front sits in the base
  cells, which is a correlation on two points.
* **The carrier transport is off on every rung**, and no rung was run with
  `Molecular reaction heat: False`, so section 14.8.5's channel split was not
  repeated under block I.
* **Only the He/H = 0.3 rung was run on a refined base grid**, as in 14.8.
* **`log10 Mdot` is quoted to the two decimals the run prints**, and the whole
  ladder spread is 0.04 dex, so the flatness claim is again a statement that
  the differences are near what is printed.

The figure is `docs/lower_atmosphere_figs/fig_ladder_blockI.png`: `f(H2)` and
`T` against `r/r_base` with the section 14.8 state of each rung dotted and the
two digitized published profiles as points; the velocity of the first twelve
cells, which is section 152's signature -- the dotted curves swing from
-10^3 to +10^2 cm/s across cells 1 to 3 and the solid ones do not; and
`rho v r^2` normalized to its wind value against `r/r_base - 1` on a log axis,
which puts the base cells where they can be read: the dotted curves oscillate
between -0.2 and 1.3 inside `r - 1 < 1e-3` and the solid ones are flat there.
All seven rungs are `info = 0`, so every rung is solid. It is drawn by
`make_fig4.py` in `.../scratchpad/ladder4/`; it is not regenerated by
`docs/lower_atmosphere_figs/make_figures.py`.

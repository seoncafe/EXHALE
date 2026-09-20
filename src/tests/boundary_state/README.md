# Boundary-state tests

What a residual's lower boundary is a function of, and whether every consumer
of one evaluation reads the same one. Built for item D5b-1 of
`docs/PLAN_20260918_rev2.md`; the measurements behind it are
`docs/lhs1140b_stationary_D5a_20260918.md` and
`docs/lhs1140b_stationary_D5b1_20260918.md`.

Run them:

```
src/tests/boundary_state/run.sh                        # every test
src/tests/boundary_state/run.sh installed_composition  # one of them
EXHALE_EXE=<binary> EXHALE_TEST_OUT=<dir> src/tests/boundary_state/run.sh
```

Names: `installed_composition`, `declared_inputs`, `molecular_seed`,
`contact_law`.

`contact_law` is a Fortran driver and links the production objects, so it
needs an object directory as well: `EXHALE_OBJDIR=<dir>`, default `build/`.

## `installed_composition`

`boundary_of_the_installed_composition.sh`.

**The statement.** One composition, one boundary, one residual. The base face
state, the ghost cell averages and the ghost conserved state are functions of
the interior cell averages and of the installed composition, through the
caloric map the ghost continuation reads; `Apply_BC` derives them and caches
them in `BC_Apply`, and `Rec_BC` puts the cached face state straight into the
left slot of the base face. A residual whose sources, pressure map and sound
speeds belong to one composition and whose cached face state belongs to
another is the residual of no single gas.

**The measure.** The largest relative difference between the cached boundary
and the boundary `base_boundary_states` gives for the state installed at the
instant of the call, written by `boundary_state_trace_mod` under
`EXHALE_BOUNDARY_TRACE`. Zero exactly where the statement holds.

**The two runs.** A row that measures zero because nothing was measured
passes silently, so the same measure is taken a second time with
`EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition`, which leaves
the boundary of the composition the route was entered with standing. That is
the order D5b-1 corrects, and on a molecular base the measure must then be
strictly positive.

**The fixture.** `backup/regression/carrier_model_a_newton/IC/`, the
hot-Uranus molecular carrier state, evaluated with `Restart intent:
stationary evaluate` on a copy; the regression directory is never written to.
It is a relaxation snapshot and not a stationary solution, so its rows are far
from zero and its certification refuses; neither is asserted. A molecular base
is required: on an atomic mixture the ghost continuation carries no
composition, so the measure is zero in both orders and the second run would
assert nothing.

**Verdict at the item's landing** (MEASURED, D5b-1 section 7): both rows PASS;
the entry text's order reads 4.2856e-07 where the invariant asks for zero.

**Reported, not gated.** The cell-1 continuity row at the residual and after
the face-flux report, in floors of that row's own rounding floor. The report
installs a boundary of its own on a copy of the conserved array and the cached
face state carries no validity yet; making it carry one is D5b-2, and until
then the movement is on the record rather than in a verdict.

## `declared_inputs`

`boundary_from_the_declared_inputs.sh`. Built for item D5b-2; the
measurements are `docs/lhs1140b_stationary_D5b2_20260918.md`.

**The statement.** The lower boundary is a function of the physical conserved
state and composition of the column, the prescribed reservoir, the radiation
context and the model options, and of nothing else. The ghost rows of a
restart file are the operation's own output and are not read; the ghost's
molecular partition is solved against the ghost's own ionization balance
rather than prescribed from the composition the state was entered with; the
prescribed reservoir and the ghost's own particle count are two named
quantities and the solved one is never written back; every read of the cached
face state belongs to the composition installed at the read; the boundary
model and the reservoir travel with the state; and the two evaluation sites
of a solve derive the boundary after the sweep.

**The rows.** Eleven assertions, in that order. The two runs the first three
compare differ in exactly the two lower ghost species rows of the restart pair
(the second is given the composition of the first physical cell, rescaled to
the ghost's own heavy-particle count), and must agree in every data line of
both products and in the signed base continuity rows the trace writes at
sixteen digits.

**The fixture.** `backup/regression/carrier_model_a_newton`, on copies; a
molecular base is required, since the closure and the two counts have no
subject without one. Nothing is asserted about the fixture's own rows or its
certification: it is a relaxation snapshot.

**Verdict at the item's landing** (MEASURED, D5b-2): eleven of eleven PASS on
the delivered build and eleven of eleven FAIL on a build of the entry text.

## `molecular_seed`

`molecular_seed_of_an_atomic_state.sh`. Built for item D9fix; the
measurements are `docs/lhs1140b_stationary_D9fix_20260919.md`.

**The statement.** The molecular seed the catalog builds from a certified
atomic state is built: the lower ghost's base handoff partition
`x_H2 = x2 (1 - x_ion)` closes against the ghost's own ionization balance to
1e-10 in `x_H2`, and the conversion 2 H -> H2 conserves the hydrogen nuclei
of every cell, the lower ghosts included, whose composition load_IC takes
from the molecular reservoir and which therefore already carry H2.

**The rows.** Two assertions: the closure stop is not printed, and the
largest relative change of the H nuclei the seed report prints is at most
`molecular_seed_closure_tol` = 1e-12.

**The fixture.** The certified generation of
`LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13` that the 2026-09-19
catalog refresh converted, read from its published (read-only) directory, and
the input of `molecular_scalar_gj1132_kzz1e9/HeH2.13` written by the script.
On this ghost the substitution map has slope 0.74 (MEASURED), which is why
the rows need this state and not the hot-Uranus fixture, where substitution
closes in 2 to 4 passes.

**Verdict at the item's landing** (MEASURED, D9fix): two of two PASS on the
delivered build and two of two FAIL on the entry text (the closure stops at
30 passes, last move 1.19e-6, so no seed report is written).

## `contact_law`

`contact_law_at_the_base.f90`. Built for item D2b of
`docs/PLAN_20260918_rev2.md`; the measurements are
`docs/lhs1140b_stationary_D2b_20260918.md`.

**The statement.** The lower boundary is a CONTACT between the lower
atmosphere and the first interior cell, closed in three steps in this
order: the two sides are matched acoustically by the interior's outgoing
relation (C-) closed with the reservoir's pressure; the direction of the
contact is read off the matched state; the entropy and composition of the
face are then upwinded on that direction. At inflow and at a
pressure-balanced contact at rest the reservoir owns the level, the lower
atmosphere being a heat bath on the time scales of a stationary solution;
only a reverse flow carries the interior's trace out. The branch takes no
intermediate value and has no width. Each transport operator carries its
own base condition and none of them follows the contact's direction.

**The rows.** Twenty-two assertions in six groups, on columns this driver
builds with the production grid, gravity and molecular equation of state:

- a hydrostatic column on the reservoir's own isentrope, at rest: the face
  velocity, the base mass flux and the distance from the face to the
  reservoir all at the rounding floor, and the derivation repeating to the
  last bit;
- a pressure-balanced stationary contact whose interior stands at twice the
  level's temperature, so the two isentropes are a factor two apart: the
  face is the reservoir's, no flow is launched, and the branch weight is
  exactly 0;
- a pressure mismatch of either sign at zero cell velocity: the branch
  follows the sign of the matched face velocity, and the outgoing acoustic
  invariant of the interior is carried through the face unchanged;
- the two transport base conditions read on either side of a branch flip:
  the conduction base temperature is the level's own and does not move, the
  element-diffusion Dirichlet composition does not move, and the conductive
  heat flux through the base face at zero bulk velocity is MEASURED and
  printed rather than assumed zero;
- the face state `Apply_BC` installs against the one the condition returns;
- the face state against the one remaining regularization width, the
  supersonic-outflow handover, changed by four decades.

**The fixture** is built in the driver, so the row needs no run directory
and no state file. `L35`'s "rest" fixture, a wind with its momentum
deleted, is NOT a hydrostatic column and its base mass flux of 1.2751e-5 is
not this row's reference; the column here is continued from the face along
the reservoir's own isentrope, which is the state the condition is well
balanced on.

**Verdict at the item's landing** (MEASURED, D2b): twenty-two of
twenty-two PASS on the delivered build. On a build of the entry text three
rows do not compile at all, the transport operators having had no base
condition of their own, and of the nineteen that remain seven FAIL: the
stationary contact takes the average of the two isentropes instead of the
reservoir (0.25 of the face density) and its branch weight reads 1/2, both
signs of the pressure mismatch read 1/2 instead of 0 and 1, and the face
state moves by 1.5e-6 and 1.8e-6 when the width is changed by two decades
in each direction.

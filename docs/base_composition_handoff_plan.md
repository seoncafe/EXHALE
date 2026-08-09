# Plan: carry the photochemical H2/H partition into EXHALE

- Date: 2026-08-09
- Status: proposal, nothing implemented
- Motivation: `docs/vulcan_photochem_comparison.md`
- Touches: `input_read.f90`, `composition.f90`, `parameters.f90`,
  `vulcan_to_base.py`, `docs/input_schema.md`

## 1. The problem

The Tier-3 pre-step runs a photochemical kinetics code to find the composition
at the base of the escape model. The result that matters is the hydrogen
dissociation state: `lower_atmosphere_coupling.md` records that photochemistry
puts `q_H = 0.23` at 1 ubar for HD 189733 b where chemical equilibrium gives
`0.020`, about 11x more atomic hydrogen, and states this as the reason the tier
exists.

That number does not reach the code. `read_base_inp`
(`input_read.f90:942-981`) accepts four keys —

```
T_base -> T0     r_base -> R0     HeH_base -> HeH     Kzz_base -> he_kzz
```

— and `src/utils/vulcan_to_base.py` writes `q_H2` and `q_H` as comments.

Meanwhile EXHALE *does* use a base H2 fraction. `input_read.f90:779-789`:

```fortran
if (molecular_base) then
   qmb  = q_h2_equilibrium(1.0d-6, T0)
   x2mb = 2.0d0*qmb*(1.0d0 + HeH)/(1.0d0 + qmb)
   ntot_bc = ntot_bc - 0.5d0*x2mb/(1.0d0 + HeH)
endif
```

`q_h2_equilibrium` (`lower_column.f90:57-71`) is the Visscher/Koskinen
**chemical-equilibrium** fit — precisely the approximation the Tier-3 pre-step
was built to replace. So the code asks the right question at the right place and
answers it from the wrong source, while the right answer sits in a comment two
files away.

Both quantities are the same thing: `q_h2_equilibrium` returns the mixture
volume mixing ratio `n_H2/(n_H2+n_H+n_He)` (`lower_column.f90:82-83`), which is
the definition `vulcan_to_base.py` extracts. **No unit or convention conversion
is needed.**

## 2. Size of the effect

`ntot_bc` is the base particle count per (H+He) nucleus, in units of `n0`. It is
1 for a neutral atomic H/He base. The `molecular_base` correction subtracts
`q/(1+q)`.

| source of q_H2 | q_H2 | subtracted from `ntot_bc` |
|---|---|---|
| equilibrium fit (what runs today) | 0.8384 | 0.4560 |
| VULCAN / NCHO | 0.6181 | **0.3820** |
| Photochem / NCHO | 0.7170 | 0.4176 |
| Photochem / Zahnle | 0.8063 | 0.4464 |

Using VULCAN's photochemical value instead of the equilibrium fit changes
`ntot_bc` by **+0.074 per nucleus, about 7%**. For comparison, the four keys the
handoff carries today agree between codes to 0.02%; this is ~350x larger.

`ntot_bc` propagates to:

- `Apply_BC.f90:65` — the ghost-cell pressure, `W_in(3,index) = ntot_bc + dp_bc`
- `set_IC.f90:66,78,141,152,242` — base sound speed and the hydrostatic /
  transonic initial condition
- `input_read.f90:796` — in `Base BC: pressure` mode, `n0 = base_p_ubar/(kb*T0*ntot_bc)`,
  so the base **density** moves by the same ~7%

One caution against expecting a large wind response: Photochem's own network
gives 0.8063, within 1% of the equilibrium fit's 0.8384. The 7% is a VULCAN-vs-
equilibrium difference, and per the comparison memo it is the *reaction network*
that produces it, not the code. Whoever runs this should expect the size of the
effect to depend on the network chosen, not on having "switched on
photochemistry".

## 3. Three routes, and which to take

| route | what it does | scope |
|---|---|---|
| **1. EOS base particle count** | photochemical `q_H2` replaces the equilibrium fit in `molecular_base` | small, self-contained, recommended first |
| **2. Tier-2 boundary condition** | `q_H2` constrains the base cell of the solved molecular network | needs a design decision |
| **3. molecular/metal mixing ratios** | `q_H2O`, `q_CO`, ... and metal release | blocked, see §6 |

Route 1 changes an existing code path from an approximate source to a measured
one. It introduces no new physics and no new switch. Route 2 introduces a
conflict of authority that has to be settled first. They are independent:
`molecular_base` is documented as "EOS-only ... the chemistry stays atomic"
(`parameters.f90:135-139`), so route 1 can land without touching Tier-2.

## 4. Route 1 in detail

### 4.1 Call order is already right

`read_base_inp` is called at `input_read.f90:693`; the `molecular_base` block is
at 779. A value read from `base.inp` is therefore available where it is needed
with no reordering.

### 4.2 Where the policy should live

`composition.f90:100-112` declares itself the **single source of the base
composition policy** — `mass_per_H`, `ntot_bc`, `rho_bc` all flow from
`comp_mass_per_H` / `comp_ntot_bc` / `comp_rho_bc`, "so the policy cannot
disagree between code paths (the §3.4 root cause)".

The current `molecular_base` block violates that: it mutates `ntot_bc` inline in
`input_read.f90` *after* `comp_ntot_bc()` has set it. The fix should move the H2
correction into `composition.f90` rather than add a second inline branch:

```fortran
! composition.f90
real*8 function comp_ntot_bc()
   comp_ntot_bc = 1.0d0
   if (eos_include_metals .and. thereis_metals) then
      comp_ntot_bc = (1.0d0 + HeH + sum(melem_ab))/(1.0d0 + HeH)
   endif
   if (molecular_base) comp_ntot_bc = comp_ntot_bc - h2_bound_fraction()
end function

! H nuclei bound into H2 at the base, per (H+He) nucleus, from the
! photochemical handoff when one was supplied and from the
! chemical-equilibrium fit otherwise.
real*8 function h2_bound_fraction()
   real*8 :: q
   if (q_h2_base .gt. 0.0d0) then
      q = q_h2_base                       ! base.inp, photochemical
   else
      q = q_h2_equilibrium(p_base_bar, T0)
   endif
   h2_bound_fraction = q/(1.0d0 + q)
end function
```

`q/(1+q)` is the existing expression simplified: `0.5*x2mb/(1+HeH)` with
`x2mb = 2q(1+HeH)/(1+q)` reduces to it exactly. Keep the algebra in whichever
form reproduces the current arithmetic bitwise for the equilibrium branch —
verify with a metals-off, `Molecular base: True` run before and after.

### 4.3 New key

`read_base_inp` gains one branch, in the style of the four existing ones:

```fortran
else if (index(line,'q_H2_base') .gt. 0) then
   str = get_word(line,2);  read(str,*) q_h2_base
```

with `real*8 :: q_h2_base = -1.0d0` in `parameters.f90` (negative = absent, so
the equilibrium fit stays the default and old `base.inp` files behave exactly as
now).

### 4.4 Converter change

`vulcan_to_base.py` currently writes

```
# PHOTOCHEMICAL base state (VULCAN): q_H2=6.1810e-01 q_H=2.3820e-01 ...
```

It should write the comment *and* the key:

```
q_H2_base 6.1810e-01
```

`docs/compare_vulcan_photochem.py` (`base` action) writes the same file for
either source and should be updated together, so the two stay in step.

### 4.5 A consistency trap to close

`q_h2_equilibrium(1.0d-6, T0)` hard-codes the handoff pressure at 1 ubar, while
`vulcan_to_base.py` takes `--pbase` (default 1e-6 bar). If someone hands off at
a different pressure, the fit and the handoff silently describe different
levels. Two options:

- have `base.inp` carry `p_base` and use it in the fit call, or
- have `read_base_inp` reject a `base.inp` whose `p_base` is not 1 ubar.

The first is better; either is better than the present silence.

## 5. Route 2: Tier-2 boundary condition

With `Molecular chemistry: True`, H2 is a solved species (`isp_H2 = 34`, bsp
position 7) along with H2+, H3+, HeH+, and `molecular_base` is auto-enabled
(`input_read.f90:754-758`). The base partition is then an *output* of the
network, not a parameter, and imposing `q_H2` becomes a boundary condition that
competes with it.

The decision to make, stated plainly: **below 1 ubar, is the photochemical code
the authority or is EXHALE's own network?** Arguments both ways:

- the pre-step code has the fuller network, real photolysis cross sections and a
  proper radiative transfer, and it is solved on a domain that extends down to
  1 bar where the composition is anchored;
- EXHALE's network is solved consistently with the wind, the local radiation
  field and the temperature the wind solve produces, and pinning it at the base
  can fight the flow.

A middle option: use `q_H2` from the handoff only to build the initial
condition, and let the network relax it. That is cheap, cannot fight the solve,
and still removes the equilibrium fit from the *starting* state.

Nothing here should be built until that choice is made. Note also the indexing
warning in `composition.f90:117-121`: `bsp_mass` is indexed by bsp position and
the shortcut used for the six atomic species does **not** extend to `isp_H2`.

## 6. Route 3: what stays out

`vulcan_to_base.py` also records `q_H2O`, `q_CO`, `q_CH4`, `q_CO2`, `q_NH3`,
`q_HCN`. These are not actionable yet: EXHALE's metal set is atomic
(C/N/O/Mg/Si/Ca/Na/K/S/Fe) and `input_read` refuses `Molecular chemistry`
together with trace metals. Neither VULCAN nor Photochem releases atomic metals
at all, so `metals.inp` stays user-supplied regardless — that part of the
Lavvas model has no public analogue. Leave these as comments.

## 7. Validation

### 7.1 What must not move

Default-off falls out of the design rather than needing a switch: with no
`q_H2_base` key the equilibrium branch runs, so **every existing run is
byte-identical**. Confirm rather than assume:

- `backup/regression/run_check.sh check` — the two WASP-121 b cases, which have
  `Molecular base` off, must stay byte-identical;
- a `Molecular base: True` run without `q_H2_base` must be byte-identical, which
  is what tests the §4.2 refactor (moving the correction into `composition.f90`).
  The Tier-2 gate cases are the ones to use: `data_g1a`, `data_g1m`, `data_g2`
  under the gate staging directory, 12000-step convention.

### 7.2 What should move, and by how much

With `q_H2_base` present, on HD 189733 b:

- `ntot_bc` from 0.5440 to 0.6180 for a metals-free base (where it is 1 before
  the correction, which is the Tier-2 case since `Molecular chemistry` and trace
  metals cannot be combined) — +13.6% on the value, the subtraction itself
  changing by 0.074 per nucleus. With metals in the EOS budget the starting
  value is `(1+HeH+sum melem_ab)/(1+HeH)` instead and the shift is the same
  0.074 in absolute terms.
- in `Base BC: pressure` mode, `n0` down by the same factor
- the wind response is the open number — it is what the change is for

Report it the way the comparison memo does: mass flux over the code's own `du`
window, not an outer-half average.

### 7.3 Planet choice

**Do not validate on HD 189733 b.** It does not converge — `du` wandered between
0.77 and 1.86 over 87000 steps, the base-breathing item in `TO_BE_DONE.md` (A) —
so a 7% base-density change cannot be separated from the oscillation. Use a
planet that reaches a steady state (HD 209458 b, or the WASP-121 b regression
cases), which means running the pre-step for that planet first: VULCAN needs its
T(p), Kzz and stellar spectrum, and the traps in
`vulcan_photochem_comparison.md` §"Traps" apply.

## 8. Risks

1. **The refactor in §4.2 is the risky part, not the new key.** Moving the H2
   correction into `comp_ntot_bc` changes when it is applied relative to the
   metal terms. Bitwise reproduction of the equilibrium branch is the gate.
2. **`molecular_base` is auto-enabled by `thereis_mol`**, so route 1 silently
   affects every Tier-2 run. That is intended, but it means the Tier-2 gates
   move as soon as a `q_H2_base` is supplied.
3. **The effect size depends on the network, not on "using photochemistry".**
   Photochem's default set lands within 1% of the equilibrium fit. A run that
   adopts photochemistry and sees no change has not necessarily done anything
   wrong.
4. `q_H2_base` and `HeH_base` are not independent — both come from the same
   solution — but only one of them is currently checked for consistency with the
   elemental abundance. Worth a validation print at read time.

## 9. Suggested order

1. Add `q_H2_base` to `read_base_inp` + `parameters.f90`; leave the
   `molecular_base` block where it is. Verify byte-identity with no key present.
2. Move the correction into `composition.f90` (§4.2). Verify byte-identity again
   on the Tier-2 gates — this step should change nothing at all.
3. Update `vulcan_to_base.py` and `compare_vulcan_photochem.py` to emit the key;
   update `docs/input_schema.md` (a new row in the `base.inp` section) and
   `docs/lower_atmosphere_coupling.md` (Tier-3 row: the partition is now
   consumed, not just recorded).
4. Close the `p_base` trap (§4.5).
5. Run the pre-step for a converging planet and measure the wind response.
6. Only then consider route 2, and settle the authority question first.

# Plan: carry the photochemical H2/H partition into EXHALE

- Date: 2026-08-09
- Status: **Route 1 implemented and measured 2026-08-10** (§10 implementation,
  §11 wind response on HD 209458 b). §11.7 shows the measured response is a
  property of the 1 microbar handoff level rather than of photochemistry as
  such, and that a handoff deep enough to reach T_eq puts the molecular-base
  EOS correction outside its usable range; §11.8 diagnoses the Newton-finish
  failure. Routes 2 and 3 remain proposals; §5 records the direction Route 2
  should take.
- Motivation: `docs/vulcan_photochem_comparison.md`
- What comes after this plan: `docs/oxygen_chemistry_new_plan.md` (the plan of
  record, superseding the ordering in `docs/oxygen_chemistry_options.md`) —
  widening the handoff beyond `q_H2_base`, computing the partition instead of
  importing it, and the choice of photochemistry code now that `q_H2_base`
  reaches the wind. Note its verified finding on semantics: `q_H2_base` is a
  base EOS anchor (particle count and pressure normalization), not an H2
  composition pin — `set_IC` seeds H2 from the equilibrium fit independently.
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

> [2026-08-15: superseded by §10. `read_base_inp` now accepts six keys — the
> four above plus `q_H2_base` (the photochemical H2 volume mixing ratio, the
> number this memo says does not reach the code) and `p_base` (the pressure
> level the handoff describes). The `molecular_base` correction quoted just
> below no longer lives inline in `input_read.f90`: it is
> `comp_ntot_bc`/`h2_mixing_ratio_base` in
> `src/modules/functions/composition.f90`, which uses `q_H2_base` when it is
> set and falls back to the `q_h2_equilibrium` fit otherwise.]

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

A middle option was proposed here — use `q_H2` from the handoff only to build
the initial condition and let the network relax it — and it **does not survive
contact with the code** (checked 2026-08-10). The molecular network is not
integrated in time: `ionization_equilibrium.f90:544-571` root-finds the local
equilibrium of all 7-8 species in every cell at every step, and the previous
state enters only as the hybrd1 initial guess. A seeded composition is
therefore erased on the first step; the seed only selects which root is found
in the cells where hybrd1 is bistable (which is exactly what the existing
chemical-equilibrium retry seed at lines 555-566 is for). "Seed and relax" and
"impose a boundary condition" are the same two options as before, with nothing
in between.

The direction Route 2 should take instead: the reason EXHALE's network wants a
more molecular base than the photochemistry gives is that the network is
missing the physics that dissociates H2 — Lyman-Werner photodissociation was
then on the Tier-2 open list. Adding that is a fix to the network; pinning
`q_H2` at the base is a patch over its absence. Prefer the fix. (It was added
on 2026-08-13, and the next subsection records what it did.)

Nothing here should be built until that choice is made. Note also the indexing
warning in `composition.f90:117-121`: `bsp_mass` is indexed by bsp position and
the shortcut used for the six atomic species does **not** extend to `isp_H2`.

### Route 2 executed, 2026-08-13: the missing term is now in, and it is not the answer

Lyman-Werner photodissociation was implemented in the coupled network — new key
`Stellar LW flux [erg/cm2/s]:` (default 0), module
`src/modules/lower_atmosphere/lyman_werner.f90`, Draine & Bertoldi (1996)
calibration and eq. (37) self-shielding of the star-ward H2 column, 0.4 eV of
heating per dissociation (Black & Dalgarno 1977). Full description, derivation
and validation: `lower_atmosphere_coupling.md`, section "H2 photodissociation in
the Lyman-Werner bands".

Result on the hot-Uranus gate (metals off, both A and B Newton-converged,
`info = 0`, band flux 343 erg cm^-2 s^-1 for the HD 209458 orbit): **the base
composition does not move.** q_H2 at the base goes from 0.8618 to 0.8612 against
the 0.75 the case's `base.inp` carries — 0.5% of the gap. Lyman-Werner does
become the largest single H2 loss at the base (5.9e-11 s^-1 against 2.5e-11 for
H+ + H2), but the base partition is a formation-destruction balance in which
n_H/n_H2 scales as the square root of the destruction rate, so reaching
q_H2 = 0.75 would need a destruction rate 4100x larger. The unshielded band
supplies 6.0e-5 s^-1, more than enough; the base column N_H2 = 4.3e21 cm^-2
suppresses it by 1e6. The band physically cannot reach 1 microbar.

The conclusion recorded above — that the network is missing the physics that
dissociates H2 — was right in kind and wrong in identity. The photochemical
codes' extra atomic H at 1 microbar comes from catalytic cycles on O, OH and
H2O (`lower_atmosphere_coupling.md` section 1), and EXHALE's H2/H2+/H3+/HeH+
network has nowhere to put those species. **Route 1 therefore stays the only way
EXHALE can carry the photochemical base partition**, and it should be described
as a deliberate modeling choice rather than as a patch over a missing term.
What Lyman-Werner does change is the top of the molecular layer, where the
column has thinned: the H2 fraction at r = 1.10 falls by 168x, the H2->H front
moves in by 0.0018 R_p and Mdot rises by 0.005 dex.

## 6. Route 3: what stays out

`vulcan_to_base.py` also records `q_H2O`, `q_CO`, `q_CH4`, `q_CO2`, `q_NH3`,
`q_HCN`. These are not actionable yet: EXHALE's metal set is atomic
(C/N/O/Mg/Si/Ca/Na/K/S/Fe), and its own molecular network is H2/H2+/H3+/HeH+,
so a VULCAN molecule has nowhere to land. (The parser also refused
`Molecular chemistry` together with trace metals when this was written; since
2026-08-13 the two are solved in one system, but that changes nothing here.)
Neither VULCAN nor Photochem releases atomic metals
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
so a 7% base-density change cannot be separated from the oscillation.
[2026-08-15: no longer a restriction. HD 189733 b reaches a Newton-grade steady
state (`info=0`) after the JFNK line-search fix and the beta(tau)/CHIANTI-guarded
cooling of 2026-08-11; see `docs/hd189_base_checkerboard.md` §10 and
`docs/newton_scaling_and_base_wall.md`.] Use a
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

Steps 1-4 were carried out on 2026-08-10 (§10). Step 5 was carried out on
2026-08-10 for HD 209458 b (§11). Step 6 is open.

## 10. Implementation notes (2026-08-10)

### What was added, and where

| file | change |
|---|---|
| `src/modules/init/parameters.f90:142-155` | `q_h2_base = -1.0d0` (negative = no photochemical value supplied) and `p_base_bar = 1.0d-6` (handoff level [bar]) |
| `src/modules/files_IO/input_read.f90:982-987` | `read_base_inp` branches for `q_H2_base` and `p_base`, echoed like the four existing keys |
| `src/modules/files_IO/input_read.f90:992-1009` | validation echo after the read loop (§8-4): the photochemical `q_H2` is printed next to the `HeH` in effect, and a `p_base` away from 1 microbar is flagged. Print only — an inconsistent handoff is never a reason to refuse to run |
| `src/modules/functions/composition.f90:149-178` | `h2_mixing_ratio_base()` (photochemical value if supplied, chemical-equilibrium fit at `(p_base_bar, T0)` otherwise) and `h2_bound_fraction()` (the particles the H2 binding removes, per (H+He) nucleus) |
| `src/modules/functions/composition.f90:144` | `comp_ntot_bc` subtracts `h2_bound_fraction()` after the metal term — the §4.2 move |
| `src/modules/files_IO/input_read.f90:776-788` | the inline `molecular_base` block is gone; `input_read` only echoes which H2 source was used |
| `src/modules/files_IO/write_setup_report.f90:215-216` | both new variables in `parse_dump.txt` |
| `src/utils/vulcan_to_base.py:134-135`, `docs/compare_vulcan_photochem.py:287-288` | both converters write `q_H2_base` and `p_base` as read keys; the comment lines are unchanged |

`h2_bound_fraction` keeps the original arithmetic verbatim —
`x2 = 2q(1+HeH)/(1+q)`, capped at 1, then `0.5*x2/(1+HeH)` — rather than the
`q/(1+q)` simplification sketched in §4.2. The two agree only while the cap is
inactive, and reproducing the equilibrium branch bitwise was the gate.

`p_base_bar` also closes the §4.5 trap: the fit is now evaluated at the level
the handoff describes, and the default 1e-6 reproduces the former hard-coded
literal exactly.

Deliberately **not** changed: `set_IC.f90:176` and
`ionization_equilibrium.f90:556` also call `q_h2_equilibrium`, but at each
cell's *local* pressure, not at the base — they build an H2 profile and a
solver seed, not a base particle count, so the base value does not belong
there.

### Bitwise gates

| gate | command | result |
|---|---|---|
| step 1, WASP-121 b matrix | `backup/regression/run_check.sh check` | PASS, both cases (`wasp_full` 13488 steps, `wasp_he23off` 13482) |
| step 1, Tier-2 molecular | `data_g2` input in a scratch directory, `OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=12000` | `Hydro_ioniz.txt` and `Ion_species.txt` identical to the pre-change binary, whole file |
| step 2, Tier-2 molecular | same, post-refactor binary | identical to the same pre-change run; Mdot log 10.58, the documented gate value |
| step 2, WASP-121 b matrix | `backup/regression/run_check.sh check` | PASS, both cases |

The pre-change outputs were produced with the unmodified binary before any
edit, so no version control was needed to obtain the reference.

`parse_dump.txt` gains two lines, so `run_parse_corpus.sh check` reports a diff
for every case until `parse_golden` is re-snapshotted. That corpus was already
stale before this change by four lines (`legacy_hhe_rates`, `use_sec_ion`,
`use_he_rec_coupling`, `he_h_charge_exchange`); refreshing it is a separate,
deliberate act.

### The new path exercised

`data_g2` (hot Uranus, `Molecular chemistry: True`, T0 = 1140 K, He/H = 0.0793)
with a `base.inp` carrying `q_H2_base 0.6181` and `p_base 1.0e-6`:

```
   base.inp: q_H2(base) ->  0.61810
   base.inp: p_base ->  1.00E-06 bar
   base.inp: photochemical q_H2 = 0.61810 is used with He/H = 0.07930
 (input_read) Molecular base: q_H2(base, photochemical) = 0.618 -> ntot_bc = 0.618
```

`ntot_bc` = 0.6180087757 against 0.5439643595 on the equilibrium branch
(q_H2 = 0.83836 at 1 microbar, 1140 K): the subtraction moves by 0.07404 per
nucleus, the +0.074 of §2 and §7.2. This was a plumbing check of a few hundred
steps, not a science run; the wind response (§7.2, step 5) has not been
measured.

### Route 2

Settled separately and recorded in §5: the "middle option" there — seed the
initial condition from the handoff and let the network relax it — does not
survive EXHALE's chemistry architecture. The molecular network is solved as a
local equilibrium in every cell at every step, so a seeded composition is
erased immediately; the seed only selects the Newton branch in bistable cells.
The recorded direction for Route 2 is therefore to add the missing H2
photodissociation physics to the network rather than to impose a boundary
condition on it.

## 11. Step 5: wind response on HD 209458 b (2026-08-10)

Everything below is measured in this working copy unless it is labeled as a
choice or a reading. Run directories: `vulcan_work/hd209_vulcan/` (the
photochemistry) and `vulcan_work/hd209_wind_response/{eqfit,photo}/` (the A/B).

### 11.1 The photochemical pre-step

The VULCAN installation at `EXHALE/VULCAN/` was copied to
`vulcan_work/hd209_vulcan/` and run there, so the shared tree keeps the
configuration and the output of the HD 189733 b comparison. The copy's
`vulcan_cfg.py` differs from that configuration in exactly seven
planet-specific lines:

| key | HD 189733 b (run A) | HD 209458 b |
|---|---|---|
| `atm_file` | `atm/atm_HD189_Kzz.txt` | `atm/atm_HD209_Kzz.txt` |
| `sflux_file` | `atm/stellar_flux/sflux-HD189_Moses11.txt` | `atm/stellar_flux/Gueymard_solar.txt` |
| `out_name` | `HD189.vul` | `HD209.vul` |
| `r_star` | 0.805 | 1.155 |
| `Rp` | `1.138*7.1492E9` | `1.36*6.9911E9` |
| `orbit_radius` | 0.03142 | 0.0480 |
| `gs` | 2140. | 1008.9 |

Everything else is held at the values that produced the converged HD 189733 b
run: `network = thermo/NCHO_photo_network.txt` with `atom_list = ['H','O','C','N']`
(the §"Traps" trap 1 — the atom list matches the network's atoms),
`sl_angle = 48 deg`, `nz = 150`, `P_b = 1e9`, `P_t = 1e-2` dyn/cm^2,
`ode_solver = 'Ros2'`, `use_photo = True`, `use_ion = False`, `diff_esc = []`.

Elemental abundances: `use_solar = True`, so VULCAN reads
`fastchem_vulcan/input/solar_element_abundances.dat` (Lodders 2009) and ignores
the `O_H`/`C_H`/`N_H`/`He_H` lines of the config (trap 2). The values in that
table are C 8.4434, H 12.00, He 10.9864, N 7.9130, O 8.7826 dex, i.e.
C/H = 2.776e-4, He/H = 9.6917e-2, N/H = 8.185e-5, O/H = 6.062e-4.

**Choices, with their reasons.**

- *Stellar spectrum.* No HD 209458-specific UV spectrum ships with VULCAN.
  `Gueymard_solar.txt` (the solar spectrum) was adopted. HD 209458 is a G0
  dwarf, T_eff = 6065 K in `HD209458b/input.inp`; the same file is what
  VULCAN's own `cfg_examples/cfg_HD209.txt` uses; `vulcan_driver.py`'s
  `pick_sflux()` returns exactly this file for 5300 < T_eff < 6800 K; and
  Moses et al. (2011), the source of the T(p) file, modeled HD 209458 b with a
  solar spectrum. This is a substitution, not a measurement of this star.
- *Gravity and radius.* R(1 bar) = 1.36 R_J (the value
  `examples/13_lower_atmosphere/README.md` uses for this planet) and
  M_p = 0.720 M_J (from `HD209458b/input.inp`), with R_J = 6.9911e9 cm — the
  constant `parameters.f90` and `vulcan_to_base.py` both use — give
  g_s = 1008.9 cm/s^2. VULCAN's `Rp` was set to the same 9.508e9 cm so that
  VULCAN's g(z) and the converter's hypsometric integration use one gravity.
  Upstream's `cfg_HD209.txt` instead has Rp = 1.38 x 7.1492e9 cm and
  g_s = 936 cm/s^2; the difference is a radius convention plus a slightly
  different mass.

**Run health.** `MPLBACKEND=Agg OMP_NUM_THREADS=1 python3 -u vulcan.py`
reached steady state in **1206 steps** with a **simulated elapsed time of
1.40e8 s** and 418 s of CPU time — the simulated time is the health indicator
that trap 1 is about, and it advanced normally. The log contains **no**
"Element conservation is violated" line; the reported total atom loss is
H 2.04e-4, O 4.52e-4, C 2.08e-4, N 2.39e-4, and the negative-solution counter
is 0. The only warning is `P_b and P_t assgined in vulcan.cfg are out of range
of the input. Constant extension is used.`, because the supplied T(p) spans
9.998e8 to 2.152e-2 dyn/cm^2 while the grid asks for 1e9 to 1e-2; the extension
is at the two ends of the domain and two decades away from the handoff level.

**One thing worth flagging about this planet.** In VULCAN's shipped Moses+2011
profiles, the 1 microbar handoff level is inside HD 209458 b's thermosphere:
T(1 dyn/cm^2) = 2332 K there, against 863 K in the HD 189733 b file. The
handoff therefore starts EXHALE from a base that the lower-atmosphere model
already heated, and EXHALE then computes its own heating on top. This is a
property of the input T(p) and of the choice of handoff level, not of Route 1,
but it is the reason the base state below is so much hotter than the
`Equilibrium temperature: 1450.0` that `HD209458b/input.inp` otherwise uses.

### 11.2 The handoff file

`python3 src/utils/vulcan_to_base.py vulcan_work/hd209_vulcan/output/HD209.vul
vulcan_work/hd209_wind_response/photo --mp 0.720 --r1bar 1.36` (defaults
`--pbase 1e-6`, `--kzz 1e9`) wrote:

```
T_base    2332.25
r_base    1.49400
HeH_base  0.096952
Kzz_base  1.000e+09
q_H2_base 4.236598e-01
p_base    1.000e-06
```

with `q_H2 = 4.2366e-01, q_H = 4.4987e-01, q_He = 0.1258, mu = 1.827` and
`q_H2O = 2.207e-04, q_CO = 3.600e-04, q_CH4 = 3.11e-15, q_CO2 = 2.54e-08,
q_NH3 = 8.20e-11, q_HCN = 2.10e-11` as comments.

**Photochemical against equilibrium at the same (p, T).** Evaluating
`q_h2_equilibrium` (`lower_column.f90:57-71`) at (1 microbar, 2332.25 K) in
Python gives:

| | q_H2 | `h2_bound_fraction` | `ntot_bc` (C/N/O in the EOS budget) |
|---|---|---|---|
| chemical-equilibrium fit | 0.006293 | 0.006254 | 0.9944995 |
| VULCAN / NCHO | 0.423660 | 0.297585 | 0.7031685 |

The photochemical value is **67x larger** than the fit, and `ntot_bc` moves by
**-0.2913 per nucleus, -29.3%**. Both the size and the **sign are opposite to
HD 189733 b** (§2: +0.074, +7%). The reading that fits the two cases is that
the sign follows the base temperature: at 863 K equilibrium keeps more H2 than
the kinetics do, while at 2332 K equilibrium dissociates H2 almost completely
and vertical mixing from below keeps far more of it than equilibrium allows.
That is an interpretation of two data points, not a demonstration.

### 11.3 The A/B setup

Two directories under `vulcan_work/hd209_wind_response/`, each holding a copy
of `HD209458b/input.inp` with `Molecular base: True` appended, a copy of
`HD209458b/metals.inp`, and the `base.inp` above. The two `base.inp` files
differ in one line: `eqfit/base.inp` has `q_H2_base` removed, `p_base` kept in
both so the fit is evaluated at the level the handoff describes. Same
executable, `OMP_NUM_THREADS=16` for both. The startup echoes confirm the
intended path:

```
eqfit:  (input_read) Molecular base: q_H2(base, chem.eq. fit)  = 0.006 -> ntot_bc = 0.994
photo:  (input_read) Molecular base: q_H2(base, photochemical) = 0.424 -> ntot_bc = 0.703
```

`Base BC` is the default density-anchored mode (`n0` fixed by
`Log10 lower boundary number density: 14.00`), so `ntot_bc` enters as the
ghost-cell pressure at fixed `rho_bc`, not through `n0`.

### 11.4 What the runs did

Both runs spent their first ~90000 steps in a damped base oscillation — `du`
fell from ~1e4 to ~0.5 by step 45000-55000, rose again to a few tens near
step 65000, and then settled — before switching PLM -> WENO3 (at step 47792
for eqfit, 53399 for photo) and entering a monotone decay. The two traces have
the same shape, so the oscillation is a property of this hot base and not of
the H2 partition.

Neither run reached the JFNK hand-off by marching. `du` decayed as a slow power
law of the step count (local exponent between about 0.9 and 1.4 over the last
several 10^5 steps), so neither the
`du < 1e-2` hand-off nor the du-plateau criterion (relative `du` change < 1e-6
held for 2000 steps) was met before the `count_max = 1e6` hard cap. Both runs
therefore stopped at the cap, symmetrically:

| | steps | final `du` | log10 Mdot [g/s] | mean rho v r^2 over r >= 2 R_p [g/s/sr] | window spread |
|---|---|---|---|---|---|
| eqfit | 1000000 | 1.3819e-2 | **10.04** | 1.7280e+09 | 1.37e-2 |
| photo | 1000000 | 1.1620e-2 | **9.98** | 1.5287e+09 | 1.16e-2 |

The window is the code's own convergence range, r >= `Escape radius` = 2 R_p
(95 of 504 cells); the spread is (max - min)/mean of `rho*v*r^2` there, the
same quantity the `[diag]` line prints.

**Measured difference: -0.060 dex in Mdot (-12.9%), -11.5% in the mass flux
over the du window, against a residual non-flatness of 1.2-1.4% in that same
window.** The effect is about eight times the spread, which is the comparison
§7.2 asks for.

Both runs were then restarted from their own end states (`Load IC? True`,
`Solver: Newton 5.0e-2`) to force the JFNK finish. It engaged at step 2 in both
and **failed in both** with `info=2` ("no new best residual for 15 iterations"),
returning its starting iterate (||R|| = 6.467e-3 for eqfit, 5.292e-3 for
photo). The code then re-armed du-based marching, with secondary ionization now
switched on, and both runs stopped on the du plateau at a *worse* `du` than
they started from:

| | steps | final `du` | log10 Mdot [g/s] | mean rho v r^2 [g/s/sr] | window spread |
|---|---|---|---|---|---|
| eqfit | 54915 | 1.0744e-1 | 9.98 | 1.4901e+09 | 1.01e-1 |
| photo | 47437 | 1.0082e-1 | 9.93 | 1.3254e+09 | 9.54e-2 |

Same sign and nearly the same size (-0.05 dex, -11.1% in the window flux), but
here the spread is 10% and the difference is not separable from it. This second
phase is reported as a consistency check, not as the measurement.

**The caveat this leaves.** The standing rule is that a `du`-threshold stop
without the Newton finish is not quantitative for Mdot, and neither of these
runs was Newton-finished. What can be said is narrower and, on this evidence,
still useful: at a matched, well-converged marching state (`du` ~ 1.2-1.4e-2,
window spread ~1.3%) the two runs differ by 11-13% in mass loss, in the
direction the base particle count predicts, and the same difference survives a
change of numerical path (the restart, with different physics active).

### 11.5 Where the difference lives, and the sign

The base. With the density-anchored base BC the ghost pressure is
`W_in(3) = ntot_bc + dp_bc` at fixed `rho_bc`, so the base temperature should
scale with `ntot_bc`, and it does, to six digits:

| | ghost-cell T [K] | ratio |
|---|---|---|
| eqfit | 2317.68 | |
| photo | 1638.73 | 0.707057 |
| `ntot_bc` ratio 0.7031685/0.9944995 | | 0.707058 |

The photochemical solution keeps far more H2 at 1 microbar than equilibrium
does at 2332 K, so more H nuclei are bound into molecules, `ntot_bc` is lower,
the base pressure and hence the base temperature are lower, and the wind is
weaker. The difference is confined to the innermost cells: at r = 1.0035 R_p
the two runs are at 1968 K and 1590 K, and by the r >= 2 R_p window they differ
only through the flux level. **The measured direction matches the expected
one** (a larger H2 subtraction lowers `ntot_bc` and weakens the wind); note
that §7.2 phrased the expectation for HD 189733 b, where the photochemical q_H2
is the *smaller* of the two and the sign is therefore reversed.

### 11.6 Two things found on the way

- **`vulcan_driver.py` cannot configure the in-tree VULCAN as it stands.** Its
  `subs` table matches literal strings in `VULCAN/vulcan_cfg.py`, and one of the
  seven, `out_name =  'HD189-photo.vul'`, no longer exists: the file was edited
  for the HD 189733 b comparison and now reads `out_name = 'HD189.vul'`. The
  driver exits with `ERROR: cfg anchor missing`. The other six anchors still
  match. This is why the run above was configured by hand rather than through
  the driver. Not fixed here.
  [2026-08-15: fixed since. `src/utils/vulcan_driver.py` matches anchored
  regular expressions (`^out_name\s*=.*$`, ...) instead of literal config
  lines, so an edited `vulcan_cfg.py` no longer breaks it.]
- **`vulcan_driver.py` mixes two Jupiter radii.** It computes
  `gs = G Mp / (r1bar * RJ)^2` with `RJ = 6.9911e9` imported from `run_lower`,
  but substitutes `Rp = <r1bar>*7.1492E9` into the config. VULCAN uses `Rp`
  for g(z) = g_s (Rp/(Rp+z))^2, so the two disagree by 4.5% in radius and 9% in
  the implied mass. Not fixed here.
  [2026-08-15: fixed since. `vulcan_driver.py` substitutes `Rp` with the same
  `RJ = 6.9911e9` it uses for `gs`.]

### 11.7 Deeper handoff at p_base = 1e-4 bar (2026-08-10)

§11.1 flagged that the 1 microbar handoff sits inside HD 209458 b's
thermosphere. This section separates the two things §11 measured together: the
choice of *handoff level* and the choice of *q_H2 source* at that level.
New run directories: `vulcan_work/hd209_wind_response/{eqfit_deep,photo_deep}/`.
No VULCAN rerun -- every handoff below comes from the same converged
`vulcan_work/hd209_vulcan/output/HD209.vul`.

**The input T(p).** `VULCAN/atm/atm_HD209_Kzz.txt` (Moses et al. 2011) spans
2.152e-8 to 999.8 bar in 183 levels. Its temperature minimum is **1266 K at
p = 0.0198 bar**; above that level T rises monotonically to 5891 K at the top,
so there is **no isothermal region and no second minimum** anywhere in the
handoff range. Interpolated in log p from the VULCAN output:

| p [bar] | T [K] | n_tot = p/kT [cm^-3] | r(p) [R_J] |
|---|---|---|---|
| 1e-3 | 1538.7 | 4.71e15 | 1.4098 |
| 1e-4 | 1850.6 | 3.91e14 | 1.4311 |
| 1e-5 | 2146.1 | 3.38e13 | 1.4583 |
| 1e-6 | 2332.3 | 3.11e12 | 1.4940 |

r(p) is `vulcan_to_base.py`'s hypsometric integration with M_p = 0.720 M_J and
r(1 bar) = 1.36 R_J.

**The level the stated criterion picks is 1e-3 bar, and it does not run.**
1e-3 bar is the deepest of the four levels considered (the profile itself runs
on down to 999.8 bar), and the only one whose temperature is near
the planet's equilibrium temperature (1538.7 K against the
`Equilibrium temperature: 1450.0` of `HD209458b/input.inp`, +6%), and its
radius, 1.4098 R_J, reproduces that file's `Planet radius: 1.401 R_J` to 0.6%.
Reading those two agreements together, the production configuration of this
planet looks like an implicit 1e-3 bar base. Both members of the pair were
built there and both terminated with a NaN at r = 1.067 after 1846 steps, with
negative HII/HeI densities and the time step collapsed to dtu = 1.2e-6. A
repeat with `IC mode: auto` did not crash but did not settle: at the
30000-step cap it stood at du = 3.4e7 and reported `Mdot = NaN`. Three
shallower levels run normally from the same cold isothermal IC (30000 steps,
photochemical handoff): du = 10.77 at 1e-4 bar, 11.57 at 3e-5, 7.50 at 1e-5,
all comparable to the shallow pair at the same step count (5 to 53).

**Why the deep level fails, in the code's own terms.** The base BC anchors the
mass density at `rho_bc` and sets the ghost pressure to `ntot_bc`, so the base
temperature the model runs with is

```
T(ghost) = T_base * ntot_bc / [(1 + HeH + sum X_metal)/(1 + HeH)]
```

which reproduces all four measured ghost temperatures exactly (2317.68,
1638.73, 838.80, 837.90 K) and generalizes the ratio §11.5 recorded. The
molecular-base correction removes the H nuclei bound into H2 from the particle
count but not from the mass, so a strongly molecular base is handed to EXHALE
as a *cold* base:

| p_base [bar] | T_base [K] | q_H2 (VULCAN) | ntot_bc | T(ghost) [K] | lambda = GM m_H/(k T r) |
|---|---|---|---|---|---|
| 1e-3 | 1538.7 | 0.8356 | 0.5455 | 838.8 | 133.5 |
| 1e-4 | 1850.6 | 0.8020 | 0.5557 | 1027.6 | 107.4 |
| 1e-5 | 2146.1 | 0.6563 | 0.6045 | 1297.3 | 83.5 |
| 1e-6 | 2332.3 | 0.4237 | 0.7032 | 1640.0 | 64.5 |
| production, no `base.inp` | 1450.0 | - | 1.0000 | 1450.0 | 77.8 |

At 1e-3 bar the escape parameter is 133, nearly twice the production value: the
base is far more tightly bound than anything this configuration is used to
integrate, and the cold isothermal IC is correspondingly far off. This is a
property of representing a molecular base inside an atomic-chemistry code, not
of the handoff. With rho anchored, matching the true base *pressure* -- what
the code does, and the right choice for the hydrodynamics, since it gives the
correct pressure scale height and sound speed -- necessarily gives the wrong
base *temperature*, because the model's mean molecular weight is that of an
atomic gas. `vulcan_to_base.py` writes the matching warning into every
`base.inp` it produces ("molecular base persists photochemically -- Tier-2 ...
is required"). Tier-2 is the consistent resolution, and it was not reachable
here until 2026-08-13: the parser rejected `Molecular chemistry` together with
a `metals.inp`, and the A/B carries metals. That restriction is gone --- the
metal stages are now solved in the same system as the molecular network
(`System_HeH_mol_metals`) --- so the deep handoff can be run in the consistent
configuration; see the smoke test recorded in `Update_EXHALE.md`.

**The pair therefore runs at p_base = 1e-4 bar**, the deepest level that
integrates with the shallow pair's numerics unchanged (cold isothermal IC, same
`input.inp` with `Molecular base: True`, same `metals.inp`, same executable,
`OMP_NUM_THREADS=16`). `vulcan_to_base.py <vul> <dir> --mp 0.720 --r1bar 1.36
--pbase 1e-4 --kzz 1e9` wrote:

```
T_base    1850.59
r_base    1.43107
HeH_base  0.096963
Kzz_base  1.000e+09
q_H2_base 8.020378e-01
p_base    1.000e-04
```

with `q_H2 = 8.0204e-01, q_H = 3.7860e-02, q_He = 0.1592, mu = 2.314` and
`q_H2O = 3.695e-04, q_CO = 4.559e-04, q_CH4 = 9.806e-12, q_CO2 = 5.230e-08,
q_NH3 = 3.161e-09, q_HCN = 2.786e-08` as comments. `eqfit_deep/base.inp` is the
same file with the `q_H2_base` line removed. Startup echoes:

```
eqfit_deep:  T0 -> 1850.6 K,  R0 -> 1.4311 R_J,  He/H -> 0.09696,  p_base -> 1.00E-04 bar
             Molecular base: q_H2(base, chem.eq. fit)  = 0.794 -> ntot_bc = 0.558
photo_deep:  (same)  q_H2(base) -> 0.80204
             Molecular base: q_H2(base, photochemical) = 0.802 -> ntot_bc = 0.556
```

**The photochemical-vs-equilibrium gap closes with depth.** Evaluating the
Visscher/Koskinen fit (`lower_column.f90:57-71`) at each level's own (p, T),
with C/N/O in the EOS budget and HeH = 0.09696:

| p_base [bar] | q_H2 photochemical | q_H2 equilibrium | ratio | ntot_bc photo | ntot_bc eq | difference |
|---|---|---|---|---|---|---|
| 1e-6 | 0.42366 | 0.00629 | 67x | 0.70317 | 0.99450 | -29.3% |
| 1e-5 | 0.65628 | 0.24392 | 2.7x | 0.60452 | 0.80467 | -24.9% |
| 1e-4 | 0.80204 | 0.79376 | 1.010x | 0.55568 | 0.55824 | **-0.46%** |
| 1e-3 | 0.83555 | 0.83767 | 0.997x | 0.54555 | 0.54495 | +0.11% |

Two readings, both tentative. First, the 67x disagreement §11.2 reported is not
a general property of photochemistry against equilibrium for this planet; it is
what the 1 microbar level does. By 1e-4 bar the two agree to 1%, and at 1e-3
bar the sign reverses. Second, the reason appears to be temperature: the
equilibrium fit dissociates H2 steeply above ~2000 K while vertical mixing
holds the photochemical q_H2 near its deep value, so the two curves separate
exactly where the thermosphere begins. If that is right, the *level* and the
*q_H2 source* are not independent choices -- the q_H2 effect exists only at
levels that are already thermospheric.

**What the runs are compared against.** There is no Newton-finished HD 209458 b
reference in this working copy to anchor an absolute Mdot on. Measuring the
stored profiles directly (4 pi rho v r^2 twenty cells from the top, halved for
`Rate/2 + Mdot/2`, R_p = 1.401 R_J):

> **2026-08-11.** That statement no longer holds, and the table row for
> `HD209458b/output` below is stale. `HD209458b/` has since been re-converged
> from a restored solar-composition `metals.inp` with the corrected JFNK line
> search and the metal-line trapping / coronal-cutoff cooling, reaching
> `info = 0`, `||R|| = 7.690e-04`, `log10 Mdot = 9.31`. That is the current
> production anchor for this planet. The absolute gaps quoted in the readings
> below (0.21 dex against 9.49, 0.49-0.54 dex against the 1 microbar handoff)
> were measured against the old stored run and are not restated here, because
> the handoff runs themselves predate the same two changes; re-anchoring them
> means re-running the handoff pair, not re-arithmetic. What the section is
> being cited for -- that the wind response is a property of the handoff
> *level* rather than of using photochemistry, measured as a difference between
> two runs sharing everything else -- is a differential result and is not
> affected.

| stored run | log10 Mdot | T(ghost) | window spread | note |
|---|---|---|---|---|
| `HD209458b/output` | 9.61 | 1450.0 | 1.4e-1 | its `EXHALE.out` reports L_EUV 27.83 / L_X 26.39 and `Reconstruction method: PLM`, against 27.93 / 27.20 and PLM+WENO3 in the current `input.inp` |
| `HD209458b/output_metals_on` | 9.49 | 1445.5 | 4.5e-2 | same vintage; the metals-on comparison |
| `HD209458b/output_metals_off` | 10.43 | 1450.3 | 3.4e-2 | metals-off |
| `benchmarks/hd209` | 10.59 | 2333.9 | 1.9e-2 | has its own 1 microbar `base.inp` (T_base 2333.94, r_base 1.50345, no `q_H2_base`) and `He_diffusion: True` |

The "production ~9.3" figure this study set out to compare against is not
reproducible here; the nearest metals-on number is **9.49**, from a run whose
XUV input and reconstruction differ from the current configuration and whose
own window spread is 4.5%.

**Four-run comparison, matched marching states.** Same measurement as §11.4:
`log10 Mdot` from the code's own report, mean `rho v r^2` and its
(max - min)/mean over the convergence window r >= 2 R_p.

| run | p_base | steps | final du | log10 Mdot | mean rho v r^2 [g/s/sr] | window spread | ntot_bc | T_base [K] | T(ghost) [K] | r_base [R_J] |
|---|---|---|---|---|---|---|---|---|---|---|
| eqfit | 1e-6 | 1000000 | 1.382e-2 | **10.037** | 1.7280e+09 | 1.37e-2 | 0.99450 | 2332.2 | 2317.7 | 1.4940 |
| photo | 1e-6 | 1000000 | 1.162e-2 | **9.984** | 1.5287e+09 | 1.15e-2 | 0.70317 | 2332.2 | 1638.7 | 1.4940 |
| eqfit_deep | 1e-4 | 300000 | 4.233e-2 | **9.857** | 1.1295e+09 | 4.14e-2 | 0.55824 | 1850.6 | 1032.3 | 1.4311 |
| photo_deep | 1e-4 | 300000 | 4.234e-2 | **9.856** | 1.1287e+09 | 4.14e-2 | 0.55568 | 1850.6 | 1027.6 | 1.4311 |

**The same four, Newton-finished** (`Solver: Newton 5.0e-2` restarts from each
run's own end state, with `Secondary_ionization: False` -- see §11.8 for why
that key is needed):

| run | JFNK | residual reached | log10 Mdot | mean rho v r^2 | window spread |
|---|---|---|---|---|---|
| eqfit 1e-6 | info=2 | 2.041e-3 | 10.029 | 1.6965e+09 | 7.80e-3 |
| photo 1e-6 | info=2 | 1.348e-3 | 9.978 | 1.5108e+09 | 6.86e-3 |
| eqfit 1e-4 | **info=0** | 9.278e-4 | 9.705 | 8.0771e+08 | 5.96e-3 |
| photo 1e-4 | **info=0** | 5.716e-4 | 9.703 | 8.0368e+08 | 6.65e-3 |

**The measured separation.**

| pair | d(log10 Mdot) | d(window flux) | window spreads |
|---|---|---|---|
| 1e-6, marching | -0.054 | -11.5% | 1.4% / 1.2% |
| 1e-6, Newton-finished | -0.050 | -10.9% | 0.78% / 0.69% |
| 1e-4, marching | -0.000 | -0.1% | 4.1% / 4.1% |
| 1e-4, Newton-finished | -0.002 | -0.5% | 0.60% / 0.67% |

**Readings.**

- The wind response §11 measured is a **property of the handoff level, not of
  using photochemistry.** At 1e-4 bar the photochemical and equilibrium
  handoffs differ by -0.5% in window flux against a 0.6-0.7% residual
  non-flatness -- inside the noise -- where the 1 microbar pair differs by
  -11% against 0.7-0.8%. The -0.46% difference in `ntot_bc` predicts exactly
  this, and it is the same sign as at 1 microbar, just 60x smaller.
- The deeper handoff **does** move the absolute mass loss toward the stored
  metals-on number: 9.70 against 9.49, a 0.21 dex gap, where the 1 microbar
  handoff gives 9.98-10.03, a 0.49-0.54 dex gap. Some of the 1 microbar excess
  is therefore attributable to the level -- but not, on this evidence, to base
  *temperature* alone: the 1 microbar photochemical case already runs a 1640 K
  base against the equilibrium case's 2318 K, and the two differ by only 0.05
  dex. The base radius moves with the level too (1.4940 -> 1.4311 R_J), and the
  escape parameter is not monotone across the set (64.5, 107.4, 77.8 for
  1e-6 photo, 1e-4, production).
- The base region stays out of steady state in **every** state measured,
  including the two that reached ||R|| < 1e-3: the innermost live cell carries
  an inward `rho v r^2` of 1.2e4 (shallow eqfit, marching) to 5.7e4 (deep
  photo, info=0) times the wind value, and 188 to 226 of the ~268 cells below
  r = 1.2 carry inward flux in four of the five states examined (111 of 268 in
  the fifth, shallow photo Newton-finished). A converged volume-weighted
  residual norm does not imply a flux-flat base.

**Recommendation.** The deeper handoff **does not supersede** the shallow one as
a recommended configuration, and neither is recommended for quantitative work
on this planet as things stand. What the evidence supports is narrower:

1. For the specific question "does the photochemical q_H2 change the wind?",
   the answer depends entirely on the handoff level, and the 1 microbar answer
   (-11%) should be quoted as such, not as a general photochemistry-vs-
   equilibrium effect.
2. Any level deep enough for T(p_base) to approach T_eq drives the
   molecular-base EOS correction outside its usable range (base temperature
   depressed to 0.55 x T_base, lambda = 133, NaN at 1e-3 bar).
3. The consistent configuration is Tier-2 molecular chemistry at a deep
   handoff, which the parser forbade together with metals until 2026-08-13.
   That restriction has been lifted, so the quantitative deep handoff on a
   molecular-base planet is now a matter of running it; the pair above was
   built before the lift and still carries the 1e-4 bar level.

---

### 11.8 Newton-finish diagnosis (2026-08-10)

The shallow pair never produced a Newton-grade Mdot. This section says what
blocks it, from the existing logs plus four short restarts. It is a diagnosis;
no fix is attempted.

The four restarts are kept as
`vulcan_work/hd209_wind_response/{eqfit,photo,eqfit_deep,photo_deep}_newton/`.
Each is its parent run with `Load IC? True`, `Solver: Newton 5.0e-2` and
`Secondary_ionization: False`; its Load-IC input was a copy of the parent's
final `Hydro_ioniz.txt` / `Ion_species.txt`, which is why those two files are
not kept a second time inside the restart directories.

**What the marching would have needed.** Sampling `du` every 2000 steps out of
`run_march.log`: after step ~2e5 both traces are strictly monotone (0 sign
changes in 399 samples) and follow a power law, local exponent
d ln du / d ln count = -1.05 (eqfit) and -1.27 (photo) over the last half
decade. Before step ~2e5 they are violent, `du` reaching 3.3e4 at step 18000.
Extrapolating the power law from the 1e6-step end states (1.38e-2, 1.16e-2):
the `du < 1e-2` JFNK hand-off would have fired at about 1.36e6 (eqfit) and
1.13e6 (photo) steps -- 2.0 h and 43 min of further marching -- but the
`du < 1e-3` marching target needs about 1.2e7 and 6.9e6 steps, i.e. 68 h and
38 h at the measured 5.5 h per 1e6 steps. Marching alone is not a route to the
target on this planet.

**The base-cell momentum residual is where JFNK always ends up, everywhere.**
Every `(JFNK)` iteration line in both restart logs reports the worst relative
residual at `k=2` (momentum) in cells `j = 1-9`, `r = 1.000-1.002`. That is not
specific to this configuration: across the eight JFNK runs archived in
`backup/regression/` -- **none of which has a `base.inp`** -- the last
iteration reports `k=2` in all eight, at `j=1, r=1.000` in six, `j=2` in one
and `j=21, r=1.004` in one, and it does so in the three that *converged*
(`jfnk_hd189` info=0 at 9.46e-4, `ptc_warm/tight` info=0 at 6.04e-5,
`solver_newton_cold` info=0 at 5.42e-4) as well as in the five that did not.
The base momentum is the terminal residual of every EXHALE steady solve;
whether JFNK reports info=0 or info=2 is whether it gets that residual under
the target. This is `TO_BE_DONE.md` item (A), quantified there on WASP-121 b at
the same signature.

**The profiles carry the matching signature.** In the last snapshot of both
1e6-step marches the innermost live cell is a strong inflow spike against a
zero-velocity ghost -- v = -701 cm/s (eqfit) and -837 cm/s (photo) at
r = 1.00019, where `rho v r^2` is inward and 1.2e4 to 1.7e4 times the wind
value -- and `dv/dr` changes sign 8 to 10 times inside r < 1.2 while changing
sign nowhere above. The base region has also not settled: between the marching
end state and the restart end state 55000 steps later the same cell moves from
-701 to -718 cm/s and from 1839 K to 1731 K.

**The staged secondary-ionization flip fires at the Newton hand-off.** Neither
1e6-step march ever activated secondary ionization -- `run_march.log` contains
no `secondary ionization activated` line, because the staged trigger
(`EXHALE_main.f90:801`) needs a convergence condition the march never met. The
restarts then activated it *at step 2, in the same block that starts the Newton
finish* (`EXHALE_main.f90:838`). JFNK was therefore asked to solve a system the
marching state had never been converged against, which also explains the
post-failure behavior: marching resumed at du = 1.4e-2 and ran *away* to
du = 4e-1 over the next 8000 steps before crawling back.

**A 2x2 separates the two effects.** Four 3000-step restarts, each from its own
run's end state, `Solver: Newton 5.0e-2`, differing only in the base level and
in whether `Secondary_ionization: False` is present:

| base | secondary ionization at the hand-off | eqfit | photo |
|---|---|---|---|
| 1e-6 bar | ON (default, staged) | info=2, no improvement from 6.467e-3 | info=2, no improvement from 5.292e-3 |
| 1e-6 bar | OFF | info=2, 6.467e-3 -> **2.041e-3** | info=2, 5.292e-3 -> **1.348e-3** |
| 1e-4 bar | ON (default, staged) | info=2, 1.770e-2 -> 1.714e-2 | info=2, 1.770e-2 -> 1.714e-2 |
| 1e-4 bar | OFF | **info=0, 9.278e-4** | **info=0, 5.716e-4** |

The two 1e-4 bar entries agree to every printed digit because the deep pair's
`ntot_bc` values differ by 0.46%. The marching that follows the failed solves
reaches du = 1.07e-1 / 1.01e-1 (1e-6, flip on) and 4.42e-2 / 4.42e-2 (1e-4,
flip on) at step 3000, against 7.84e-3 / 6.89e-3 (1e-6, flip off); the two
flip-off deep runs stopped at step 2 because JFNK had already met its target.

**Causal statement.**

1. *Established by the 2x2.* The staged secondary-ionization activation firing
   simultaneously with the Newton hand-off blocks the solve in **both**
   configurations: with it on, JFNK improves ||R|| by at most 3% and the
   subsequent marching degrades; with it off, JFNK improves ||R|| by 3-4x at
   the 1 microbar base and converges outright at the 1e-4 bar base. It is a
   necessary condition for the failure reported in §11.4.
2. *Established, given (1).* The residual floor that remains depends on the
   base state. The thermospheric 1 microbar base floors at 1.3-2.0e-3, above
   the 1e-3 target; the 1e-4 bar base reaches 5.7-9.3e-4, below it. So the hot
   base does contribute, but only as a factor ~2 in the achieved residual, and
   only once the flip is out of the way.
3. *Not implicated.* The molecular-base / `base.inp` code path itself. The same
   terminal signature appears in the archived runs that have no
   lower-atmosphere handoff at all.
4. *Unchanged by any of this.* The base region is not flux-flat even in the
   info=0 states (§11.7, last reading). Item (A) of `TO_BE_DONE.md` -- explicit
   viscosity -- remains the standing fix.
   [2026-08-15: the "standing fix" attribution is withdrawn. §11.10 traced the
   residual floor to the JFNK diagonal scaling and the stagnation watchdog, not
   to a missing viscous term; the fix that made these configurations converge
   was the solver scaling plus the line-search repair of 2026-08-11. Explicit
   viscosity and conduction do exist
   (`src/modules/time_step/viscous_conduction.f90`, keys `Viscosity:` /
   `Conduction:`) but they are not what closed this.]

**Is the deep configuration usable for quantitative work?** Numerically, yes:
it is the only one of the four that reaches ||R|| < 1e-3, and its A/B pair sits
at a 0.6-0.7% window spread. Physically it is still compromised, for the reason
§11.7 gives -- its 1028 K model base is an artifact of applying an EOS-only H2
correction to a base that is 80% molecular. The numbers below 1.2 R_p from that
pair should not be used.

**What this does to the §11 numbers.** They hold, and they are better
constrained than §11 could claim. At the Newton-improved states the 1 microbar
pair differs by **-0.050 dex, -10.9% in the window flux, against a 0.7-0.8%
residual non-flatness** -- about 14x the spread, where the marching pair gave
8x -- and pushing convergence that far moved each run's own Mdot by less than
0.01 dex (10.037 -> 10.029, 9.984 -> 9.978). One caveat to carry: because the
marches never fired the staged trigger and the productive restarts disable it,
**every number in §11, §11.7 and §11.8 is a secondary-ionization-OFF number**;
a default-settings production run would include it.

### 11.9 Ordering fix for the Newton hand-off, and what remains (2026-08-10)

Following the §11.8 finding, the hand-off block in `EXHALE_main.f90` was
changed: when the Newton trigger fires while the staged secondary ionization
is still pending, the code now flips the coupling on, re-arms the marching
stops exactly as the staged activation does (with the `N_stall` hold), and
postpones JFNK until the trigger next fires -- so Newton is only ever asked to
solve the system the marched state was relaxed against. The three default
regression cases never reach the hand-off and stayed byte-identical
(`run_check.sh check` PASS after the change).

Two tests on the deep pair, both with default (ON) secondary ionization and
the new binary, `EXHALE_MAXSTEPS=60000`:

| start state | flip | JFNK engaged | result | afterwards |
|---|---|---|---|---|
| `photo_deep` march end (du 4.2e-2) | step 2 | step 2002 | info=2, \|\|R\|\| stuck at 1.84e-2 | du climbs to 0.146 by step 40719, Mdot 9.80 |
| `photo_deep_newton` converged state (du 1.9e-3, \|\|R\|\| 5.7e-4) | step 2 | step 2002 | info=2, \|\|R\|\| 2.80e-3 | du plateau 0.103 at step 41375, Mdot 9.66 |

Kept as `vulcan_work/hd209_wind_response/photo_deep_newton_fix/` and
`.../photo_deep_secion_cont/`.

Reading, stated tentatively: the sequencing now behaves as designed, and the
second test shows that even the best available starting state does not survive
2000 steps under the activated coupling -- the state drifts to a residual
floor of a few 1e-3 (base momentum rows again) and long marching settles on
the familiar du ~ 0.10-0.15 base oscillation. So the ordering was necessary
but is not sufficient here: with secondary ionization active, this
configuration appears to have no steady state reachable by either marching or
JFNK, which is the `TO_BE_DONE` item (A) base-momentum problem amplified by
the coupling, not a hand-off artifact. The secondary-ionization-OFF caveat at
the end of §11.8 therefore still stands for every quantitative number in §11;
lifting it waits on item (A), for which the recorded direction is the explicit
base viscosity work.

> [2026-08-15: the last clause is superseded by §11.10 below — the floor was the
> JFNK scaling and watchdog, not a missing base viscosity.]

### 11.10 The §11.9 residual floor was the solver, not the base (2026-08-10)

The reading at the end of §11.9 — that `photo_deep_secion_cont` has no steady
state reachable by marching or JFNK, and that the `2.80e-3` floor is the
`TO_BE_DONE` item (A) base-momentum problem — is **withdrawn**. It was traced
instead to the JFNK diagonal scaling and to the stagnation watchdog; full
account and numbers in `docs/newton_scaling_and_base_wall.md`.

Two things were measured. First, the base momentum row is solvable: `R(2,1)` is
the remainder of a four-order cancellation, `6.5e-5` of the gravity term, and a
2.9 ppm change of the ghost pressure nulls it. Second, the residual the merit
was actually dominated by lived at `r = 1.32-1.37`, where `rho v` changes sign
at the hand-off, and it appeared there only because the momentum scale was
floored at `1e-6` of the *base* `|rho v|`. The watchdog then aborted after 15
iterations without a new best `||R||`, while every variant of this solve needs
27-38 iterations before it first improves on the warm start.

With the local scale `D_mom = rho(|v| + c_s)` and a watchdog on consecutive
failed line searches, the same configuration converges: `info = 0`,
`||R|| = 5.053e-04`, 59 outer iterations, `log10 Mdot = 9.47`. The `9.66-9.70`
values quoted in §11.8-11.9 came from marching stopped on `du` after the JFNK
failure; `9.47` is the first Newton-converged number for this configuration and
is the one to carry forward. The secondary-ionization-OFF caveat on §11's other
numbers is unaffected — those runs never fired the staged trigger.

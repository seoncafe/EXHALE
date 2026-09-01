# How EXHALE's chemistry solve compares with Photochem and VULCAN, and what to take from them

**Purpose.** This memo compares selected chemistry algorithms in EXHALE,
Photochem, VULCAN, and the Equilibrate library. They do not solve nominally the
same problem: EXHALE imposes local steady reaction balances inside a
hydrodynamic update, Photochem and VULCAN integrate coupled kinetics and
transport in time, and Equilibrate minimizes thermodynamic free energy. The
comparison is therefore about specific numerical and rate-handling techniques,
not about interchangeable solvers.

Measured 2026-08-31, against:

- **VULCAN** — `VULCAN/` (chemical kinetics, 18 reaction-network files under
  `thermo/`).
- **Photochem** — `photochem/` v0.9.0 and the source present under
  `photochem/src/`; its build declares Equilibrate 0.2.6 and CVODE 5.7.0 as
  dependencies.

**Where the Equilibrate source is.** The expanded dependency tree is **not**
under `EXHALE_v1.00/photochem/`; it is under the workspace-level clone, at
`../photochem/_skbuild/linux-x86_64-3.11/cmake-build/_deps/equilibrate-src/`.
A review that looked only inside `EXHALE_v1.00/` reported it missing and
marked the Equilibrate citations of section 5 unverifiable; they were
re-read at that path on 2026-08-31 and every line quoted below is verbatim
from it. All Equilibrate line numbers refer to
`_deps/equilibrate-src/src/equilibrate_cea.f90`.

**A caution about provenance.** The two files under
`photochem/src/dependencies/patches/` are this project's own local changes,
not upstream Photochem, and nothing here cites them as independent evidence.

## Review outcome: corrections required

Implementation review found the following material problems in the original
memo:

| original claim | review finding | required correction |
|---|---|---|
| VULCAN does not carry `H3+` | `H3_p` is present in `thermo/all_compose.txt`, but no checked reaction network uses it | distinguish composition metadata from an active mechanism |
| the stated R15/R12 equilibrium table | the numbers do not reproduce from the stated pressure, He/H, rate functions, and definition of `q_H2` | use the corrected table and show the conversion formula |
| Photochem freezes composition above a quench level | the gas-giant extension applies quenching only to the initial guess, then integrates the selected photochemical domain | do not use this as precedent for freezing EXHALE cells |
| CEA's abundance-weighted update can be added to EXHALE acceptance | an update test is not a reaction-balance test; EXHALE already rechecks every active and held row | use it only as an optional iteration diagnostic |
| CEA trace-species limits are local additions to MINPACK `hybrd` | `hybrd` exposes one trust radius and no component bounds | this requires a different step controller or variable transformation |
| `u_shift = 1` supplies `h = eta max(1,abs(u))` | it supplies `h = eta abs(1+u)` except at exact zero, so the collapse moves to `u = -1` | implement an actual absolute floor in a constrained-specific difference routine |
| EXHALE cannot use time integration | no cost or accuracy comparison establishes that categorical claim | retain local equilibrium as the current design, but treat integration or pseudo-transient continuation as benchmarkable alternatives |
| all reverse rates in both reference codes come from thermodynamics | the code applies thermodynamic reversal only to reactions declared or indexed as reversible | narrow the statement to configured reversible pairs |

The cell-387 failure that motivated the comparison has also been resolved by a
missing He+ + H channel, as recorded in
`docs/charge_exchange_cancellation_limit.md`. Numerical improvements remain
useful, but they are no longer remedies for that failure.

---

## 1. The headline

**Neither code validates the complete EXHALE H/He ion-molecule network. The
transferable results are narrower: thermodynamic reversal for reactions that
are true microscopic reverse pairs, falloff-rate support, independently
verified derivatives, and careful treatment of trace species.**

The rate comparison is incomplete because the ion reactions of interest are
absent from the checked active mechanisms (section 2). Photochem and VULCAN
reach steady states by time integration, whereas EXHALE currently solves local
algebraic balances (section 4). That difference prevents a direct solver swap,
but it does not establish that time integration or pseudo-transient
continuation is unaffordable without a benchmark.

---

## 2. Rate cross-validation: the reactions are not there

### 2.1 He <-> H charge exchange is absent from both

**VULCAN** carries the ion species `He_p` and `HeH_p` in its composition table
(`VULCAN/thermo/all_compose.txt:195-196`) but in no reaction network. Helium is
a spectator, and the networks state it explicitly:

```
653  [ He -> He                           ]  0.00E+00     0.000       0.0
```
(`VULCAN/thermo/NCHO_photo_network.txt:362` — a self-reaction with rate
coefficient zero.) Photoionization cross sections for helium ship with the
distribution (`thermo/photo_cross/He/`), but no network in this tree uses them.

The checked **Photochem** mechanism
(`photochem_clima_data/.../reaction_mechanisms/zahnle_earth.yaml`) declares 98
gas species and 16 particle species, with no charged species. Helium appears
as an element and a gas species and in no reaction. This is a statement about
that mechanism file, not a general limitation of the Photochem parser.

So neither `He+ + H -> He + H+` nor `He + H+ -> He+ + H` exists in the checked
active mechanisms. The non-radiative electron-capture channel EXHALE adopted on
2026-08-31 (Kingdon & Ferland 1996, from Zygelman et al. 1989) has no
independent implementation here to check against, and neither does the question
of whether the radiative and non-radiative channels should be separated. **The
absence is a measured result for the local mechanism files**: search returns
`He_p` in the VULCAN composition catalog but not in its network files, and the
Photochem mechanism has no charged-species declaration.

### 2.2 H3+ chemistry is absent from the checked active mechanisms

VULCAN's composition catalog does carry `H3_p`
(`VULCAN/thermo/all_compose.txt`), just as it carries `He_p` and `HeH_p`, but
none of the checked network files uses it. The `zahnle_earth.yaml` mechanism
contains no charged species. EXHALE's H3+ formation, destruction and cooling
(`src/modules/lower_atmosphere/mol_rates.f90`, `h3p_cooling.f90`) therefore
have no active counterpart in these mechanisms. Presence in a composition
catalog must not be reported as either an implemented reaction network or
complete absence from the code.

### 2.3 The dissociation-equilibrium fit is ours alone

Searching both trees for Visscher 2006, Koskinen 2022, or the numeric
signatures of the fit (`23672`, `1.9845`, `6.2645`) returns nothing. Neither
code computes an H2/H partition from a closed-form fit at all:

- **Photochem's equilibrium-initialization path** minimizes Gibbs free energy
  with the **NASA CEA** algorithm (Gordon & McBride 1994), through the
  Equilibrate library, using NASA9 / NASA7 / Shomate polynomials. Its kinetic
  steady state is obtained by time integration, as described in section 4.
- **VULCAN** calls **FastChem** when `ini_mix = 'EQ'`; it also supports saved,
  tabulated, and constant initial compositions (`VULCAN/build_atm.py`).

EXHALE's `q_h2_equilibrium` (`src/modules/lower_atmosphere/lower_column.f90`)
is a closed-form shortcut with no counterpart in either code. It is not used
only as a continuation seed: callers also use it for lower-column structure,
initial conditions, the base particle count, and the molecular-basin retry.
Replacing it would therefore change several physical initialization and
boundary paths. The cost of a free-energy calculation in those paths has not
been measured here.

---

## 3. Thermodynamic reverse rates and falloff

**Both codes can derive reverse rates from thermodynamics for reactions that
their mechanisms designate as reversible.** This does not apply to photolysis,
one-way kinetic channels, or VULCAN reactions after its `# reverse stops`
marker.

VULCAN pairs each reaction with its reverse and divides by the equilibrium
constant (`VULCAN/op.py:302-304`):

```python
var.k_fun[i] = lambda temp, mm, i=i: var.k_fun[i-1](temp, mm)/chem_funs.Gibbs(i-1,temp)
var.k[i] = var.k[i-1]/chem_funs.Gibbs(i-1,Tco)
```

where `Gibbs(i,T)` is generated automatically from NASA polynomial free
energies for the indexed reverse pairs. Photochem does the same for reactions
parsed as reversible in
`photochem/src/photochem_common.f90:180-183`:

```fortran
Dg_forward = gibbP_forward - gibbR_forward ! DG of the forward reaction (J/mol)
rx_rates(j,i) = rx_rates(j,n) * &
                (1.0_dp/exp(-Dg_forward/(Rgas * var%temperature(j)))) * &
                (k_boltz*var%temperature(j)/1.e6_dp)**(m-l)
```

with the trailing factor converting a molar-concentration equilibrium constant
to number density. The structure enforces thermodynamic consistency for those
configured pairs; it does not validate whether a mechanism author correctly
classified two physical channels as a reversible pair.

**EXHALE fits both directions independently, and for at least one pair that is
a defect.** The H2 three-body recombination and its collisional dissociation
are the *same channel* run forwards and backwards, so detailed balance does
apply to them — unlike the He <-> H pair, where the two published fits describe
a radiative and a collisional channel and are correctly independent
(`docs/molecular_chemistry_audit_he_rich.md`). EXHALE has:

- `H + H + M -> H2 + M` : `8.0e-33 (300/T)^0.6` (Ham et al. 1970),
  `mol_rates.f90`
- `H2 + M -> H + H + M` : `1.5e-9 exp(-48350/T)` (Baulch et al. 1992),
  `mol_rates.f90`

**Fixed 2026-09-01, and the rest of this section is the record of the defect,
not of the current code.** `mol_rates.f90` now carries a single recombination
coefficient `k3b_H_H_to_H2` = `2.8e-31 T^-0.6` (Cohen &amp; Westberg 1983,
recommended over 50-5000 K, the range that covers the molecular layer; Ham
et al. measured only 77-300 K), and R12 is its thermodynamic reverse,
`k3b_H_H_to_H2(T)/keq_H_H_to_H2(T)`, with `K_eq` evaluated from the H2
rovibrational partition function and D0. The two directions are therefore
consistent by construction. Sources, validity ranges and three independent
checks of `K_eq` are recorded at the R12 comment in `mol_rates.f90`; the
tables below were measured with the old pair and are kept as the diagnosis.

Their ratio is an equilibrium constant, and it can be compared directly with
EXHALE's own `q_h2_equilibrium` seed. Evaluated at p = 1e-6 bar with
He/H = 0.0793, treating only H, H2 and He (M cancels in the ratio):

| T [K] | q_H2 from the Visscher/Koskinen fit | q_H2 from the R15/R12 rate ratio |
|---|---|---|
| 1500 | 0.8239 | 0.8359 |
| 1800 | 0.5851 | **0.5040** |
| 1930 | 0.3194 | **0.2297** |
| 2100 | 0.0710 | **0.0441** |
| 2500 | 0.0013 | 0.0009 |

The earlier kinetic column did not reproduce from the stated assumptions. The
corrected values above use `N = p/(k_B T)`, `f = n_He/n_H,nuclei`, atomic H
density `a`, and H2 density `b`:

```text
b = (k_R15/k_R12) a^2
N = (1 + f) a + (1 + 2 f) b
q_H2 = b/N.
```

The two estimates differ through the transition, but this is not by itself a
solver defect. A continuation seed is expected to move toward the kinetic
root. The discrepancy becomes actionable only if a controlled run shows that
it prevents continuation or materially changes the selected physical branch.
The estimate also omits H2+, H3+, HeH+, and the other network channels.

A second, smaller finding sits alongside it. Both reference codes carry a
high-pressure limit for H2 formation; EXHALE carries only the low-pressure
term:

| code | low-pressure limit | high-pressure limit | source |
|---|---|---|---|
| EXHALE | 2.45e-31 T^-0.6 | **none** | Ham et al. 1970 |
| VULCAN | 2.70e-31 T^-0.6 | 3.31e-6 T^-1.0 | `NCHO_photo_network.txt:367` |
| Photochem | 2.696e-31 T^-0.6 | 1.0e-12 | `zahnle_earth.yaml:1180-1184` |

The two reference low-pressure coefficients agree with each other to 0.2% and
EXHALE's is 0.91 of them, a different source but the same magnitude and
exponent. The substantive difference is the missing falloff. The ratio of the
listed high-pressure limits is temperature dependent,
`(3.31e-6 T^-1)/(1e-12) = 3.31e6/T`, not a fixed factor of 3500. Neither value
should be copied without its source, validity range, third-body efficiencies,
and falloff convention. First measure the reduced pressure for EXHALE's base
states using a physically supported falloff model.

---

## 4. Why their overall strategy is not a direct replacement

**Both codes reach a steady state by integrating in time, not by solving for
equilibrium.** Photochem uses SUNDIALS **CVODE** with BDF, a banded linear
solver and a user-supplied Jacobian; its `find_steady_state` is a loop over
robust time steps, converged when the mixing ratios stop changing or a nominal
1e17 s is reached. VULCAN uses a second-order **Rosenbrock** integrator
(Verwer et al. 1997) with a banded direct solve.

Neither checked implementation uses EXHALE's radiation-field continuation.
Time integration provides a trajectory from the initial state, but it should
not be described as a continuation parameter in the numerical-analysis sense
without qualification. A text search for a few method names is also not proof
that no related globalization or pseudo-transient strategy exists.

**EXHALE cannot replace the present solve with a column integration without a
coupling redesign.** Its chemistry is currently a local algebraic closure
inside a hydrodynamic update. A time-dependent chemistry state would require a
choice among operator splitting, a coupled implicit update, or a
pseudo-transient local solver, together with error and conservation controls.
That is a substantial change. It is not proven unaffordable: a reduced local
benchmark on the difficult cells would be needed to decide cost and robustness.

### 4.1 The derivative strategy, which does transfer as a direction

The available VULCAN and Photochem default paths do not use finite differences
for their chemistry Jacobians. The earlier Equilibrate inspection reported an
analytic CEA matrix, but that dependency source is unavailable for this review:

| implementation | derivative |
|---|---|
| Equilibrate (equilibrium) | analytic — the CEA matrix is the exact linearization |
| VULCAN | symbolic, generated by sympy (`VULCAN/make_chem_funs.py:668`) |
| Photochem (default) | forward-mode automatic differentiation, `autodiff = .true.` |

Photochem retains a finite-difference path as an opt-out, and it carries
**exactly the pathology EXHALE measured**: a purely relative step with no
absolute floor,

```fortran
R(j) = var%epsj*abs(wrk%usol(i,j))
```

so the perturbation collapses as the unknown approaches zero. Photochem's
answer was not to fix the step rule but to stop differencing — automatic
differentiation is the default, and the documentation notes that the tolerance
must be tightened to suit it.

EXHALE currently sets `v = 1 + ln(n/n_ref)` and passes
`epsfcn = fd_eta^2`, with `fd_eta = 1e-4`. This improved the saved cell's
derivative comparison, but it does **not** implement the claimed absolute-floor
rule. MINPACK forms

```text
h = fd_eta abs(v) = fd_eta abs(1 + ln(n/n_ref))
```

except when `v` is exactly zero. The perturbation therefore still collapses
when `n/n_ref` is near `exp(-1)`; the shift moves the singular location rather
than removing it. The correct local repair is a constrained-specific
difference routine using `h = eta max(1,abs(log(n/n_ref)))`, or an equivalent
absolute floor. A complete analytic or automatically differentiated Jacobian
remains the durable option. The tree carries `hybrd` but not `hybrj`, so that
path requires a new solver interface or adaptation of the existing Newton
driver.

---

## 5. Candidate techniques and their limits

The Equilibrate snippets in this section were quoted from an expanded build
tree that is no longer present, so their surrounding logic could not be
rechecked. More importantly, a free-energy minimizer's update rules do not
automatically define valid acceptance rules for a kinetic steady-state solve.

### 5.1 Weight the convergence test by abundance

Equilibrate does not test the logarithmic update itself. It tests the update
weighted by the species' share of the total
(`equilibrate_cea.f90:1760`):

```fortran
IF (n_spec(i_reac)*ABS(delta_n_gas(i_reac))/SUM(n_spec) > 0.5d-5) THEN
   gas_good = .FALSE.
```

A trace species whose logarithmic update is large may reasonably delay an
*iteration-step* convergence test less than an abundant species. This rule
must not be added to EXHALE's acceptance test. A small abundance does not imply
a negligible production or loss flux, electron contribution, heating term, or
coupling to another row. EXHALE correctly accepts the constrained candidate
only after `full_network_residual` rechecks every active row, including rows
removed by the trace holdout and restored afterward.

The claimed factor `1.41` change in floor-cell occupancy has no reproducible
artifact or definition beside this memo and is not used as evidence. If an
abundance-weighted update norm is tested, retain the existing reaction and
element residual requirements unchanged.

### 5.2 Limit the step of trace species separately

Equilibrate applies a different step limit above and below a mole fraction of
1e-8 (`equilibrate_cea.f90:1719-1727`; `SIZE = 18.420681 = ln(1e8)`,
`9.2103404 = ln(1e4)`):

- above 1e-8: one step may not move a logarithm by more than 2, i.e. a factor
  of about 7.4, with a limit five times tighter on the change in total moles;
- below 1e-8 and increasing: the step is cut so the species cannot overshoot a
  mole fraction of 1e-4 in a single step.

This is a reasonable trust-region design candidate, but it is not a local
addition to the present call to MINPACK `hybrd`: that interface accepts one
scaled Euclidean trust radius and no component bounds. Implementing the rule
requires a solver with bound or component-step control, a transformed variable,
or an outer step filter with a consistent model update. EXHALE already removes
species below `eps_hold` from the Newton system, restores them from their own
balance rows, readmits them when necessary, and rechecks the full network.
Measure an additional benefit over that existing treatment before adding a
second trace-species policy.

### 5.3 Removing a degree of freedom changes what must still be checked

Equilibrate removes ionized species from the unknowns entirely below 750 K
(`equilibrate_cea.f90:805-810`):

```fortran
if (self%ions) then
   if (temp > 750) then
      N_atoms_use = self%N_atoms
   else
      self%remove_ions = .true.
      N_atoms_use = N_atoms_in
   end if
```

Photochem's gas-giant extension does **not** freeze the evolved composition
above the quench level. It applies estimated quenched profiles to the initial
guess when `initial_cond_with_quenching` is true, places the bottom of the
photochemical grid deeper than the deepest estimated quench point, and then
integrates the selected domain. VULCAN offers a manual convergence exclusion,
`conver_ignore`, which drops named species from the convergence test — used for
species that have no sink and therefore no steady state at all.

These are not precedents for silently accepting an underdetermined EXHALE
reaction row. A frozen-layer approximation would change the physical model
from local kinetic steady state to prescribed thermochemical equilibrium. It
requires an explicit timescale criterion, a validity range, continuity across
the handoff, and validation of charge and element fluxes. The current
trace-species holdout is narrower: it eliminates only small intermediates from
the Newton step and still requires their restored balance rows to pass.

### 5.4 Break the degeneracy of species at the floor

When a species density underflows, Equilibrate does not reset every such
species to the same value but offsets each by its own index
(`equilibrate_cea.f90:1576`):

```fortran
n_spec(i_reac) = 1d-13*(1d0 + 1d-6*DBLE(i_reac))
```

The original memo inferred that this index offset prevents a singular matrix,
but the cited Equilibrate source context is unavailable in the present tree.
In EXHALE it would also perturb otherwise symmetric seeds without addressing
the actual derivative formula. Current seed floors scale with the relevant
element density, and species below `eps_hold` are normally held out of the
Newton system. Do not add index-dependent physical densities unless a
reproducible rank test shows that identical floors cause a singularity and the
final root is invariant to the offsets.

---

## 6. Recommended order of work

1. **Resolve R12/R15 thermodynamic consistency.** DONE 2026-09-01. They are
   the same three-body channel in opposite directions, so physical detailed
   balance is the acceptance criterion. The recombination direction was
   selected (it is the one with a published range covering the layer) and the
   dissociation derived from a statistical-mechanics `K_eq`; see the R12
   comment in `mol_rates.f90`.
2. **Measure whether falloff matters.** Evaluate the reduced pressure across
   the molecular base using a documented falloff model and third-body
   efficiencies. If the low-pressure expression is valid throughout, record
   the range and make no numerical change.
3. **Replace the shifted relative difference with a real absolute floor.** Add
   a constrained-specific Jacobian routine using a scale-aware step and compare
   it with central differences through `cce_probe`. This is narrower than
   changing the shared `fdjac1` path. Retain the full reaction and conservation
   acceptance checks.
4. **Benchmark a complete derivative path.** A stoichiometry-driven analytic
   Jacobian or automatic differentiation should be compared against the probe
   oracle and representative full runs before changing the solver interface.
5. **Consider additional trace-step control only after measurement.** The
   current hold-place-readmit-full-recheck path already handles trace species.
   Any new component limit must demonstrate a benefit and cannot weaken
   reaction-balance acceptance.
6. **Treat a frozen layer as a physical model change.** Adopt it only with an
   explicit timescale criterion and a verified handoff, not as a convergence
   exception.

### What not to do

- **Do not claim that the local VULCAN and Photochem mechanisms validate the
  He <-> H rates.** They do not contain those reactions. This does not rule out
  validation from primary literature or another evaluated mechanism.
- **Do not replace `q_h2_equilibrium` only to imitate a reference code.** It is
  used by several initialization, boundary, and retry paths, so replacement
  requires a deliberate interface and physics review. Cost and accuracy should
  be measured on those callers first.
- **Do not move EXHALE's chemistry to stiff time integration merely to imitate
  another code.** A pseudo-transient or operator-split prototype is reasonable
  only if a focused benchmark shows a robustness or physical-model advantage.

---

## 7. Cross references

- `docs/molecular_chemistry_audit_he_rich.md` — the He <-> H rate pair and why
  its two directions are correctly independent.
- `docs/charge_exchange_cancellation_limit.md` — the diagnosis that led here.
- `docs/supersonic_molecular_base.md` — the problem this work belongs to.
- `docs/Update_EXHALE.md` sections 113-114 — the acceptance rule and the
  constrained solver.
- `src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90`,
  `src/modules/lower_atmosphere/mol_rates.f90`,
  `src/modules/lower_atmosphere/lower_column.f90` — the EXHALE side of every
  comparison above.

---

## 8. Review and validation scope

This review inspected the current EXHALE call sites and solver implementation,
the locally available VULCAN source and mechanisms, the locally available
Photochem source, and the checked `zahnle_earth.yaml` mechanism. It also
recomputed the R15/R12 table directly from the rate functions and stated
composition assumptions. The absent expanded Equilibrate source prevented a
new source-level check of the quoted CEA snippets; section 5 marks those claims
accordingly.

A focused `make cce_probe` build compiled the relevant EXHALE Fortran sources,
including `constrained_equilibrium_probe.f90`, but the final link failed. The
selected Conda LAPACK library requires `libgfortran.so.3`, which is not
available in the active compiler environment, and the linker reported missing
`GFORTRAN_1.0` symbols. Therefore no new `cce_probe` executable result was
measured in this review. A full hydrodynamic or regression run was not
performed because this review changes documentation only and does not modify a
code path.

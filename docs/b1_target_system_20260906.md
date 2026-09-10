# B1: the target governing system of EXHALE

**Status: accepted by the user 2026-09-06 as the Phase 2 target specification
(PLAN_20260906_rev2 Step B1).** The fourteen choices that were open when the
draft was written are decided as recommended; section 10 records each decision
with the question it answers, and every target and acceptance test below is
stated unconditionally on those decisions. The draft this was accepted from is
`docs/b1_target_system_20260906_draft.md`.

Step B1 of `docs/PLAN_20260906_rev2.md`, the Phase 2 gate of
`docs/development_plan_20260905_rev3.md` section 4.2. This document states
**what the code shall solve**. It is the counterpart of
`docs/d0_governing_system_20260906.md`, which states what the code solves
today; every target below names the D0 consistency item (`C<n>`) it replaces.

Scope. The parts of B1 that do not depend on the B1a active-equation
inventory were written first; the three sections that needed it (the
independent species space of T1.6, the validity states of the rollback
contract, the verification rows by configuration and the boundary equations by
cell class) are filled from
`docs/b1a_active_equation_inventory_20260906.md` sections 2 to 5, together
with the two B1 contents that depended on the same inventory: the oxygen
reaction and excitation-energy ledger (section 1.6) and the rule for the
`_adv` derived products (section 8).

**Line numbers drift, and the routine name is the citation.** The tree is
under edit by parallel steps of the same plan, and the files this document
cites moved within the hour it was written: the `ioniz_eq` call of the
marching loop shifted by 150 lines in `EXHALE_main.f90` between the first
draft of section 7 and 10:29 on 2026-09-06, and the module-level declarations
of `ionization_equilibrium` moved by about 88 lines when Step A1 landed. The
variable and routine names are the durable citation and are what an
implementation brief must grep for. **On 2026-09-07 (item DOCS-LINES) every
citation in this document was re-anchored to the file and the routine, or to
the labelled block of a long routine, read against the live tree.**

**Provenance.** Every code statement is READ from the tree at
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` on 2026-09-06. Nothing here was MEASURED: no build, no run, no test. Every
published equation is quoted by its number from the publisher PDF under
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/references/`.

**Owners.** Every target carries `Owner: (unassigned)`. The specification is
accepted, so the advisor assigns them with the implementation briefs.

---

## 1. Energy convention and the formation/excitation ledger

Replaces **D0 C1** (the temperature-preserving composition projection of the
marching source stage, `EXHALE_main.f90`; removed 2026-09-06, item B3c) and
**D0 C3** (the two-iteration semi-implicit
energy update with a lagged cooling, `solve_energy_semi_implicit`,
`energy_semi_implicit.f90`; replaced by the bracketed, residual-controlled
solve of item B2 on the same date). Answers D0 open questions 1 and 2, in the direction
recommended by `docs/To_be_determined_by_user_recommend_20260906.md`
section 5.3 rows "Energy variable" and "Composition projection".

### 1.1 The evolved variable

**T1.1.** The third conserved row stays thermal plus kinetic energy per
volume,

```text
U(3) = (1/2) rho v^2 + u_th ,
u_th = (3/2) n_mon k_B T + n_H2 [ (3/2) k_B T + u_rv(T) ] ,
```

with `n_mon` the number density of every particle treated as monatomic and
`u_rv(T)` the rovibrational internal energy per bound H2 molecule. This is
the variable the hydrodynamic fluxes, the caloric equation of state and the
Riemann solver already use (`caloric_eos.f90`, module header), and keeping it means
no flux, no EOS inversion and no boundary condition is redefined. A
formation-inclusive total energy is also a valid choice, but it requires all
of those to be reformulated together, and this specification does not take
that route.

### 1.2 The formation/excitation reservoir

**T1.2.** Define a second, non-transported energy density

```text
u_form = sum_s n_s eps_s ,
```

where `eps_s` is the formation plus excitation energy of one particle of
species `s` measured from a single declared reference state. The reference
is: every element as a neutral, ground-state, free atom at rest, at the zero
of the internal-level ladder. With that reference,

| species | `eps_s` | content |
|---|---|---|
| H I, He I, O I, C I, neutral metals | 0 | reference |
| e- | 0 | the ionization energy is carried by the ion |
| H II | `I(H)` | ionization potential |
| He II, He III | `I(He I)`, `I(He I) + I(He II)` | cumulative |
| metal ion of stage `k` | sum of the `k` ionization potentials from the neutral ground state | cumulative; the zero is the declared reference, never the element's dominant state at the base (decision 1) |
| He 2^3S | `E(2^3S)` above He I ground | metastable excitation |
| H(n=2) | `E_21` above H(1s) | excitation, both `2s` and `2p` |
| H2 | `-D0(H2)` | bound below two free H atoms |
| H2+, H3+, HeH+ | formation energy from the reference atoms | binding plus ionization |
| OH, H2O, CO | formation energy from the reference atoms | binding |

**T1.3 (the non-overlap rule).** Energy already held in the equation of
state is not also in the reservoir, and the reverse. The rule is enforced by
fixing one zero per internal ladder:

- `u_rv(T)` in `T1.1` is the rovibrational energy of **bound** H2 above its
  own `v = 0, J = 0` level. The module header of `caloric_eos.f90` states this explicitly
  ("Dissociation is not a heat capacity and is not in here: it is a chemical
  source term").
- `eps(H2) = -D0(H2)` is measured from `v = 0, J = 0` as well:
  `mol_rates.f90` carries `D0_H2_cm = 36118.11 cm^-1` labelled
  `D0(v=0,J=0)`.

The two therefore share one zero and cannot double count. The same
statement is required, and is the acceptance condition, for every other
species with a level ladder: exactly one of the EOS and the reservoir owns
each level's energy, and the specification names which.

**T1.4 (H(n=2) and He 2^3S).** Both are excited states with an explicit
population in the code (`n2s_arr`, `n2p_arr`, declared in `parameters.f90`;
`nheiTR`). Their excitation energy is in the reservoir, not in `u_th`: a
metastable atom carries `(3/2) k_B T` of translational energy like any other
particle and `eps_s` of stored excitation.

### 1.3 The local source-step identity

**T1.5.** For a local source step at fixed volume and fixed mass density,
with no transport and no mechanical work,

```text
Delta u_th + Delta u_form = integral_over_dt Q_external ,
```

with `Q_external` the net energy exchanged with the radiation field and any
other reservoir external to the cell's material system, under the local
boundary declared in `development_plan_20260905_rev3.md` section 4.2 item 3.
This is the schematic form of `To_be_determined_by_user_recommend_20260906.md`
section 5.3, adopted here as an exact requirement on the discrete step.

**T1.5 AS THE CODE IMPLEMENTS IT (B3c decision, instrumented by T15,
2026-09-07).** The form above is stated for a `Q_external` that carries only
exchanges with the radiation field. This code's `heat` and `cool` are not that
`Q_external`: they are the NET THERMAL source. Photoionization deposits
`h nu - E_th` and never `h nu`, so the potential the photon paid never passes
through the thermal pool; collisional ionization removes `E_th` from the
electrons as the `coio` term of `eval_cool`; radiative recombination lets the
potential leave as a photon and keeps only the electron's kinetic energy; and
every collisional reaction heat of the molecular, Penning, associative,
Lyman-Werner and oxygen channels is a difference of the one species
formation-energy table. With those sources the reservoir is already accounted
for, and the identity the discrete step has to satisfy is

```text
Delta u_th = integral_over_dt (heat - cool) ,
```

with `Delta u_form` a REPORTED quantity beside it and never a term of it:
adding it would count every collisional transfer twice and charge the gas for
the ionization energy the photons paid. Reaching the form stated above instead
is a redefinition of `heat` and `cool` across `util_ion_eq.f90`,
`Cool_coeff.f90` and `electron_energy_degradation.f90`, not an addition to the
row; that redefinition has not been made. An ATTEMPTED STEP is this local
source step plus transport, so the instrument that gates it
(`attempted_step.f90`) carries a MEASURED transport contribution on the right:
the marching loop marks the thermal energy on both sides of the source step,
and the residual is the closure of the source step alone.

The published precedent for this convention in a code of the same class is
Schulik and Booth (2023, MNRAS 523, 286) eqs. (53) to (55): their
photoionization heating and cooling terms carry `E_HI`, the average energy
injected per ionization, explicitly among the products, and the paper states
of those three terms that "It can be verified by summing the three terms
that the net heating rate is just the energy injected per ionization minus
the total cooling." Their evolved variable is likewise the internal energy
`rho_s e_s` (their eq. 42), not a formation-inclusive total.

### 1.4 The coupled temperature-composition source step

**T1.6.** The composition and the temperature of a cell are advanced by
**one** solve, not by a composition solve at frozen `T` followed by a
pressure rebuild followed by an energy solve.

Unknowns: `T^{n+1}` and the independent species densities `{n_i^{n+1}}` on
the independent space of T1.7.

Equations:

```text
(a)  n_i^{n+1} - n_i^n = dt * R_i({n^{n+1}}, T^{n+1}, field)        for each independent i
(b)  u_th(T^{n+1},{n^{n+1}}) - u_th(T^n,{n^n})
       + sum_s eps_s ( n_s^{n+1} - n_s^n )  =  dt * Q_ext({n^{n+1}}, T^{n+1}, field)
```

Constraints: element totals per nucleus, charge neutrality, and the total
mass condition of section 2, all applied at `t^{n+1}`. Row (a) is replaced
by an algebraic closure only for a species eliminated by a justified fast
process, which T1.7 enumerates.

`Q_ext` contains only exchanges with the radiation field and declared
external reservoirs: photoabsorption deposited in the cell, radiation
emitted and lost from the cell (section 3), and any specified external
reservoir flagged under B6. Every collisional reaction heat is **absent**
from `Q_ext`: it is already in the `sum_s eps_s Delta n_s` term of row (b),
which is the point of the ledger. A reaction heat added on top of the
reservoir term is a double count and the test of 1.5 catches it.

**Consequences, stated so they are not rediscovered later.**

- The pressure reset `comp_p_from_T` followed by `W_to_U` in the marching
  source stage of `EXHALE_main.f90` disappears as an independent update (done
  2026-09-06, item B3c). Pressure after
  the source step is whatever `T^{n+1}` and `{n^{n+1}}` give through the
  EOS, and the energy row is the constraint that produced them, so no
  unsourced `Delta u = (3/2) k_B T Delta n_part` can appear.
- The two-iteration structure of `energy_semi_implicit.f90` is replaced by
  a residual-controlled solve (B2): sources are evaluated at the returned
  temperature and composition, and the final residual of row (b) is tested.
- Splitting stays allowed. What is forbidden is a split whose pieces do not
  each satisfy their own instance of `T1.5`.

**T1.7 (the independent space, by configuration).** READ from
`docs/b1a_active_equation_inventory_20260906.md` sections 2.4, 2.5 and 5. The
local solve of T1.6 carries, in each physical cell, one temperature unknown and
`N_eq` species unknowns, `N_eq` being the stage-fraction count the
configuration selects (`input_read.f90`, `input_read`, the `N_eq` assignment
block). Everything outside that
list is eliminated by a named closure, and the constrained Jacobian of T1.6 is
assembled **after** every elimination, so it is square of dimension
`1 + N_eq` in a cell and carries no redundant row.

The eight configurations are the ones B1a section 5 generates from the
regression matrix and the examples:

| # | configuration | source case | independent species unknowns of the cell solve | `N_eq` | balances outside the cell solve |
|---|---|---|---|---|---|
| 1 | pure H/He | `examples/tutorial_nometals` | `x_HII, x_HeII, x_HeIII, x_HeITR` | 4 | none |
| 2 | He 2^3S and metals | `backup/regression/wasp_full`, `wasp_he23off` | the four above plus two stage fractions for each of C, N, O, Mg, Ca, Na, Fe | 18, or 17 with the triplet row absent | none; H(n=2) is a lagged closure, not an unknown |
| 3 | hot Uranus molecular | `backup/regression/mol_base_handoff` | `x_HII, x_HeII, x_HeIII, x_H2, x_H2p, x_H3p, x_HeHp, x_HeITR` | 8 | none |
| 4 | carriers | `backup/regression/mol_carrier` | the same eight, the H2 row replaced by the identity `x_4 - x_H2^fix` | 8 | one transported H2 balance |
| 5 | oxygen | `backup/regression/oxygen_chemistry` | the same eight plus `x_OH, x_H2O` plus two stage fractions for each of C, N, O | 16 | four transported balances (H2, OH, H2O, CO); CO carries no kinetic row (section 5) |
| 6 | diffusion | `backup/regression/mol_diffusion` | the eight molecular rows | 8 | one elemental balance for `X = rho_He/rho` |
| 7 | profile | `backup/regression/lower_profile` | `x_HII, x_HeII, x_HeIII, x_HeITR` plus two stage fractions for each metal element the profile carries | `4 + 2 n_melem` | the He/H elemental balance and one trace-element balance for each metal element, the latter with no counter-flux today (`C29`, section 4.4) |
| 8 | H+ transport | `backup/regression/hp_front` | the eight molecular rows, with both the H2 row and the H+ row replaced by identities | 8 | two transported balances (H2 and H+); the two base ghosts still re-solve H and H+ locally (`C18`) |

The eliminations, and what each one is:

| eliminated quantity | closure | imposed or checked | code site |
|---|---|---|---|
| the neutral stage of each element | the element total, neutral stage by difference | imposed by the layout of the fractions | `ion_system_HeH_mol` (`System_HeH_mol.f90`); `metal_fractions` (`ion_residual_core.f90`) |
| `n_e` | charge neutrality, `n_e = sum_s z_s n_s` | imposed as an identity; never an unknown and never a row | `calc_ne` (`utilities.f90`); `ion_system_HeH_mol` (`System_HeH_mol.f90`) |
| `rho` | none in the target: the density is supplied to the step and preserved (T2.1) | checked against the species mass sum (T2.2) | overwritten from `calc_rho` (`ionization_equilibrium.f90`) until item B3c on 2026-09-06 adopted T2.1 |
| O(1D) | local steady state on a fast process: one sink (O6) and a 1.6e-3 s chemical lifetime, so the flux through O6 equals the production `oj4 n_H2O` | imposed; the only fast-process elimination the code justifies at its site | `System_HeH_mol.f90`, the note above `oxygen_carrier_rows` |
| CO inside the cell solve | the transported value, clipped to the carbon and oxygen the cell has; the kinetic destruction row lives in the carrier operator, and the thermal equilibrium ceiling that stood here was deleted on 2026-09-06 (item CEILING-DEL) | imposed after the solve | `ioniz_eq` (`ionization_equilibrium.f90`); `carrier_source` and `limit_to_element_budget` (`diffusive_photochemistry.f90`) |
| H(n=2), both `2s` and `2p` | a 2x2 statistical equilibrium solved outside every solver, one outer pass before the state it closes | imposed, lagged | `n2_populations` (`excited_hydrogen.f90`), called through `excited_H_update` before `ioniz_eq` in the `coupled_source` loop of `EXHALE_main.f90` |
| a transported carrier inside the cell solve | the identity row `x_i - x_i^fix` carrying the transport operator's value | imposed | the end of `ion_system_HeH_mol` (`System_HeH_mol.f90`); `ioniz_eq` (`ionization_equilibrium.f90`) |

The three constraints, and what each removes:

- **element totals per nucleus** remove one unknown for each element, by
  writing its neutral stage as the difference. They are imposed, so they are
  not rows of the Jacobian and their residual is exactly zero by
  construction; what is checked is the nucleus count itself
  (`element_census.f90`, checked and never corrected).
- **charge neutrality** removes `n_e`. It is an identity at every iterate,
  not a row, so it contributes no equation and no Jacobian column.
- **total mass** removes nothing and adds nothing: under T2.1 it is an input
  to the step and a checked constraint (T2.2), not an equation of the cell
  solve. This is the one place where the target changes the space, because
  today the density leaves the sweep rewritten.
- a **transported** species keeps its identity row in the local system
  (decision 11), so that the Jacobian keeps its dimension and every other row
  sees the imposed value; its balance is a grid-coupled equation owned by B4
  and is not a local unknown. One layout serves every configuration, the row
  scaling and the acceptance classes do not change with the transport flags,
  and the identity row costs one trivially conditioned row. Removing the row
  and its column would give a smaller local system whose dimension changes
  with the flags, and it is not adopted.
- the **constrained continuation** inverts the roles: the unknowns become
  `u_k = ln n_k`, the element totals become explicit residual rows and charge
  neutrality holds identically at every iterate
  (`network_balance_rows` and `element_conservation_rows`,
  `constrained_chemical_equilibrium.f90`). The target
  requires the energy row (b) of T1.6 to be carried in that branch as well, on
  that same reduced space measured in `ln n_k`. A continuation branch that
  solves the composition alone re-introduces the `C1` projection through the
  back door, whatever section 1 does to the main path.
- **He 2^3S is an unknown** (`tr_triplet_row`, `ion_residual_core.f90`)
  and **H(n=2) is not**. Both carry their excitation energy in the reservoir
  of T1.4 either way; the difference is only which of them the Jacobian sees.

**T1.8 (one instantaneous evaluation, and the integrated budget beside it).**
Review 2 section 5.2. "One source evaluation" in T1.6 and in the final-state
assembly means one consistent instantaneous evaluation of the state and of its
radiation field, and it serves two things: the source outputs written beside
that state, and the endpoint checks (the residual of row (b), the acceptance
conditions of the physical step of A0). It does **not** mean that a substepped
or multistage update reports its energy budget from that endpoint evaluation.
Such an update keeps its **time-integrated** budget, assembled at its own time
levels, and what is reported as the step's heating, cooling and reaction
budget is that integral, never the endpoint rate multiplied by `dt`.

Both quantities use one definition of each physical process. Each carries, in
the output, its time level, its units, and whether it is an instantaneous rate
or an integral over the step; they appear as separate columns and are never
reconciled by redefining one as the other. The identity of T1.5 is stated on
the integrated budget, which is the quantity that must balance over the step;
the endpoint evaluation is what the residual test and the written sources see.

**Validity range.** `T1.5` holds for a fixed-volume, fixed-mass local step.
It does not constrain the transport operators, which are covered by B4:
there the requirement is that the thermal and chemical energy fluxes use the
same species face fluxes as the species continuity equations, so that their
sum is the material energy flux.

**Owner:** (unassigned).

### 1.5 Acceptance tests

**AT-1a (closed reacting cell).** One cell, no radiation field, no
transport, `Q_ext = 0` by construction, initialized off chemical
equilibrium and integrated to equilibrium. Required: `u_th + u_form` is
constant to round-off at every accepted step, over a run that crosses a
dissociation front (`n(H2)` from 0.8 to 0.0 of the H nuclei) and over a run
that crosses an ionization front. Any step of the current code fails this
by the `C1` projection amount, which is the RED reference.

**AT-1b (isolated photon event).** One monochromatic photon of energy
`h nu` absorbed by one H I atom, with every other process switched off.
Required: the ledger assigns `I(H)` to `Delta u_form`, `h nu - I(H)` to
`Delta u_th` through the photoelectron, and nothing else, with
`Delta u_th + Delta u_form = h nu` to round-off. With secondary ionization
active the photoelectron energy partitions further into additional
ionization potentials (reservoir), heat (thermal) and prompt radiation
(leaves the cell, so it appears in `Q_ext` with a negative sign), and the
test asserts that the three shares sum to `h nu - I(H)` and that no
recipient appears twice. The published statement of the same requirement is
Schulik and Booth (2023) section 3.4, where the C2Ray time-averaged rates
are adopted precisely so that "the number of photons absorbed in the cell
for heating is identical to the number of atoms that are photoionized", and
where a separately tracked high-energy optical depth prevents "double
counting of photons in cases of mixed ionizing and non-ionizing
absorption".

**AT-1c (recombination).** One recombination event returns `I(H)` from the
reservoir; the fraction radiated away appears in `Q_ext` and the rest in
`u_th`; the split is the one the adopted recombination cooling function
implies and is stated at one site.

### 1.6 The oxygen reaction and excitation-energy ledger

Replaces **D0 C25**: the oxygen cycle has photolysis heat but no collisional
reaction-energy ledger, the retained O(1D) electronic excitation is absent from
the heat sum, and the He(2^3S) associative branch heat is deposited nowhere.
B1a section 4.1 records that this state has no flag, counter or header field
today, so it is an unconditional property of any run with
`Oxygen chemistry: True`.

READ. The oxygen network is nine channels, listed at
the head of `System_HeH_mol.f90` with the identifiers of `docs/a2_reaction_audit.md`
section 2: `O1  OH + H2 -> H2O + H` and its reverse `O1r`, `O2  O + H2 -> OH + H`
and its reverse `O2r`, the three water photolysis branches `O3  H2O + hv -> OH + H`,
`O4  H2O + hv -> H2 + O(1D)` and `O5  H2O + hv -> O + H + H`, the hydroxyl
photolysis `O7  OH + hv -> O + H`, and `O6  O(1D) + H2 -> OH + H`, which is the
single sink that lets O(1D) be eliminated (T1.7). Of these, only the photolysis
channels deposit energy today: `heat_fuv` sums the excess of each absorbed band
photon over the bond energy for H2O and OH
(the `heat_fuv` block of `ioniz_eq`, `ionization_equilibrium.f90`, with the
O(1D) statement in the comment above it). The four collisional channels O1, O1r,
O2, O2r deposit nothing, and the code states at the same site that the 1.96 eV
of electronic excitation carried out of O4 by O(1D) is deliberately kept out of
the photolysis heat and that its release in O6 "is NOT deposited by this
network". The module header of `molecular_reaction_heat.f90` records the third gap: the
associative branch `He(2^3S) + H -> HeH+ + e` releases 8.1 eV and "is NOT
deposited anywhere", left out because the metastable density in the molecular
layer makes the channel 1e-13 of the sum there.

**T1.9 (the ledger, not a list of reaction heats).** Every oxygen species that
the run carries gets an `eps_s` in the table of T1.2, measured from the same
reference as every other species: O I and C I at zero, OH, H2O and CO at their
formation energy from the reference atoms, O(1D) at its electronic excitation
above O I. With those entries, the energy of each of the nine channels is the
difference of the reservoir contents of its reactants and products and needs no
separate oxygen heat term, which is the construction
the module header of `molecular_reaction_heat.f90` already states for the H2/He network
("Not as a list of reaction enthalpies, a list can disagree with itself, but
from ONE table of species enthalpies"). The requirement is therefore not a new
term but the extension of the one table to the oxygen and carbon species, and
the deletion of any channel-specific oxygen heat that would then be a second
count.

Three things the extension must state, because each is a place where the
present code loses energy:

1. **An elimination transfers the reservoir with the nuclei.** O(1D) stays
   eliminated by the fast closure of T1.7 (decision 13): its local steady
   state is exact, one sink and a 1.6e-3 s lifetime, so promoting it to an
   unknown would add a row and a column to every oxygen cell and change no
   answer. What the elimination must carry is the energy: the flux through O6
   carries both the oxygen nucleus and `eps(O(1D))` to the products at the same
   rate, and `AT-1d` (ii) is the guard that it does. A closure that moves the
   nuclei and drops the excitation energy is exactly the `C25` gap, and the
   transfer rule is general for every fast-process elimination, not a statement
   about oxygen.
2. **A photoevent names its recipients.** For O4 the absorbed photon pays the
   bond energy, `eps(O(1D)) - eps(O I)` goes to the reservoir and the remainder
   to the fragments' translation, exactly as `AT-1b` requires of a
   photoionization. The FUV branches O3, O5 and O7 follow the same rule with
   their own products.
3. **The associative He(2^3S) branch is a reaction of the same table.** With
   `eps(He 2^3S)` and `eps(HeH+)` present, its 8.1 eV appears as the difference
   of reservoir contents and is deposited without any channel-specific term. It
   stops being an omission the moment the reservoir exists, which is why it is
   listed here and not as a separate target.

**T1.10 (the exclusion rule and its observable).** Until T1.9 exists, a
configuration with `Oxygen chemistry: True` is excluded from validated physical
results, because its energy balance is incomplete by construction and no
tolerance can be stated for the gap. The state is **static**: it follows from
the configuration flag `thereis_oxychem`, not from a runtime event, so the
observable is the flag itself, recorded in the B6 validity state of the run and
carried into the output header and into every derived product (section 8). That
is what B1a section 4.1 asks for under "a static oxygen active implication".
Once T1.9 exists the exclusion is lifted and the observable becomes the
measured one: the closed-cell identity `AT-1d` below, run with the oxygen
network active.

**Acceptance test AT-1d (the oxygen ledger).** (i) A closed cell with the
oxygen network active, no radiation field, initialized off equilibrium:
`u_th + u_form` constant to round-off at every accepted step, over a run in
which the H2O and OH fractions change by an order of magnitude. RED today by
the whole collisional oxygen heat, which is absent. (ii) A single O4 event
followed by its O6 partner returns `eps(O(1D))` to the products and nothing
else, so that the pair deposits the same energy as the direct route
`H2O + hv -> OH + H` at the same photon energy. RED today, since the 1.96 eV
disappears. (iii) A single `He(2^3S) + H -> HeH+ + e` event deposits 8.1 eV,
RED today. (iv) With `Oxygen chemistry: True` and the ledger absent, the run
reports the static validity state and is refused certification, which is the
T1.10 half of the test.

---

## 2. Mass in the chemistry step

Replaces **D0 C2** (`ioniz_eq` declared the density `intent(inout)` and
overwrote it with `calc_rho` of the equilibrium composition, both in
`ionization_equilibrium.f90`, and the marching source stage of
`EXHALE_main.f90` pushed that into the conserved mass row; the code's own
`EXHALE_UPDATE_MAP` reports it as a `chem` column absent from the steady
residual, `update_map_end_step`). **Adopted 2026-09-06 (item B3c):** the
density is now `intent(in)` and `calc_rho` builds the optional `rho_recon`
check. Answers D0 open question 3 in the direction
recommended in section 5.3, row "Mass density".

**T2.1.** Chemistry preserves the mass density supplied to the local source
step:

```text
rho^{n+1} = rho^n   across the source step,   d rho / dt |_chem = 0 .
```

The physical content is that chemical reactions rearrange nucleons and
electrons among species and create no mass. It is not an approximation and
needs no tolerance in the equation.

**T2.2.** The species mass sum is **verified** against the supplied density,
not used to redefine it:

```text
| sum_s n_s m_s - rho | / rho  <=  tol_mass ,
```

with `m_s` from the one species mass table. **The declared convention
(decision 2) is that the electron mass is carried with its ion**: `m_s` for an
ion of stage `k` is the neutral atomic mass, and the free electrons carry no
mass in the sum. An ionization then moves no mass between species, the mass
sum is independent of the ionization state, and the omitted quantity is
`n_e m_e`, so `tol_mass` follows from the convention as `n_e m_e / rho`
evaluated on the state, which is at most `m_e/m_H` of order `5e-4` in a fully
ionized hydrogen gas and smaller everywhere else. The convention is stated at
exactly one site, `tol_mass` is derived from it and is not tuned until an
existing state passes. A violation is a failure of the composition solve
(reported through B6), never a correction applied to `rho`.

**T2.3.** Reconstruction of `rho` from the composition is a diagnostic. The
steady equations gain **no** chemical mass source: the fix is to remove the
term from the production map, not to add a matching fictitious term to the
residual, so that the marching map and the stationary map have the same
fixed point.

**Interface consequence.** The `n_io` argument of `ioniz_eq` becomes
`intent(in)` (decision 3, confirmed 2026-09-06). That is a public-interface
change inside the source tree: the callers in the marching loop, the steady
residual and the post-process all pass a density that the routine may no
longer write, and each of them keeps the density it supplied. The
reconstructed density stays available as a separate returned diagnostic for
T2.2.

**Owner:** (unassigned).

**Acceptance test AT-2.** (i) The update map with the chemistry stage active
reports an identically zero `chem` column in the mass row, for a
configuration in which the composition changes strongly over the step
(hot-Uranus molecular base). RED today by the `C2` amount. (ii) The
assembled stationary equations contain no chemical mass source term:
a source inspection test over `steady_residual.f90`. (iii) `AT-2c`: a cell
whose composition is driven through H2 dissociation returns
`| sum_s n_s m_s - rho | / rho <= tol_mass`, and an injected species mass
table error of one part in `10^3` makes it RED.

---

## 3. Local energy ownership of radiation, Lyman-alpha included

Replaces **D0 C7**: the collisional excitation cooling function charges the
full `E_21 = 10.2 eV` per H I excitation as escaping radiation with no
escape probability (`lambda_coex_HI`, `Cool_coeff.f90`, the Cen-type fit
`7.5e-19/(1+sqrt(T/1e5)) exp(-118348/T)`, summed into `coex` by
`eval_cool_cells` in `util_ion_eq.f90`), in the same cells where `lya_rt`
computes an escape probability far below one, while the collisional
de-excitation heat returns `E_21` for each de-excitation (`excited_H_update`,
`excited_hydrogen.f90`,
`Hdx_arr = ne * E21_erg * (c2s1s n2s + c2p1s n2p)`).

### 3.1 The rule

**T3.1 (local ownership, rev 3 section 4.2 item 3).** "Escaping" means
leaving the **local material system of the cell**, not leaving the
atmosphere. Emission is a loss of the emitting cell; absorption is a gain of
the absorbing cell. For an explicit radiation field, the field's energy
change and its boundary flux enter the domain budget. For an eliminated
diffuse field or a local recycling approximation, the corresponding **net**
material exchange is derived, and the same photon is never reintroduced as
an independent source.

### 3.2 What beta is in this code

READ, `lya_rt.f90`, the module header and `jlya_escape_prob`. The code builds

```text
<N>      = (4 sqrt(6)/pi^2) [(tau + tau_dn)/2] Phi(tau_dn/(tau + tau_dn))
beta_esc = 1 / (1 + <N>)
beta_sob = min(1, (1 - exp(-tau_S)) / tau_S) ,     tau_S = C_sob n_HI / |dv/dr|
beta_tot = beta_esc + beta_sob (1 - beta_esc)
```

with `tau` the line-center optical depth to the outer boundary and `tau_dn`
the depth to the planet-ward face of the trapping slab: `tau` itself for the
default reflecting boundary (the emitting plane is then a symmetry plane), the
rest of the column when `lya_bottom_absorber` declares the bottom a pure sink.
`<N>` is the mean number of scatterings before escape of the static
plane-parallel damping-wing slab, Neufeld (1990) eq. (3.27) at zero continuum
destruction, which for the mid-plane source is Harrington (1973) eq. (40),
`<N> = 0.909316 tau`; `Phi` is its source-position factor, `Phi(1/2) = u_2`,
Catalan's constant. `<N>` counts absorptions, so `beta_esc = 1/(1 + <N>)` is
the escape chance per emission, exact in the thin limit as well as in the
wings. `beta_tot` is therefore the probability that a line photon created in
the cell **leaves the resonance**, by the static slab channel or by the
velocity-gradient channel. It is not a probability of leaving the atmosphere,
and the escape-probability closure of the same routine uses it as
`Jbar = S (1 - beta_tot)` (the `Jint` assembly of `jlya_escape_prob`). The two FACES of the slab are
not two channels: they share one escaping population, split by Neufeld
eq. (2.25) in the ratio `tau_dn/(tau + tau_dn)` star-ward.

That distinction is what makes the present arrangement a double charge: the
Cen cooling function already assumes every collisional excitation ends in an
escaping photon, so `1 - beta_tot` of them are charged to the gas although
they are re-absorbed locally, and the fraction of those that end in
collisional de-excitation are additionally credited back by `Hdx_arr`. The
same photon appears in two independent terms.

### 3.3 The target: one 2p budget, one exchange term

**T3.2.** For a cell in which an explicit `H(n=2)` population is solved, the
material energy exchange with the Lyman-alpha line is derived from that
population's own balance and from nothing else. Let `n_2p` satisfy the
steady balance already assembled in the code
(the `Jint` denominator of `jlya_escape_prob` in `lya_rt.f90`, and
`excited_hydrogen.f90`):

```text
n_2p [ A_2p1s beta_tot + D_2p ]  =  P_coll + P_rec + P_pump ,

P_coll = C_1s2p(T) n_e n_1s                (collisional excitation)
P_rec  = alpha_2p(T) n_e n_p               (recombination cascade into 2p)
P_pump = n_1s B_12 Jbar                    (absorption from the local field)
D_2p   = n2p_destruction_rate(T, n_e, ...) (collisional de-excitation,
                                            n = 2 photoionization,
                                            l-mixing followed by two-photon decay)
```

The destruction set is the one Huang et al. (2023, ApJ 951, 123) section 2.7
defines: "The scattering Lya photons are considered destroyed if H(2p) is
photoionized by Balmer continuum photons or deexcited by electron collisions
instead of radiative decay. If the collisional l-mixing to H(2s) is followed
by photoionization, collisional deexcitation, or two-photon decay, the Lya
photons are also considered destroyed." The code already implements
`D_2p` as `n2p_destruction_rate` and already feeds `beta_tot` and `D_2p` into
the same denominator (`jlya_escape_prob`), so the target reuses one
population, not a second one.

**T3.3 (the exchange term).** The net material energy loss to the
Lyman-alpha line, per unit volume, is

```text
Q_Lya = E_21 [ A_2p1s beta_tot n_2p  -  P_pump ]
```

with every other 2p channel appearing where it belongs and nowhere else:

- `P_coll` removes `E_21` from `u_th` and adds it to `u_form` through
  `eps(H(n=2))` of section 1. It is **not** a cooling term.
- collisional de-excitation returns `E_21` from `u_form` to `u_th`. It is
  **not** a heating term added on top; it is the reverse of the line above,
  inside the same ledger.
- `n = 2` photoionization moves `eps(H(n=2))` plus the absorbed photon
  energy into the H II reservoir and the photoelectron, under section 1's
  rules; the Lyman-alpha field loses the photon that made the excitation and
  that loss is already in `beta_tot`'s denominator.
- two-photon decay emits into a continuum, not into the line; it is a loss of
  `E_21` from `u_form` to the radiation field and appears in `Q_ext` with its
  own label, not inside `Q_Lya`.
- `P_pump` is the absorption of the local field, stellar plus internal. It
  is the gain half of the exchange and must appear with the same `Jbar` the
  escape-probability closure used, or the closure is not self-consistent.

`Q_Lya` then enters `Q_ext` of `T1.5`, and the H I collisional excitation
term of the bulk cooling function drops its `1s -> 2p` share, keeping only
the transitions for which no explicit level population exists.

**Decision 4 (2026-09-06).** That remaining term stays the Cen-type fit as it
stands, used for the higher levels only, with its Lyman-alpha share removed
and the removal stated at the code site; it is not replaced by a new fit. The
fit was never published as a Lyman-alpha-excluded quantity, so the removal is
an approximation of the target and is recorded as one: the `1s -> 2p` share of
the fit is subtracted with the same rate coefficient the 2p balance uses for
`P_coll`, so that one collisional excitation rate serves both and the
subtraction cannot drift from the term it removes. Where no explicit `H(n=2)`
population exists (`use_excited_H` false), the full fit is kept and the
configuration carries the `C7` validity state of section 9.1.

**T3.4 (explicitly: a multiplier is not the target).** Replacing
`coex_HI` by `beta_tot * coex_HI` is **not** the target specification. It
gets the escaping fraction right and leaves the trapped fraction with no
owner: the trapped photon is re-absorbed, its energy re-enters the gas
through collisional de-excitation or leaves through a destruction channel,
and both are already terms of the 2p balance. A multiplier on the old
cooling term keeps `Hdx_arr` as an independent heating source and therefore
keeps the double count that `C7` names.

The published models bracket this. Salz et al. (2016, A&A 586, A75) treat
Lyman-alpha cooling with an escape-probability radiative transfer (their
section 3.1: "Radiative transfer is approximated by the escape probability
mechanism (Castor 1970; Elitzur 1982)"), and their section 4 states the
physical dependence the target must reproduce: "Lya cooling is further
increased by the now smaller neutral hydrogen column density above the
emission layer, which increases the escape probability for the radiation."
Huang et al. (2023) instead solve the field by Monte Carlo and take the
level populations from it (their eq. 11, `n_2p = n_1s (g_2p/g_1s)
(c^2 / 2 h nu^3) JLya`), iterating with the hydrodynamic model. EXHALE's
closure is the first kind with the second kind's destruction bookkeeping;
`T3.2` and `T3.3` are what makes it a ledger rather than a multiplier.

**Validity range.** The damping-wing slab solution requires
`(a tau)^(1/3) > 10` (Neufeld 1990 section Va); `lya_wing_domain_record`
counts the cell visits that fail it, and there the closure is carried by its
thin limit `beta -> 1`. The Sobolev channel requires the Sobolev length to be
short compared with the gradient scale of the source function. The static and
Sobolev channels are combined as independent escape routes, which is an
approximation and is stated as one at the code site.

**Owner:** (unassigned).

### 3.4 Acceptance test AT-3 (two-cell radiation exchange)

Two adjacent cells, an explicit line field, every other process off. Cell 1
is given a 2p population and emits; the optical depth between them is set so
that a specified fraction `f` of the emitted photons is absorbed in cell 2
and the rest leaves the domain. Required, at every accepted step:

```text
Delta (u_th + u_form)_1 + Delta (u_th + u_form)_2 + E_escaped = 0
```

to round-off, with `E_escaped` the boundary flux the field reports. Two
limits pin it: with `f = 1` the pair is closed and the material total is
unchanged; with `f = 0` the material loses exactly `E_escaped`. A build in
which the collisional excitation is charged as escaping cooling **and** the
de-excitation returns `E_21` as an independent heat source fails both limits,
which is the RED reference. Extension: repeat with the bottom-absorber
boundary declared, where `beta_bot` sends photons into a specified external
reservoir; the reservoir's gain is then a labelled term of `Q_ext`, not a
disappearance.

---

## 4. Charged transport closure

Replaces **D0 C28** (element diffusion in advective, not conservative form),
**D0 C29** (zero net diffusive mass is structural for the H/He binary and is
not enforced for the trace metals, which move against a fixed hydrogen
background with no counter-flux, `solve_trace_element_in_hydrogen`,
`binary_element_diffusion.f90`),
and the transport half of **D0 C21**. Answers D0 open questions 4 and 13, in
the direction recommended in section 5.3 rows "Charged transport" and
"Trace-metal diffusion". Follows rev 3 section 4.2 item 5, which requires
that the closure be **derived** for the EXHALE model and not attributed to
the paper.

### 4.1 The published equations

Koskinen et al. (2022, ApJ 929, 52) Appendix B. The species diffusion
equation with the electric term is their eq. (B9),

```text
(1 + Lambda_s) dx_s/dr + ( x_s - rho_s/rho ) dln(p)/dr - n_s e_s E / p
    = - sum_{t /= s} ( x_s x_t / D_st ) ( w_s - w_t )                    (B9)
```

under the condition of their eq. (B10),

```text
sum_s rho_s w_s = 0 ,                                                    (B10)
```

with `x_s` the volume mixing ratios, `Lambda_s = K_zz / D_s`, `K_zz` the eddy
diffusion coefficient, `D_st` the mutual diffusion coefficient, `e_s` the
electric charge, `E` the radial electric field, and the net molecular
diffusion coefficient given by their eq. (B11), `1/D_s = sum_{t /= s}
x_t / D_st`. Their `w + w_s` is the species velocity and `w` the center of
mass velocity (their eq. B4). The paper states the assumptions attached:
"We assume that the planetary magnetic field is negligible and ignore
thermal diffusion since the coefficients for the latter are poorly known at
the relevant temperatures."

**The paper states no current condition.** (B10) is a mass condition. The
closure that fixes `E` is derived below for EXHALE's own model and is not
attributed to Koskinen et al.

### 4.2 The current condition, derived

**T4.1 (the stated condition).** EXHALE imposes no external current and
carries no magnetic field. Take charge neutrality, which the species solve
already enforces exactly,

```text
sum_s z_s n_s = 0 ,                                                       (N)
```

so the charge density `rho_q` vanishes identically. Charge conservation then
gives, in spherical symmetry,

```text
d rho_q / dt + (1/r^2) d( r^2 J ) / dr = 0   =>   J(r) = C / r^2 ,
J = sum_s z_s e n_s ( w + w_s ) = e w sum_s z_s n_s + sum_s z_s e n_s w_s
  = sum_s z_s e n_s w_s      (by N).
```

`C` is fixed by the boundary: no current is injected at the inner boundary
and none is drawn at the outer one, so `C = 0` and

```text
sum_s z_s e n_s w_s = 0 .                                                 (J0)
```

**(J0) is the adopted closure of the EXHALE model**, derived from
neutrality, charge conservation, spherical symmetry and a stated
no-injected-current boundary. It is not a universal identity: a model with an
imposed current `I` through the shell replaces `C = 0` by
`J = I / (4 pi r^2)`, and the closure becomes an inhomogeneous condition
whose solution shifts `E` by the field required to drive `I` against the
total Coulomb friction. Any future magnetospheric coupling that supplies `I`
therefore changes this one equation and nothing else in the derivation.

### 4.3 The ambipolar field

**T4.2.** Write the electron member of the (B9) family. With `m_e` negligible
compared with every ion mass, electron inertia and electron gravity drop, and
the electron momentum balance reduces to

```text
n_e e E = - d p_e / dr + R_e ,
```

with `R_e` the friction force density on the electrons from the ions. Under
(J0) the ion and electron drifts carry no net current, so `R_e` is the
resistive term `eta J` evaluated at `J = 0` and vanishes at the level of the
closure. Hence

```text
e E = - (1 / n_e) d p_e / dr = - k_B T dln( n_e T ) / dr .                (F)
```

This is exactly the field the code computes:
`settling_coefficient` (`binary_element_diffusion.f90`) builds
`eEf = -kb_erg * TK * dln(n_e T)/dr` by central difference, one-sided at the
ends, gated by `he_ambipolar`; the routine's own header states the same
expression. The target therefore **adopts the
field the code already has** and supplies the derivation it was missing,
together with the conditions under which it holds.

**Omitted order of (F).** Three terms are dropped and each is named:

1. electron inertia, `O(m_e/m_i)` relative to the retained pressure term;
2. the thermal force, `alpha_T dlnT/dr`, which the code carries as a separate
   optional term with `he_alphaT` defaulting to zero
   (the `he_alphaT` term of `settling_coefficient`,
   `binary_element_diffusion.f90`), and which Koskinen et al.
   (2022) Appendix B also drop with a stated reason;
3. the resistive term, identically zero under (J0) and nonzero only in a
   model with an imposed current.

The target requires these three to be printed as the closure's stated domain
and evaluated as a diagnostic, not assumed small.

### 4.4 The diffusion velocities and the background response

**T4.3.** For the EXHALE species set, (B9) is written for every species that
diffuses, neutrals and ions together:

- neutrals: H I, He I, He 2^3S, H2, OH, H2O, CO, neutral metals (`e_s = 0`,
  so the electric term is absent, which is why the present carrier operator
  can omit it: the module header of `diffusive_photochemistry.f90` states that "All four
  transported species are neutral, so eq. (2) of the design has no
  (Z e E)/(k T) term for any of them; the field is not neglected, it is
  absent");
- ions: H+, He+, He++, H2+, H3+, HeH+, metal ions of every stage, each with
  its `e_s` and the field (F);
- electrons: eliminated by (N) and (J0), which is what produced (F).

The system is solved for `{w_s}` by inversion of the friction matrix on the
right-hand side of (B9), as Koskinen et al. (2022) Appendix B prescribe
("The diffusion velocities are obtained by matrix inversion from (B9) under
the condition that (B10)").

**T4.4 (the background response, the D0 open question).** (B10) is the
closure that makes the complete diffusive mass flux vanish. It is not
automatic for a trace species solved against a frozen background. Today,
`solve_trace_element_in_hydrogen` (`binary_element_diffusion.f90`)
moves a metal element through a fixed hydrogen background with no
counter-flux, so `sum_s rho_s w_s /= 0` by the amount the metals carry. The
target requires the background response

```text
sum_s rho_s w_s = 0    =>    for a single background carrier b :
    rho_b w_b = - sum_{X in trace} rho_X w_X ,
```

that is, the hydrogen (and, where helium is resolved as a second major
component, the pair) carries the exact counter-flux, distributed over the
background components in proportion to their share of the friction so that
(B9) is satisfied for those components too.

**Decision 5 (2026-09-06).** The target solves the trace elements **jointly
with the H/He binary in one matrix inversion**, which satisfies (B10) by
construction and derives the counter-flux with the charged closure (J0) rather
than assuming it. The sequential solve against a frozen background is retained
only as an interim, and only with its residual
`| sum_s rho_s w_s | / ( sum_s rho_s |w_s| )` measured, bounded at a stated
order and reported with a stated validity domain; an interim result carries
that measurement or it is not a validated result.

**T4.5 (conservative form).** The element and carrier transport equations are
written in conservative flux form on the spherical cell volume, with the same
face fluxes as the mass equation. That is B4's business; it is stated here
because `C28` (advective form, helium mass conserved only to hydro truncation
in a transient) and `C21` (cell-centered derivative, `r_center^2 dr`, base
direction from a wind-region average) are diffusion items and their target is
this one.

**T4.6 (the common-velocity approximation and its domain).** Setting
`w_s = 0` for a species, which is what EXHALE does today for H+ (D0 open
question 4: the proton has zero molecular diffusion and is neither
diffusing nor consistently comoving), is admissible only where the species
is advection-dominated. The criterion is the Peclet number on the local
gradient scale,

```text
Pe_s = L |w| / ( D_s + K_zz ) >> 1 ,   L = | dln x_s / dr |^{-1} ,
```

and the omitted relative error in the species mixing ratio is `O(1/Pe_s)`.
Koskinen et al. (2022) Appendix B make the converse statement about their
own approximation: "The diffusion approximation does not place constraints
on the bulk (center of mass) velocity w that can be either subsonic or
supersonic. It does place constraints on the magnitude of ws, but it can be
shown that a violation of the conditions in which the approximation is valid
is very unlikely to occur and generally requires minor species velocities to
be significantly faster than the bulk flow velocity." The target requires
`Pe_s` to be evaluated in each cell for every species held at `w_s = 0` and
reported through B6 where it falls below a stated threshold.

**Decision 6 (2026-09-06).** In the target H+ **enters the (B9) set with its
charge**, together with every other ion, because the ambipolar field (F) was
derived from a closure (J0) that sums over all charged species and is not
consistent with holding the most abundant ion at `w_s = 0`. Until B4 lands
that set, the present `w_s = 0` treatment is retained as an interim with the
`Pe_s` diagnostic reported in every cell, and a run in that interim makes no
claim about ion separation; the change to every ionized-region profile is
expected and is not a reason to defer it, since the old profiles were produced
by a proton that was neither diffusing nor consistently comoving.

**Owner:** (unassigned).

### 4.5 Acceptance tests

**AT-4a (zero net diffusive mass).** A static isothermal column with a trace
metal and no wind. Required: `| sum_s rho_s w_s |` divided by
`sum_s rho_s |w_s|` is at or below the stated order at every cell and every
step. RED today for a metals-on run by `C29`.

**AT-4b (the current condition).** The same column with several ion stages
present. Required: `| sum_s z_s n_s w_s |` divided by
`sum_s |z_s| n_s |w_s|` is at or below the stated order, this being the
closure the field was derived from and therefore a check on the
implementation, not on the physics.

**AT-4c (the ambipolar limits).** The relative settling mass `dmeff` that
`settling_coefficient` returns reproduces the three analytic limits the
header of `settling_coefficient` names (`binary_element_diffusion.f90`):
`m_He/m_H - 1` in neutral gas, `m_He/m_H - 3/2` in an H+ plasma,
`2 m_He/3 m_H - 1` in a He++ plasma, that is 2.9715, 2.4715 and 1.6477 with
`m_He/m_H = 3.9715`.

**AT-4d (element conservation in a transient).** A diffusion-only column run
through a transient conserves each element's total nucleus count to
round-off, which the conservative form of T4.5 gives and the advective form
of `C28` gives only to hydro truncation.

---

## 5. CO: exclusion adopted, the destruction model held open

Replaces **D0 C26** (CO has no balance row and no kinetic destruction rate:
an inert transported reservoir until an algebraic equilibrium ceiling removes
it). Follows rev 3 section 4.2 item 9 and D0 open question 5. Both halves of
C26 have since been answered in the code: the destruction row is section 5.2a
(item B3b-CO) and the ceiling is deleted (item CEILING-DEL).

### 5.1 What the code did when this section was written, and what it does now

READ, as of 2026-09-06 10:29 KST. With the oxygen chemistry on, the CO
density was either the local equilibrium value
`co_equilibrium_density(n_C, n_O, T)` or, when carrier transport is on, the
transported value, and it was then clipped to the carbon and oxygen the cell
has:

```text
nCO = min( nox_eq(:,3), n_el(C), n_el(O) )         (ioniz_eq, ionization_equilibrium.f90)
n_tot(O) = max( n_el(O) - nCO, 0 ) ,  n_tot(C) = max( n_el(C) - nCO, 0 )   (same block)
```

with a second, thermodynamic ceiling in `limit_to_element_budget`
(`diffusive_photochemistry.f90`) whose cumulative record was four quantities:
`co_ceiling_applications` and `co_ceiling_cells_hit` (cell applications and
distinct cells), `co_ceiling_CO_removed` (molecules removed) and
`co_ceiling_seen` (a flag for each cell).

READ, current tree. The clip against the cell's carbon and oxygen is
unchanged. The thermodynamic ceiling and all four counters were **deleted on
2026-09-06 (item CEILING-DEL)**; `limit_to_element_budget`
(`diffusive_photochemistry.f90`) now applies conservation only, the
CO abundance is set by the destruction row of section 5.2a, and what is
recorded is the destruction model's domain (`carrier_co_domain_take`, same
file) rather than a clamp count. MEASURED at the deletion: on
`oxygen_chemistry` at 1000 steps the ceiling never fired
(`co_ceiling_cells 0`, `co_ceiling_applications 0`,
`co_ceiling_CO_removed 0.0`), so every output column is byte-identical
across it.

### 5.2 Alternative A: constructible above the helium ionization front (2026-09-06, revised the same day)

**This section is the record of the argument that led to decision 7, and its
present tense is that of 2026-09-06.** Alternative A was implemented the same
day (section 5.2a, item B3b-CO) and the ceiling it argues against was deleted
(item CEILING-DEL), so every sentence below that speaks of the ceiling as
acting or firing is describing the tree as it then stood. The measurements are
left exactly as they were made.

**First reading (morning).** With the neutral channels alone (Westley 1980
NSRDS-NBS 67: `CO + H -> C + OH`, `CO + M -> C + O + M`, barriers of 6.7 eV
and more) `tau_dest << tau_res` fails everywhere the ceiling acts, by a
factor 9 to 1e8; `tau_dest << tau_form` holds wherever it fires (detailed
balance, `tau_dest/tau_form = K_eq/Q`).

**Revised (afternoon), with the UMIST RATE22 file** (Millar et al. 2024, A&A
682, A109; `~/RT_Codes/UMIST/rate22_final.rates`, 68 CO-reactant entries
quoted in `docs/co_destruction_rates_literature_20260906.md` section 12):
entry 4068, `He+ + CO -> C+ + O + He`, `alpha = 1.60e-9 cm^3 s^-1`, `beta =
gamma = 0`, measured, accuracy A, temperature independent, every product a
species EXHALE carries. RATE22 has no `H+ + CO`, no `CO + M` and no `CO +
H2` destruction entry. Combined with the Visser et al. 2009 shielding on the
run's own columns and measured against the run's `min(tau_adv, tau_diff)` on
the `oxygen_chemistry` snapshot: `tau_dest/tau_res` is 6.7e2 at the base
(He neutral, thick columns), 1.5 at 1.03 Rp, 5.7e-4 at 1.09 Rp, 9.6e-5 at
1.14 Rp (carried by He+, `n(He+)` from 0.05 to 1.3e7 cm^-3), 7e-2 at 3.4 Rp.
**The one-sided destruction model is therefore constructible above the
helium ionization front and not below it; its domain boundary is
compositional (the He+ front), not thermal (the `K_eq` surface the present
ceiling switches on).** Caveat stated in the document: the snapshot's cells
never exceed 1456 K, so the nesting of the two surfaces where the ceiling
fires (3000 to 4000 K) is inferred from the structure and must be measured
on a converged run reaching those temperatures.

**Decision 7 as it now reads.** Exclusion (alternative B) stays in force
until T5.1 is implemented, and for a sharper reason than before: the present
ceiling is a thermodynamic switch standing in for a kinetic process whose
switching surface is not the destruction domain's boundary. T5.1 becomes an
implementable increment (B3b-CO, after B3c): `He+ + CO -> C+ + O + He` at the
RATE22 rate with its reaction energy from the one formation-energy table,
CO photodissociation as a fourth absorber on the Lyman-Werner beam with the
Visser shielding on a CO column (the `lyman_werner.f90` pattern; Tsai et al.
2021 and Moses et al. 2011 implement no CO self-shielding, so EXHALE's own
H2 treatment is the precedent), the domain stated as `tau_dest << tau_res`
evaluated on the state, and the ceiling removed. Also recorded there: CO2
is not in the species table (removes `CO + OH -> CO2 + H`), and the two
enthalpy tables' zeros are now one table (B3b-ENTH).

#### What alternative A would have needed (record)

**T5.1.** CO gains a balance row of the same form as the other carriers,

```text
d n_CO / dt + (1/r^2) d[ r^2 ( n_CO v + Phi_CO ) ] / dr = - L_CO ,
L_CO = sum_k k_k(T) n_CO n_{M_k} + J_CO n_CO ,
```

with named destruction partners `M_k` and a named photodissociation rate
`J_CO`, the formation terms deliberately omitted (that is what "one-sided"
means), and the products of every channel named so that C, O and H nuclei
close. The reaction enthalpy of each channel enters `u_form` of section 1
through the products' `eps_s`; no separate CO heat term is added.

**T5.2 (the domain argument, which is the part that makes it more than a
diagnostic).** Omitting formation is legitimate only where formation is slow
compared with the residence time, and instantaneous destruction is a
legitimate ceiling only where destruction is fast compared with it:

```text
tau_dest = n_CO / L_CO  <<  tau_res  <<  tau_form ,
tau_res  = min( L / |v| , L^2 / (D_CO + K_zz) ) ,   L the local gradient scale .
```

Rev 3 section 4.2 item 9 states the requirement this encodes: "An
equilibrium constant alone does not give the destruction timescale; the
domain argument is required for the model to be more than a diagnostic."
The two inequalities must be shown to hold, from the adopted rate
coefficients, over the whole region in which the present ceiling is active,
and the region where they fail must be reported.

**Data requirement.** T5.1 needs published `k_k(T)` with their temperature
ranges and `J_CO` with its band data. This specification does **not** supply
them:
no rate compilation for CO destruction in this temperature and density
regime was consulted here, and a specification that names rates without
reading their source would violate the citation rule. This is a blocked
part, listed as such in the report.

### 5.2a Alternative A implemented (2026-09-06, item B3b-CO)

T5.1 and T5.2 exist in the code. The CO carrier row is no longer
`src(ic_CO) = 0`: it carries

```text
S_CO = - [ k_D1 n(He+) + k_CO(N_CO, N_H2, tau_cont) ] n_CO ,
k_D1 = 1.60e-9 cm^3 s^-1     (UMIST RATE22 entry 4068, measured, accuracy A)
k_CO = <(F_LW/<hv>) sigma_CO Theta(N_CO,N_H2) exp(-tau_cont)>_cell ,
```

with `sigma_CO = 1.0160e-17 cm^2` per photon of the 912-1201 A beam
(MEASURED from the Heays et al. 2017 cross section as concatenated by
Photochem, integrated over the 912-1118 A window CO absorbs in) and `Theta`
the 12CO block of Visser, van Dishoeck and Black (2009) Table 6, log-bilinear
in the two logarithms with the edge value outside the grid. The rate is the
mean over the cell on the same composite quadrature the H2 rate uses, for
the same reason: `Theta` falls by more than three decades along the CO column
axis and by more than six along the H2 axis.

Both reaction energies are differences of the one formation-energy table:
`q(D1) = eps(He+) + eps(CO) - eps(C+) - eps(O) - eps(He) = +2.2117 eV`
(MEASURED) and `D0(CO) = 11.1157 eV` (MEASURED), the second shared with the
photodissociation channel so the two channels of one molecule cannot state
two bond energies. The photon energy of a dissociation event is the
oscillator-strength weighted mean over Visser's 37 bands, 12.8674 eV
(MEASURED from their Table 1), so the fragments keep 1.7517 eV. Both
deposits are named channels 18 and 19 of `heating_of_composition` and are
formed nowhere else.

T5.2 is the domain record: `tau_dest = 1/(k_D1 n(He+) + k_CO)`,
`tau_res = min(r/|v|, dr^2/(D_CO + K_zz))` from the code's own definitions,
and `tau_form` from RATE22 entry 8597 for the record only. A cell is out of
domain when `tau_dest > 0.1 tau_res`; the count, the worst ratio and its
radius, the count of cells above the 512 K excitation-temperature limit of
the shielding table, and the count of cells in which the He+ channel removes
He+ faster than recombination does are written to the
`output/Oxygen_chemistry.txt` header. The record is informational: the rates
are evaluated everywhere, and where the ordering fails it is the omitted
formation that fails with it, so the transported value stands.

**The He+ ledger (2026-09-06, item B3b-CO2).** D1 consumes a helium ion, and
row (2) of `mol_heh_rows` now loses it: `- k_D1 n_CO n(He+)`, with the CO a
background density of that row because the one-sided model gives CO no
formation and therefore no local equilibrium. The reaction is split across
the two operators, each reading the other's frozen value from the one rate
`rk_D1_Hep_CO`, so the heat of channel 18 is deposited for exactly the events
the species ledger performs. MEASURED on `oxygen_chemistry` at 1000 steps:
the largest share of the total heating taken by channel 18 falls from 0.672
to 0.035, `n(He+)` at the front falls by up to a factor 137, and `log10 Mdot`
stays 9.61. The products need nothing added: the free-atomic-oxygen closure
and `carrier_write_back`'s `nC_free - n_CO` already carry both nuclei, and
the element census closes on carbon and oxygen.

**The thermodynamic ceiling is deleted (2026-09-06, item CEILING-DEL).**
`limit_to_element_budget` now applies conservation and nothing else: the CO
of a cell is bounded by the free carbon and the water family by the free
oxygen, and no thermodynamic bound is imposed on a transported density. Gone
with the clamp are `co_ceiling_applications`, `co_ceiling_cells_hit`,
`co_ceiling_CO_removed` and `co_ceiling_seen`, their attempted-ledger twins,
the three accessors, the `co_ceiling_cells` header line of
`output/Oxygen_chemistry.txt`, the three `co_ceiling_*` keys of
`EXHALE_resolved.out` and the end-of-run "CO ceiling active" line. The domain
record above stands in every one of those places.

**T5.3 is lifted by construction.** The exclusion of section 5.3 is stated on
an observable, ceiling activity during accepted physical integration, that no
longer exists; what replaces it is not an exclusion but a reported domain, so
a run out of domain in some cells is a physical result with that record
attached to it. In the certification the CO items move from validity state
4.1, active unvalidated physics, to 4.2, out-of-domain closure, beside the
H3+ records and informational in the same way. 4.1 is left with no producer
in the code at all and reports "not produced".

### 5.3 Alternative B (adopted): exclusion, with the record as the observable

**T5.3.** Until T5.1's rates, products, energies and domain argument exist,
any configuration in which the CO ceiling is active during **accepted
physical integration** is excluded from validated physical results. The
observable is the cumulative activation record: `co_ceiling_seen(j)` true for
any physical cell, or `co_ceiling_CO_removed > 0`, at any accepted physical
step. The record is:

- accumulated only from accepted steps (A3 and B3a: a rejected trial leaves
  no contribution, so it never writes to the record, and by the B6 semantics
  of review 2 section 5.3 a rejected trial does not invalidate a result);
- carried into the output header and into the B6 validity state of the run;
- separated between the initialization or continuation mode and physical
  integration (A0), so that ceiling activity during relaxation does not
  exclude the physical history that follows a clean handoff;
- kept as a whole-atmosphere diagnostic label even for an excluded run, so
  that the exclusion is auditable after the fact.

**Superseded (2026-09-06, item CEILING-DEL): the ceiling and its record are
deleted, so T5.3 has no observable and the exclusion is lifted. What follows
is the record of the decision as it stood.**

**Decision 7 (2026-09-06): B is adopted.** Exclusion is the specification:
any configuration in which the CO ceiling is active during accepted physical
integration is outside the validated set, with the cumulative record as the
observable and the diagnostic label kept for audit. Alternative A stays open
as the destruction model to be written once the rate literature exists; it is
blocked on data, not on a decision, and the data requirement below states
exactly what is missing. Nothing in A is implemented on the strength of this
acceptance.

**Owner:** (unassigned).

**Acceptance test AT-5. Superseded (2026-09-06, item CEILING-DEL): its four
cases are stated on a ceiling that no longer exists.** What replaces it is a
domain report, not an exclusion: a run states, by ledger family, how many cell
visits fell outside `tau_dest <= 0.1 tau_res`, the worst ratio and where, how
many visits sat above the shielding table's excitation-temperature limit and
in how many He+ + CO led the He+ loss. `carrier_retry` asserts the two
composition rules of that record (counts sum over the families, the worst
ratio is their maximum) and that a refused attempt leaves it unchanged, which
is the deliberate opposite of the old (iii). The superseded text follows as
the record of the decision. (i) A run whose ceiling never fires reports an
empty record and is not excluded. (ii) A run whose ceiling fires once during
accepted physical integration is excluded and says which cells and how many
molecules. (iii) A ceiling firing inside a rejected trial leaves the record
empty after the trial is restored (this is a joint test with A3 and B3a).
(iv) Ceiling activity during initialization does not exclude the subsequent
physical history and is reported separately.

---

## 6. H3+ emission and H2 thermodynamics

### 6.1 H3+ emission

Replaces **D0 C24** (the collider density is clamped to the first table
column, so the cooling does not vanish as `n(H2) -> 0`, and the published
piecewise fits join discontinuously; `h3p_cooling.f90`, `h3p_nonlte_factor`
and `h3p_emission_lte`; both were replaced on 2026-09-06, item B3b-H3, in the
direction this section specifies).
Answers D0 open question 6. This section states **what the target requires of
the published fit definitions**; it does not re-derive the fits.

READ, the module header of `h3p_cooling.f90`: the code uses Miller, Stallard, Tennyson and
Melin (2013, J. Phys. Chem. A 117, 9770) Table 5 for the LTE emission per
molecule, `ln E(T) = sum_n C_n T^n` with `E` in W molecule^-1 sr^-1,
piecewise over 30-300, 300-800, 800-1800 and 1800-5000 K, and Table 6 for
the non-LTE departure factor `s(T, n_H2)` on a grid of `T` and
`log10 n(H2)` from 6 to 14.

**T6.1 (population convention).** The specification states, at one site,
which quantity the published fit is: `E(T)` is the emission of one H3+
molecule per steradian assuming an LTE internal population at `T`, so the
volumetric rate is `Lambda = n(H3+) * 4 pi * E(T) * s`, and `s` is a
departure factor referenced to that same LTE emission. Whether `s` multiplies
`E` (scaled LTE) or is reported as a separate model is decided from that
published definition, following section 5.3 row "H3+ emission": "Choose the
model from the published fit definitions and their internal-population
convention, not by which scaling preserves old cooling." The paper's own
caveat, that Table 6 values are upper limits because of the proton-hopping
rate coefficient, is carried at the code site.

**T6.2 (collider dependence).** `s` depends on the H2 collider density only.
The target forbids the present clamp of `log10 n(H2)` to the table's lowest
column at 6 (`h3p_nonlte_factor`, `h3p_cooling.f90`), because it leaves a finite emission
where there are no colliders. In the low-collider limit, where radiative
decay empties the emitting levels faster than collisions populate them, the
emission per molecule is set by the collisional excitation rate and is
therefore linear in `n(H2)`. **Decision 8 (2026-09-06): below the table's lowest column the target
extrapolates on that limit**, `s ~ s(T, 10^6) * n(H2)/10^6`, and records the
cell as evaluated below the tabulated collider range. Refusing to evaluate
would leave an active coolant undefined in exactly the region where the model
has a known analytic limit, whereas the linear extrapolation is that limit and
carries its own domain statement; the record is informational, not an
exclusion, and the emission it produces vanishes with `n(H2)` as `AT-6a` (i)
requires.

**T6.3 (temperature limits and extrapolation policy).** The four published
segments do not join, and the mismatch is the paper's: the code's own header
records that continued across 300 K the 300-800 K polynomial stands 15 to 43
per cent above the 30-300 K one over 200-300 K, and the 1800-5000 K
polynomial 2.4 per cent below the 800-1800 K one at 1800 K
(the module header of `h3p_cooling.f90`, a value the code site records; not re-measured
here). The target states the policy rather than smoothing it away: each
segment is used only on its published range; below 30 K and above 5000 K the
evaluation is clamped and the cell is marked as evaluated outside the fit's
domain through B6; the discontinuity at a segment boundary is a known
property of the published fits and is REPORTED with its measured size at
the code site and in the records. Because the evaluated cooling function
must be single-valued, the code joins the published segments by a stated
rule (B3b-H3, 2026-09-06: a linear blend in `log_e E` over the 5 percent of
temperature below each join, so that the upper published polynomial is
exact at and above the join and the lower one exact below the ramp, and
every value the paper tabulates is reproduced; at 300 K the paper's own
Table 6 value adjudicates for the upper polynomial, measured 0.002 percent
against 30 percent). That rule is the code's, labeled as such, and is not
attributed to the publication; the alternative of leaving a step in the
evaluated function was rejected because a step is not a property of the
molecule.

**Acceptance test AT-6a.** (i) Low-collider limit: `Lambda / n(H3+) -> 0`
proportional to `n(H2)` as `n(H2) -> 0`. RED today by the clamp. (ii) LTE
join: for `n(H2) >= 1e14 cm^-3`, `Lambda = n(H3+) 4 pi E(T)` within the
tolerance the specification states. (iii) Domain: a cell at 20 K or at
6000 K produces a B6 out-of-domain record. (iv) A cell below the tabulated
collider range is evaluated on the linear extrapolation of decision 8, carries
the informational record that says so, and is not excluded.

### 6.2 H2 thermodynamics

Replaces **D0 C22** (chemistry and the caloric equation of state use
different H2 internal-state models, so one thermodynamic potential does not
generate both).

READ. The EOS evaluates `u_rv` and `c_rv` from the complete observed bound
rovibrational ladder of H2 X^1 Sigma_g^+, 302 levels to 51966 K, from Roueff
et al. (2019, A&A 630, A58) table 2, with the nuclear-spin weights inside the
sum, hence the ortho/para equilibrium mixture
(`caloric_eos.f90`, `build_h2_rovibrational_table` and
`h2_rovibrational_energy_and_heat_capacity`). The chemistry evaluated its own
partition sum `q_rovib_H2` from Huber and Herzberg spectroscopic constants
(`mol_rates.f90`: `D0 = 36118.11 cm^-1`, `omega_e = 4401.213`,
`omega_e x_e = 121.336`, `B_e = 60.853`, `alpha_e = 3.062`), tabulated on
100-20000 K and interpolated. **That second source was removed on 2026-09-06
(item B3b-H2Q)** in the direction T6.4 specifies: one
`h2_partition_function` in `caloric_eos.f90` now serves both.

**T6.4 (one potential).** One partition function `Q_H2(T)`, from one level
set, with one zero of energy at `v = 0, J = 0`, generates both:

```text
u_rv(T) = k_B T^2 d ln Q_H2 / dT              (the EOS)
ln K_eq(T) = ... ln Q_H2(T) ... - D0 h c / (k_B T)   (the chemistry)
```

**Decision 9 (2026-09-06): that one level set is the Roueff et al. (2019)
observed bound rovibrational ladder**, the set the equation of state already
uses, and the Huber and Herzberg constants of `mol_rates` stop being a second
source of `Q_H2`. `K_eq` moves by the amount the two models differ and that
movement is the physics, not a regression to be suppressed; it is measured and
reported with the change. A test that pins the two present models together is
kept only while the transition is in progress and is not the ownership model.
With one source, the forward and reverse rates, the caloric equation of state
and the reaction enthalpy of section 1's reservoir are thermodynamically
consistent by construction rather than by coincidence. The zero of the ladder
is the same zero the reservoir of `T1.3` uses for `eps(H2) = -D0`, which is what
makes section 1's non-overlap rule verifiable for H2.

**Acceptance test AT-6b.** (i) Van 't Hoff: `-d ln K_eq / d(1/T)` computed
from the shared `Q_H2` reproduces the reaction enthalpy assembled from the
section 1 reservoir plus the EOS heat capacities, to a stated tolerance,
over 100-5000 K. (ii) `c_v(H2)` from the shared
`Q_H2`, which is the Roueff ladder of decision 9, is reported against both of
today's models and the movement of `K_eq` is reported with it, rather than
assumed zero: the code site already records that the observed ladder exceeds a
rigid rotor plus harmonic oscillator by 1.4 per cent at 200-900 K, 2.3 per
cent at 1300 K and 7.3 per cent at 3000 K (the module header of
`caloric_eos.f90`, code site records). (iii) A run with no molecules reproduces the constant `gamma_ad`
arithmetic exactly, the structural branch the EOS already carries
(`caloric_state_from_composition`, `caloric_eos.f90`).

**Owner:** (unassigned).

---

## 7. Rollback contract

Rev 2 section 5.1: the B1 checkpoint must cover **every** state-changing
operation before the final adoption boundary, viscosity, conduction and the
optional Shapiro filter included, or state a narrower interim contract.

### 7.1 The enumerated step

READ from `src/EXHALE_main.f90` at 10:29 KST on 2026-09-06, in the order the
main program executes them. The **adoption boundary** is the
`update_map_end_step` call of the marching loop: from there onward the loop
reads the state it has adopted and computes diagnostics and control flags.
Everything above it is inside the attempted step and must be restorable. The
row numbers below are the operation tags `attempted_step.f90` declares
(`as_op_dt` through `as_op_boundary`), which is the durable index: that module
and this table number the same fourteen operations, and a rejection names its
operation with them.

| # | tag | operation | state it changes |
|---|---|---|---|
| 0 | `as_op_dt` | `eval_dt(W,dt,dt_loc)`, update-map scaling of `dt`, `dt_loc` | `dt`, `dt_loc` (main) |
| 1 | `as_op_checkpoint` | `u_old = u` | the checkpoint anchor; the retry loop restarts every attempt from it |
| 2 | `as_op_hydro` | `retry_step` loop: three RK stages, `Apply_BC`, `positivity_limited_fluxes`, `species_advection_stage` | `u`, `u1`, `u2`, `WL`, `WR`, `dF`, `S` (main); `n_faces_flux_positivity_limited` and `n_faces_flux_positivity_limited_accepted` (`RK_rhs.f90`), `n_dt_halve`, `n_steps_dt_halved` (`EXHALE_main.f90`), and `dt`, `dt_loc` on a halving; since B4-1 also the advected element mass fractions `Yetr`, `Yetr_n` (`binary_element_diffusion.f90`) and the operator's counters `n_species_faces_bounded`, `species_face_excursion`, `species_face_identity_residual` (`species_face_flux.f90`). `Yetr` needs no checkpoint entry: `species_advection_begin_step` rebuilds it from `f_sp` at the top of every attempt, so an attempt retaken after a restore starts from the restored composition |
| 3 | `as_op_primitives` | `U_to_W`, `get_species_densities`, `comp_T_from_p` | `W`, `rho`, `v`, `p`, `nhi`...`nheiTR`, `nm`, `ne`, `n_tot`, `T` (main) |
| 4 | `as_op_diffusion` | `species_advection_project` and `element_diffusion_step` (`he_diffusion`) | `f_sp` (main), which row 4 now writes twice: once with the element totals the Runge-Kutta stages advected and once with the diffusive step's own solution; `he_fraction_over_one`, `he_fraction_under_zero`, `he_fraction_newton_steps`, `he_fraction_newton_resid` and `trace_ratio_under_zero`, all declared in the head of `binary_element_diffusion.f90` |
| 5 | `as_op_carriers` | `photochemical_transport_step` | `f_sp` carrier rows (main); the module-save state declared in the head of `diffusive_photochemistry.f90`: `n_carrier`, `carrier_solved`, the `cbg_*` saved backgrounds, the `cph_klw`, `cph_jh2o`, `cph_joh` photolysis arrays, `pct_newton_resid`, `pct_newton_resid_before_limit`, `pct_cell_constrained`, `pct_verdict`, `pct_worst_limit`, `pct_Dco`, `row_terms`, `col_scale_car`, `row_scale_car`, `headroom_car`. `save_carrier_module_state` and `restore_carrier_module_state` are the one enumeration of that set |
| 6 | `as_op_excited_H` | `excited_H_update` (`use_excited_H`) | globals `gph_balmer_HI`, `heat_balmer`, `Jlya_arr`, `n2s_arr`, `n2p_arr`, `Sproton_arr`, `Hpe_arr`, `Hdx_arr` (declared in `parameters.f90`); `jlya_rt_loaded`, `jlya_rt_grid`, `Tdiag`, `nhidiag`, `nediag`, `nhiidiag`, `taulya` (head of `excited_hydrogen.f90`); `jint_arr`, `jstar_arr` (head of `lya_rt.f90`) |
| 7 | `as_op_ioniz_eq` | `ioniz_eq(T,rho,f_sp,heat,cool,eta)` | `f_sp`, `heat`, `cool`, `eta` **and `rho`** (the `C2` write: the `n_io` dummy of `ioniz_eq` was `intent(inout)` and was overwritten from `calc_rho`, both in `ionization_equilibrium.f90`, until item B3c on 2026-09-06 made it `intent(in)`); module state declared in the head of that file: `nmol_eq`, `NH2_col_lw`, `f_shield_lw`, `k_lw_diss`, `p_lw_single`, `p_lw_absorbed`, `tr_lines_lw`, `P_H2_eq`, `nox_eq`, `n_o1d_eq`, `NH2O_col`, `NOH_col`, `j_h2o_fuv`, `j_oh_fuv`, `tau_fuv`, `heat_fuv`, `heat_chem`, `bg_cell`, `bg_cell_adopted`, `bg_cell_best`, `ieq_sweep_state_kind`, `ieq_marching_ledger`, `ieq_steady_iterate_ledger`, `ieq_steady_candidate_ledger`, `ieq_nonroot_streak`, `ieq_acc_nprint` |
| 8 | `as_op_composition` | `get_species_densities`, `comp_p_from_T`, `W` assembly, `W_to_U` | `ne`, `n_tot`, `p`, `W`, `u` (main). **This was the `C1` projection**; section 1 removes it as an independent update and item B3c did so on 2026-09-06. The tag survives and now names the point at which the composition and the energy of the coupled source step have become one state (`attempted_step.f90` says so at the declaration) |
| 9 | `as_op_energy` | `solve_energy_semi_implicit` (or the explicit `u(3,:)` update beside it) | `u` (main); `n_energy_floor_hits`, `energy_floor_first_step`, `energy_floor_last_step`, `energy_floor_cell_hits`, declared in the head of `energy_semi_implicit.f90` |
| 10 | `as_op_bc` | `Apply_BC(u)` | ghost cells of `u` |
| 11 | `as_op_conduction` | `viscous_conduction_step` with its `U_to_W`, `comp_T_from_p`, `Apply_BC` (`transport_active()`) | `W`, `rho`, `v`, `p`, `T`, `u` (main); `n_conduction_floor_hits`, `conduction_floor_first_step`, `conduction_floor_last_step`, `conduction_floor_cell_hits`, declared in the head of `viscous_conduction.f90` |
| 12 | `as_op_shapiro` | `shapiro_filter(u)`, `Apply_BC(u)` (`shapiro_eps > 0` and on the `shapiro_every` cycle) | `u` only; the filter holds no module state (`shapiro_filter`, `Apply_BC.f90`) |
| -- | `as_op_boundary` | **adoption boundary** | |
| 13 | (after `as_op_boundary`) | post-adoption reads and control state: `U_to_W`, `check_base_inflow_is_subsonic`, `get_species_densities`, `comp_T_from_p`, `du` and arming flags, residual diagnostics, and the two "Staged secondary ionization" flips further down the loop | control and diagnostic state, not the physical state; `sec_ion_active` (declared in `parameters.f90`) and `sec_ion_armed_step` are mode state and belong to A0, not to the rollback |

The update-map snapshots `u_umA`, `u_umB`, `u_umC` and `u_umD` mark the same
operator-split seams and are the existing instrumentation of this boundary.
The comment above `u_umA` in `EXHALE_main.f90`
records why a snapshot must not be placed inside an RK stage: a statement
there changes how the compiler contracts the update into an FMA, which was
measured at one ulp per step growing to 1e-11 over a converging run
(code site records, 2026-09-03; not re-measured here).

### 7.2 The contract

**T7.1.** The checkpoint saves, before operation 1, every item in the "state
it changes" column of rows 1 through 12, and restores all of them on a
rejected attempt. `u_old` alone is not a checkpoint: rows 4 through 9 and 11
mutate module-level saved arrays that `u = u_old` does not touch, so a retry
today re-enters the step with, for example, a partially updated
`bg_cell_adopted`, `ieq_nonroot_streak` and the floor-hit counters of rows 9
and 11.

**T7.2.** Counters are partitioned into three classes, and the checkpoint
treats them differently:

- **physical accumulations** (the CO record, reaction and heating budgets,
  elapsed physical time): restored on rejection, so a rejected trial leaves
  no contribution, which is the B6 semantics of review 2 section 5.3;
- **attempt statistics** (attempted steps, retries, halvings): deliberately
  **not** restored, since they count attempts;
- **diagnostic extrema** (`he_fraction_over_one`, `trace_ratio_under_zero`,
  floor-hit counters): **decision 10 (2026-09-06): kept, not discarded.** An
  extremum reached only inside a rejected trial is retained as an attempt
  statistic, labeled as such, and never becomes a property of the adopted
  state, so it is informational and cannot invalidate the accepted result. The
  reason to keep it is that it is the only trace of what the rejected
  direction did, and the reason it cannot invalidate is A0's rule that a
  rejected trial contributes nothing to the physical history. Each counter
  states which of the two families it belongs to at its declaration site.

**T7.3 (the narrower interim contract, if B3a lands before B3c).** If the
controller is first built as a restoration mechanism against the present
update, the contract it satisfies is stated as: rows 1 through 12 restored,
physical accumulations restored, attempt statistics excluded from restoration.
That demonstration certifies restoration only. It does not certify the `C1`
composition reset that row 8 still contains, and no energy tolerance is
loosened to make it pass. This is rev 2 section 5.1 as an explicit clause.

**T7.4 (the validity states of the enumerated step).** The five B6 categories
are those of B1a section 4: **(1)** active unvalidated physics, **(2)** an
out-of-domain closure, **(3)** a rejected numerical trial with no adopted
contribution, **(4)** an unbudgeted accepted correction, **(5)** a specified
external reservoir, which is informational. The table gives, for each operation
of 7.1, the categories it can produce and what happens to the record on a
restore. Every code site named is READ; the line numbers carry the drift note
of the provenance section.

| # | operation | validity states it can produce | survives a restore? |
|---|---|---|---|
| 0 | `eval_dt` and the update-map scaling | none | not applicable |
| 1 | `u_old = u` | none | the anchor itself is the restore |
| 2 | RK stages, `Apply_BC`, positivity-limited fluxes | (4) reconstruction faces scaled toward the cell average and ghost extrapolation dropped to zero gradient, in an accepted step; (3) the same repairs inside an attempt the `dt` bisection discards | no for the state records; the attempt counts stay. The code already separates the two (`calls`, `faces`, `faces in accepted steps`, `write_run_counter_report` in `EXHALE_main.f90`) and that separation is the model the rest of the table follows |
| 3 | `U_to_W`, `get_species_densities`, `comp_T_from_p` | none of its own; an inadmissible state here is a rejection, not a validity state | not applicable |
| 4 | `element_diffusion_step` | (4) the helium fraction range clip and `trace_ratio_under_zero`; (2) Blanc's law and the Chapman-Enskog first approximation, both stated at the code site and neither checked; (5) the inner Dirichlet reservoir of the operator | no. The clip counters are diagnostic extrema and follow the T7.2 policy of decision 10: kept as attempt statistics, never a property of the adopted state; the (5) record is a property of the configuration and is not produced by the step at all |
| 5 | `photochemical_transport_step` | the CO ceiling was category (4) by decision 12 and was **deleted on 2026-09-06 (item CEILING-DEL)**, so no correction of that kind is left in this operator; what stands in its place is the destruction model's domain record, which is a measurement and not a correction. Remaining: (4) the carrier element rescale of `limit_to_element_budget`; (3) a carrier trial above the headroom, which is refused rather than clamped; (2) the H2 self-shielding table clamps | no. The refusal count is an attempt statistic and stays; the domain record is deliberately not restored, since whether a state lay inside the domain has an answer whether or not the step that read it was accepted |
| 6 | `excited_H_update` | (2) the damping-wing slab closure evaluated where `(a tau)^(1/3) > 10` fails, and the Sobolev channel outside its own domain; (1) the H(n=2) lag, which is a closure evaluated one outer pass behind the state it closes and whose lag error is not measured | no; both are properties of the attempted state |
| 7 | `ioniz_eq` | (1) the oxygen ledger gap (section 1.6), the H3+ model domain, the H2 caloric and chemistry state mismatch, the metal Voronov transcription note; (2) every clamped closure of B1a section 4.2 that the sweep evaluates, the CHIANTI coronal cutoff, the Ca II and C I range clamps, the infrared band clamps, the Penning and Oklopcic-Hirata ranges, the secondary-ionization extrapolation, the Gaunt table clamp; (3) steady candidate probes, which are already tagged `ieq_state_steady_candidate`; (4) the simplex clamp accepted as class 3, the class-4 relaxation amnesty, and the `rho` rewrite until section 2 lands | the (1) states are static properties of the configuration and are not produced by a step, so a restore does not touch them; every (2) and (4) record produced inside a rejected attempt is discarded with the attempt; the probe tags are attempt statistics and stay |
| 8 | the composition projection, `comp_p_from_T` and `W_to_U` | (4) the `C1` projection itself, in every configuration, with no run record today | no. Section 1 removes the producer; until then its activity refuses certification of a physical step (A2 section 5) |
| 9 | `solve_energy_semi_implicit` or the explicit update | none since B2 (2026-09-06): the solve is residual-controlled and a floor activation is a rejected attempt (status `FLOOR`), recorded in the attempts ledger, not an accepted state; the explicit `u(3,:)` update still has no returned-state test, a (1) state of that path | no |
| 10 | `Apply_BC(u)` | (5) the base `(p, s)` reservoir, the base composition handoff and the base infrared field; (2) a supersonic base inflow, which over-specifies the characteristic boundary and is diagnosed only today; (4) the outer ghost extrapolation dropped to zero gradient | the (5) records are configuration properties and always survive; the (2) and (4) records of a discarded attempt do not |
| 11 | `viscous_conduction_step` | none since B2b (2026-09-06): a floor activation is a rejected attempt (status `FLOOR`), no clamped state is returned; the base Dirichlet anchor is a (5) reservoir exchange stated at the code site | no |
| 12 | `shapiro_filter` | (4) unconditionally: the filter alters the adopted state with no source term behind it, so an active filter is an unbudgeted accepted correction unless the energy it removes is budgeted and reported. Default off (`shapiro_eps = -1`) | no; and an active filter is a validity state of every step it runs on, not only of the step that is restored |
| 13 | post-adoption reads and control state | none: it changes control and diagnostic state, not the physical state. The staged secondary-ionization flip is mode state and belongs to A0 | not part of the rollback |

**The restore rule, stated once.** A validity state is a property of a state.
The states of an attempt that is not adopted are discarded with the attempt,
which is the review 2 section 5.3 semantics that a rejected trial leaving no
contribution does not invalidate a result. The only class that survives a
restore is the attempt record itself, category (3), and it is recorded as an
attempt statistic and never as a property of the adopted state. Category (1)
and the configuration half of category (5) are not produced by a step at all:
they follow from the flags and hold for the whole run, so a restore neither
sets nor clears them.

**Owner:** (unassigned).

**Acceptance test AT-7.** (i) An injected failure after row 5 and after row 9
restores every item of the table bit for bit, verified by hashing the saved
set before and after. (ii) The elapsed physical time after a rejected attempt
equals the time before it. (iii) No physical accumulation, and in particular
no reaction budget, contains a contribution from a rejected attempt. (iv) The attempt statistics do contain it.

---

## 8. Derived products: the `_adv` rule

Replaces **D0 C8** (the post-process temperature is solved with a different
channel set from the hydrodynamic temperature: no H3+, no molecular infrared
bands, molecular ions absent from `n_e`, and on the heating side photoheating,
the He recombination coupling and the He(2^3S)+H Penning heat only, so no
Balmer, no Lyman-Werner, no molecular chemical and no FUV heat;
`T_equation` (`T_equation.f90`), and the heating and electron-density assembly of `post_process_adv`) and **D0 C20**
(the `_adv` profiles the transit tools read are reconstructed on a
molecule-free background, and their molecular columns are the equilibrium
sweep's values; the module header of `post_process_adv.f90`, the limitation
printed into the file by `write_output`). D0 open questions 11 and 12 are the two
readings this section closes.

**T8.1 (one assembly, or no corrected product).** A derived product is
assembled from the **same source assembly as the main solution**: the same
heating and cooling channels, the same electron and particle sums, the same
species set, on the active-equation set of T1.7 for that configuration. One
definition of each process, with at most two entry points into the same code,
and a test that pins the two entry points together on a shared state.

Where a channel that is active in the run cannot be evaluated in the
post-process, the corrected product for the affected quantity is **not
written**: it is refused, the file states which quantity was refused and which
channel was missing, and the uncorrected profile keeps its own name.
**Decision 14 (2026-09-06): the suppression is by quantity, not by file.** The
columns that do not depend on the missing channel are still written as
corrected, each column carries whether it is corrected or uncorrected, and a
column that depends on a refused one is refused with it, which is why a missing
cooling channel refuses the corrected temperature and everything derived from
it. A file that mixed the two silently is the defect `C20` names; a file that
labels each column is the target. Writing a
corrected temperature whose balance omits channels the run solved is the defect
`C8` names, and it is worse than writing none, because the transit tools cannot
see the omission. A product that is refused is not a run failure; it is a
missing column with a stated reason.

**T8.2 (the product names the state it was made from).** Every derived product
carries, in its header, the A0 `mode=` field of the state it was made from
(`init`, `phys` or `certified`) and, for a physical state, that state's
`t_phys`; plus the state's identity (the `solution_id` or step count already in
the coupling header) and the B6 validity state of that state, including the
static entries of section 1.6 and section 5. A0 section 5 states the same
requirement from the other side: "Every derived product (`_adv`, breakdowns,
transit inputs) carries the same mode field as the state it was made from."
A product made from an initialization snapshot is labeled as such and is never
presented as a result; a product made from an uncertified physical state is
labeled a transient, which is a valid thing to be (review 2, F1).

**Acceptance tests.** `AT-8g`: on one state, the post-process assembly and the
main assembly are evaluated for the same cell and the same channels and agree
to round-off for every channel the configuration activates; a configuration
with H3+ or the molecular infrared bands active makes the present code RED by
the whole missing channel. A channel that cannot be evaluated produces the
labeled refusal of decision 14: the corrected temperature and everything
derived from it are refused with the channel named, and the columns that do
not depend on it are still written as corrected and are labeled as such. `AT-8h`: an `_adv` product made from an `init`
snapshot carries `mode=init` and no `t_phys`; made from a physical state it
carries `mode=phys` with the state's `t_phys`; a product whose source state
carries a static exclusion (the oxygen energy ledger; the CO ceiling until its
deletion on 2026-09-06, item CEILING-DEL) carries the same exclusion, and a
tool reading the file can tell.

**Owner:** (unassigned).

---

## 9. Verification matrix

Rev 3 plus review 2 section 7. Units, tolerances and floors are set in each
implementation brief; this table fixes what is tested and against what.

| # | target | test | acceptance | reference |
|---|---|---|---|---|
| AT-1a | section 1 | closed reacting cell, `Q_ext = 0`, across a dissociation and an ionization front | `u_th + u_form` constant to round-off per accepted step | thermodynamic identity |
| AT-1b | section 1 | isolated photon event, with and without secondary ionization | each recipient assigned once, shares sum to `h nu` | Schulik and Booth (2023) section 3.4, eqs. 53-55 |
| AT-1c | section 1 | single recombination event | radiated and thermal shares as the adopted cooling function implies, defined at one site | |
| AT-1d | section 1.6 | closed cell with the oxygen network; an O4 and O6 pair against the direct route; a single associative He(2^3S) event; the static exclusion | ledger closes to round-off; the pair deposits what the direct route deposits; 8.1 eV deposited; the exclusion reported and certification refused | code sites of `C25` |
| AT-2 | section 2 | update map `chem` column in the mass row; source inspection of the stationary assembly; species mass sum | zero chemistry source in the mass row; no chemical mass source in the stationary equations; mass sum within `tol_mass` | |
| AT-3 | section 3 | two-cell radiation exchange at `f = 0` and `f = 1`, then with the bottom absorber | material plus escaped energy closes to round-off | Huang et al. (2023) section 2.7; Salz et al. (2016) sections 3.1, 4 |
| AT-4a | section 4 | static column with a trace metal | `sum_s rho_s w_s` at the stated order | Koskinen et al. (2022) eq. B10 |
| AT-4b | section 4 | multi-stage ion column | `sum_s z_s n_s w_s` at the stated order | derived closure (J0) |
| AT-4c | section 4 | `dmeff` in three plasma limits | 2.9715, 2.4715, 1.6477 at `m_He/m_H = 3.9715` | code site header |
| AT-4d | section 4 | diffusion-only column through a transient | element totals to round-off | conservative form |
| AT-5 (superseded by item CEILING-DEL, 2026-09-06) | section 5 | ceiling inactive, active, active-in-a-rejected-trial, active-in-initialization | exclusion only for accepted physical activity; record auditable. No longer testable: the ceiling is deleted and the CO destruction domain is reported, not excluded | rev 3 item 9; review 2 section 5.3 |
| AT-6a | section 6.1 | `n(H2) -> 0`; `n(H2) >= 1e14`; out-of-range `T` | linear vanishing; LTE join; B6 out-of-domain record | Miller et al. (2013) Tables 5, 6 |
| AT-6b | section 6.2 | van 't Hoff on the shared `Q_H2`; `c_v(H2)`; molecule-free run | consistency to the stated tolerance; difference reported; constant `gamma_ad` exact | Roueff et al. (2019) table 2 |
| AT-7 | section 7 | injected failure after two different operations | complete restoration; correct elapsed time; no rejected contribution in a physical accumulation; attempt statistics keep it | rev 2 section 5.1 |
| AT-8a | B3, physical mode | temporal refinement: successively halved steps to one fixed end time | observed temporal order matches the claimed one; splitting error controlled | rev 3; D0 section 6.6 records that no such measurement exists |
| AT-8b | B3c | splitting-order check: operator order reversed at fixed end time | difference at the order the splitting claims | |
| AT-8c | B2 | returned-temperature residual on a falling cooling branch, a strongly temperature-dependent heating, an unavailable bracket, a floor encounter | final residual tested; a floor is a specified reservoir with budgeted heat or a defined failure | |
| AT-8d | B5 | species stationary residual for every active transported balance | evaluated on the independent space of T1.7 for that configuration, same source assembly as B3c, same face fluxes as B4 | B1a section 2.2 |
| AT-8d (i) | B5 | the row registry per configuration: one row and one unknown per transported balance, the band half-width `3 nvar - 1` and the coloring stride `6 nvar - 1` | the counts of T1.7's table for configurations 3, 4, 5, 8 and the oxygen-and-proton set; the two geometries the code carried as constants, `nvar = 3 -> (8, 17)` and `nvar = 4 -> (11, 23)`, reproduced by the formula | `src/tests/steady_species_rows/`, driver rows |
| AT-8d (ii) | B5 | Jacobian action of the system the solve carries: `J r` against `[F(Y + eps r) - F(Y)]/eps` of the same map, species rows measured on their own | hydrodynamic rows to 1e-4, species rows to 1e-2, and at least one Jacobian color resolved; a wrong band layout or stride is out by O(1) | `EXHALE_SPECIES_JAC_TEST=1`, log rows of the same suite |
| AT-8d (iii) | B5 | the completion flag reads the species rows of the system it solved | a solve whose species row is above `cert_tol_carrier` may not end `info = 0` | same suite, log rows |
| AT-8g | section 8 | post-process and main source assembly on one state | every active channel agrees to round-off, or the corrected product is refused with the channel named | D0 C8 |
| AT-8h | section 8 | header of a derived product made from an `init`, a `phys` and an excluded state | mode, `t_phys` and validity state carried through | A0 section 5 |
| AT-8e | all | invariants: element nucleus totals, charge neutrality, total mass, the stated current condition, material energy | closed at every accepted step | |
| AT-8f | B3c | endpoint source outputs beside a substepped or multistage energy budget | each at the correct state and time level, one definition per process | review 2 section 5.2 |

### 9.1 The rows each configuration carries

The eight configurations are those of T1.7. "n/a" means recorded as not
applicable, never skipped: A2 reports a status for every row of the matrix.

| # | configuration | rows required | rows recorded n/a | validity states that hold unconditionally today |
|---|---|---|---|---|
| 1 | pure H/He | AT-1a, AT-1b, AT-1c, AT-2, AT-6b (iii), AT-7, AT-8a, AT-8b, AT-8c, AT-8e, AT-8f, AT-8g, AT-8h | AT-3 (no explicit H(n=2) population: `Stellar Teff` and `Stellar radius` are not both given, so `use_excited_H` is false), AT-4a to AT-4d, AT-5, AT-6a, AT-6b (i) and (ii) | the `C1` projection and the `C2` density rewrite; the Lyman-alpha collisional excitation charged as escaping radiation with no escape probability and no explicit 2p population (`C7`), which is the configuration where a multiplier is not even available |
| 2 | He 2^3S and metals | the rows of 1, plus AT-3 (H(n=2) active) | AT-4a to AT-4d, AT-5, AT-6a, AT-6b (i) and (ii) | the above, plus the metal Voronov transcription note, the CHIANTI fit range with its coronal cutoff, the Oklopcic and Hirata recombination range and the Penning fit range |
| 3 | hot Uranus molecular | the rows of 1, plus AT-6a and AT-6b (i), (ii), and AT-3 where the run gives both stellar keys | AT-4a to AT-4d, AT-5, AT-8d (no transported balance) | the above, plus the H3+ collider clamp and fit joins, the H2 caloric and chemistry state mismatch, the H2 self-shielding table clamps, and the infrared band clamps where `Molecular IR bands` is on |
| 4 | carriers | the rows of 3, plus AT-8d for the H2 balance | AT-4a to AT-4d, AT-5 | the states of 3, plus the carrier element rescale |
| 5 | oxygen | the rows of 4, plus AT-1d (section 1.6) and AT-8d for the four transported balances; AT-5 is superseded (item CEILING-DEL) | AT-4a to AT-4d | the states of 4, plus the static oxygen ledger exclusion (T1.10) and the CO destruction domain record, which replaced the `AT-5` observable |
| 6 | diffusion | the rows of 3, plus AT-4b, AT-4c, AT-4d, and AT-8d for the elemental balance | AT-4a (no trace-metal arm: `He_metal_diffusion` is not set), AT-5 | the states of 3, plus Blanc's law and the Chapman-Enskog first approximation, both stated and unchecked, and the advective form of the elemental transport (`C28`) |
| 7 | profile | AT-1a, AT-1b, AT-1c, AT-2, AT-3, AT-4a, AT-4b, AT-4c, AT-4d, AT-6b (iii), AT-7, AT-8a to AT-8h | AT-5, AT-6a, AT-6b (i) and (ii) | the states of 2, plus `C28` and `C29`, the trace-metal arm with no counter-flux, which this configuration owns |
| 8 | H+ transport | the rows of 4, plus AT-8d for the H2 and H+ balances | AT-4a to AT-4d, AT-5 unless oxygen is added | the states of 4, plus `C18`, the two base ghosts re-solving H and H+ locally. `C19` is retired: since 2026-09-07 the configuration has a steady solver, the proton being an unknown of it (`Solver: Newton` with `Coupled carrier solve`, `input_read.f90`), and certification is no longer only through a marching solution |

Every configuration additionally carries the always-required items of B1a
section 5: A0 and A0-impl in physical mode, A1, A2, the exact output state of
B3c's assembly, B6 clear or the run excluded, and C1 and C3 for the input
field.

### 9.2 Boundary equations by cell class

What is **imposed** is an equation of the system; what is **checked** is a
diagnostic that may refuse a state but never alters it. The present code is
READ from B1a section 3; the target column states what changes.

| cell class | cells | imposed | checked, never imposed | target |
|---|---|---|---|---|
| base reservoir ghosts | `1-Ng .. 0` | the characteristic face condition at `r_edg(0)`: the LODI `C^-` compatibility relation plus the `(p, s)` reservoir, with the number of conditions set by the eigenvalue signs, so `0 < v < c` gives two reservoir conditions and one interior relation (`base_boundary.f90`, module header and `characteristic_base_face_state`); the ghosts are volume averages of the hydrostatic isentrope through the face state, not copies of it (`base_ghost_averages`, called from `Apply_BC_W`) | that the base inflow is subsonic, since a supersonic inflow over-specifies the boundary (`check_base_inflow_is_subsonic`, `Apply_BC.f90`) | the subsonic condition becomes a stated domain of the boundary closure: outside it the state carries a category (2) validity record and certification is refused, in place of the present diagnosis only (`C33`) |
| base reservoir ghosts, composition | same | the H2 Dirichlet inflow composition, applied as an advective inflow term when the base face flux is inflowing and a handoff states the partition (`diffusive_photochemistry.f90`, the `base_dirichlet` note and `carrier_source`); no diffusive flux at face 0; the diffusion operator's own inner Dirichlet reservoir (`element_diffusion_step`, `binary_element_diffusion.f90`) | the element and charge constraints of the ghost composition | the proton deliberately gets no Dirichlet value and the imposed partition applies only for `j >= 1`, so the two ghosts re-solve their ionization locally (`C18`). That local equilibrium is itself a closure and the target requires it to be stated as one, with its residual evaluated by A2 like any other closure, instead of being an unstated difference between the ghosts and the first physical cell |
| physical cells | `1 .. N` | the hydrodynamic rows, the local rows of T1.7, and the transported balances of the configuration; in the stationary solve `set_base_fix` may additionally anchor the first `nfix` cells by replacing their rows with `Y - Yfix` (`set_base_fix` and `apply_base_fix`, `steady_newton.f90`) | the constraints of section 2 and of T1.7: element nucleus totals, charge neutrality, the species mass sum against the supplied density | an anchored row is an imposed condition and is declared as one in the certification record, so that a certified state says how many of its cells carried an equation of the physics and how many carried an anchor |
| outer ghosts | `N+1 .. N+Ng` | free outflow as a zero-gradient copy, with linear extrapolation added under WENO3 as an accuracy device (`free_outflow_ghost`, `Apply_BC.f90`) | that the extrapolation left `rho > 0, p > 0`; where it did not, the drop back to zero gradient is counted | the drop is a category (4) record of the state it was applied to, as it already is by count; the target adds that a state whose outer boundary is not outflowing is out of the closure's domain, since the free-outflow condition imposes no incoming characteristic |

---

## 10. Decisions taken with the acceptance

Every item below was open in the draft and is **decided as recommended on
2026-09-06**. The decision is stated first; the question it answers follows as
the rationale, in the words the draft used, so that what was traded away stays
readable. The target sections named are written unconditionally on these
decisions.

1. **Reservoir zero for the metals. Decided 2026-09-06: `eps_s` for a metal
   ion is the sum of its ionization potentials from the neutral ground state,
   the same declared reference every other species uses.** The question was
   whether it should instead be measured from the element's dominant state at
   the base, which would make the reservoir smaller but tie the zero to a
   configuration. A reference that moves with the state is not a reference.
   (Section 1.2.)
2. **Electron mass convention. Decided: the species mass table carries the
   electron mass with its ion, so `m_s` for stage `k` is the neutral atomic
   mass and free electrons carry none, and `tol_mass` follows as
   `n_e m_e / rho`.** The question was whether to carry it separately, given
   that the choice sets `tol_mass` and changes the meaning of `calc_rho`.
   With the mass on the ion an ionization moves no mass between species.
   (Section 2, T2.2.)
3. **`ioniz_eq` interface. Decided: the density argument becomes
   `intent(in)`.** The question was whether to confirm a change to a routine
   called from the marching loop, the steady residual and the post-process
   before it is coded. Each caller keeps the density it supplied, and the
   reconstructed density stays available as a diagnostic. (Section 2.)
4. **The residual H I collisional-excitation cooling. Decided: keep the
   Cen-type term as it stands for the higher levels, with its `1s -> 2p` share
   subtracted using the same rate coefficient the 2p balance uses.** The
   question was whether to replace it by a fit that excludes Lyman-alpha
   explicitly. No such published fit was identified, so the subtraction is
   stated as the target's own approximation at the code site. (Section 3.3.)
5. **Trace metals in the diffusion matrix. Decided: solve them jointly with
   the H/He binary in one inversion, which satisfies (B10) by construction.**
   The question was whether to keep the sequential solve with a measured and
   bounded residual instead. The sequential solve survives only as an interim
   carrying that measurement and a stated domain. (Section 4.4, T4.4.)
6. **H+ diffusion. Decided: H+ enters the (B9) set with its charge in the
   target, and the present `w_s = 0` treatment is an interim that reports the
   `Pe` diagnostic and claims nothing about ion separation.** The question was
   whether to make the change now, which moves every ionized-region profile,
   or to defer it until B4. The closure (J0) sums over all charged species and
   is not consistent with holding the most abundant ion fixed, so the target
   cannot be written the other way. (Section 4.4, T4.6.)
7. **CO. Decided: adopt exclusion now (alternative B) and hold the
   destruction model (alternative A) open.** The question was whether to wait
   for the rate compilation before deciding. A is blocked on data, not on a
   decision: its rate coefficients, their temperature ranges and the
   photodissociation band data are still to be obtained, and nothing of A is
   implemented on this acceptance. (Section 5.)
8. **H3+ below the table. Decided: extrapolate on the linear low-collider
   limit and record the cell as evaluated below the tabulated collider
   range.** The question was whether to refuse to evaluate and mark the cell
   out of domain instead. Refusal would leave an active coolant undefined
   where the model has a known analytic limit. (Section 6.1, T6.2.)
9. **The shared H2 partition function. Decided: the Roueff et al. (2019)
   observed ladder is the single source, and `K_eq` moves by the amount the
   two present models differ.** The question was whether to keep both and pin
   them with a test, which section 5.3 of the recommendations calls useful in
   transition but not the permanent ownership model. The movement is measured
   and reported with the change. (Section 6.2, T6.4.)
10. **Diagnostic extrema across a rollback. Decided: an extremum reached only
    inside a rejected trial is kept as an attempt statistic, labeled as such,
    and never invalidates the accepted result.** The question was whether to
    discard it with the trial. It is the only trace of what the rejected
    direction did, and A0's rule already keeps a rejected trial out of the
    physical history. (Section 7.2, T7.2.)
11. **The identity row of a transported species. Decided: the row
    `x_i - x_i^fix` stays in the local system.** The question was whether to
    remove the row and its column and substitute the imposed value into the
    remaining rows, which gives a smaller and better conditioned local system
    but a dimension that changes with the transport flags. One layout for
    every configuration was preferred. (Section 1.4, T1.7.)
12. **The validity category of the CO ceiling. Decided: category (4), an
    unbudgeted accepted correction, which invalidates the state it fires in;
    the category (1) label is kept beside it for audit.** The question was
    whether it is instead active unvalidated physics, which would label rather
    than exclude a ceiling-active run. It moves carbon, oxygen and chemical
    energy with no rate and no source, and review 2 section 5.3 requires an
    unbudgeted accepted correction to invalidate. The decision is the same one
    as item 7, seen from the validity side. (Section 7.2, T7.4.)
    **Moot since 2026-09-06 (item CEILING-DEL): the ceiling is deleted, so no
    state is categorized under it. The CO items are now reported under
    validity state 4.2, out-of-domain closure, beside the H3+ records, and 4.2
    is informational and does not invalidate. Validity state 4.1 is left with
    no producer and reports "not produced".**
13. **The reservoir of an eliminated excited species. Decided: O(1D) stays
    eliminated by its fast closure, and the closure transfers `eps(O(1D))` to
    the products at its own rate.** The question was whether to promote it to
    an explicit species with its own `eps` and row. Its local steady state is
    exact, so the promotion would add an unknown and change no answer; what
    was missing was the energy transfer, and `AT-1d` (ii) now guards it.
    (Section 1.6, T1.9.)
14. **What a refused corrected product suppresses. Decided: suppression is by
    quantity, not by file.** The question was whether a `_adv` file should
    suppress only the affected columns or refuse a corrected profile
    altogether. Every column states whether it is corrected or uncorrected, and
    a column that depends on a refused one is refused with it. (Section 8,
    T8.1.)

---

## 11. What this specification does not establish

- Nothing here was measured. Every code number is READ with its line, at the
  reading time the provenance note states; every published number is quoted
  from the publisher PDF.
- The independent species space (T1.7), the validity states of the enumerated
  step (T7.4) and the boundary equations by cell class (section 9.2) are READ
  from B1a's inventory and carry its provenance, not an independent reading of
  the tree.
- Section 5.2 names the form of a CO destruction model but supplies no rate
  coefficients: no compilation for CO destruction in this regime was
  consulted, and none is quoted from memory. Decision 7 adopts exclusion, so
  nothing in the specification depends on those rates; obtaining them is what
  reopens alternative A.
- Section 1.6 states what the oxygen ledger must contain and which three gaps
  it closes; it does not supply the species enthalpies themselves. `eps` for
  OH, H2O, CO, O I and O(1D) must be taken from a published table with the
  reference state of T1.2 stated, and none is quoted here.
- Section 8 states the rule for the derived products; it does not enumerate
  which channels of the present post-process would have to be added to satisfy
  it beyond those `C8` already names.

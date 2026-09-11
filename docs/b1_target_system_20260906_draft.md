# B1 (draft): the target governing system of EXHALE

**Superseded 2026-09-06 by `docs/b1_target_system_20260906.md`, the accepted
specification; this file is the draft it was accepted from and is kept for the
record.**

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

**Line numbers drift.** The tree is under edit by parallel steps. Every line
cited here was read on 2026-09-06, but a number is a pointer into a moving
file and the routine name is the durable citation: between the first draft of
section 7 and this filling, `ioniz_eq` in the marching loop moved from
`EXHALE_main.f90:1223` to `:1244`, the composition projection from
`:1226-1238` to `:1247-1259`, the Shapiro call from `:1283-1287` to `:1306`
and the adoption boundary from `:1291` to `:1312`; the CO ceiling counters
moved from `diffusive_photochemistry.f90:613-618` to `:630-633`. The section 7
table keeps the numbers of its own reading and this note states the offset.

**Provenance.** Every code statement is READ from the tree at
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` on 2026-09-06, with
file and line. Nothing here was MEASURED: no build, no run, no test. Every
published equation is quoted by its number from the publisher PDF under
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/references/`.

**Owners.** Every target carries `Owner: (unassigned)`. The advisor assigns
them when the user accepts this specification.

---

## 1. Energy convention and the formation/excitation ledger

Replaces **D0 C1** (the temperature-preserving composition projection,
`EXHALE_main.f90:1226-1241`) and **D0 C3** (the two-iteration semi-implicit
energy update with a lagged cooling, `energy_semi_implicit.f90:190, 210-211,
272-278, 294-295`). Answers D0 open questions 1 and 2, in the direction
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
Riemann solver already use (`caloric_eos.f90:30-33`), and keeping it means
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
| metal ion of stage `k` | sum of the `k` ionization potentials | cumulative |
| He 2^3S | `E(2^3S)` above He I ground | metastable excitation |
| H(n=2) | `E_21` above H(1s) | excitation, both `2s` and `2p` |
| H2 | `-D0(H2)` | bound below two free H atoms |
| H2+, H3+, HeH+ | formation energy from the reference atoms | binding plus ionization |
| OH, H2O, CO | formation energy from the reference atoms | binding |

**T1.3 (the non-overlap rule).** Energy already held in the equation of
state is not also in the reservoir, and the reverse. The rule is enforced by
fixing one zero per internal ladder:

- `u_rv(T)` in `T1.1` is the rovibrational energy of **bound** H2 above its
  own `v = 0, J = 0` level. `caloric_eos.f90:70-73` states this explicitly
  ("Dissociation is not a heat capacity and is not in here: it is a chemical
  source term").
- `eps(H2) = -D0(H2)` is measured from `v = 0, J = 0` as well:
  `mol_rates.f90:102` carries `D0_H2_cm = 36118.11 cm^-1` labelled
  `D0(v=0,J=0)`.

The two therefore share one zero and cannot double count. The same
statement is required, and is the acceptance condition, for every other
species with a level ladder: exactly one of the EOS and the reservoir owns
each level's energy, and the specification names which.

**T1.4 (H(n=2) and He 2^3S).** Both are excited states with an explicit
population in the code (`n2s_arr`, `n2p_arr`, `parameters.f90:1239-1240`;
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

- The pressure reset `comp_p_from_T` at `EXHALE_main.f90:1230` followed by
  `W_to_U` at `:1238` disappears as an independent update. Pressure after
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
configuration selects (`input_read.f90:1822-1852`). Everything outside that
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
| the neutral stage of each element | the element total, neutral stage by difference | imposed by the layout of the fractions | `System_HeH_mol.f90:522, 528`; `ion_residual_core.f90:167` |
| `n_e` | charge neutrality, `n_e = sum_s z_s n_s` | imposed as an identity; never an unknown and never a row | `utilities.f90:201-254`; `System_HeH_mol.f90:539` |
| `rho` | none in the target: the density is supplied to the step and preserved (T2.1) | checked against the species mass sum (T2.2) | today overwritten from `calc_rho`, `ionization_equilibrium.f90:2152-2157` |
| O(1D) | local steady state on a fast process: one sink (O6) and a 1.6e-3 s chemical lifetime, so the flux through O6 equals the production `oj4 n_H2O` | imposed; the only fast-process elimination the code justifies at its site | `System_HeH_mol.f90:408-418` |
| CO inside the cell solve | the transported value, or the thermal equilibrium ceiling; no kinetic destruction row | imposed after the solve | `ionization_equilibrium.f90:716-726`; `limit_to_element_budget`, `diffusive_photochemistry.f90` |
| H(n=2), both `2s` and `2p` | a 2x2 statistical equilibrium solved outside every solver, one outer pass before the state it closes | imposed, lagged | `excited_hydrogen.f90:234-289`, called `EXHALE_main.f90:1219` |
| a transported carrier inside the cell solve | the identity row `x_i - x_i^fix` carrying the transport operator's value | imposed | `System_HeH_mol.f90:565, 574`; `ionization_equilibrium.f90:1244-1246` |

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
- a **transported** species keeps its identity row in the local system, so
  that the Jacobian keeps its dimension and every other row sees the imposed
  value; its balance is a grid-coupled equation owned by B4 and is not a
  local unknown. Whether that row is instead removed from the local space
  entirely is open question 11.
- the **constrained continuation** inverts the roles: the unknowns become
  `u_k = ln n_k`, the element totals become explicit residual rows and charge
  neutrality holds identically at every iterate
  (`constrained_chemical_equilibrium.f90:1124-1126, 1198-1244`). The target
  requires the energy row (b) of T1.6 to be carried in that branch as well, on
  that same reduced space measured in `ln n_k`. A continuation branch that
  solves the composition alone re-introduces the `C1` projection through the
  back door, whatever section 1 does to the main path.
- **He 2^3S is an unknown** (`tr_triplet_row`, `ion_residual_core.f90:81-95`)
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
`System_HeH_mol.f90:399-405` with the identifiers of `docs/a2_reaction_audit.md`
section 2: `O1  OH + H2 -> H2O + H` and its reverse `O1r`, `O2  O + H2 -> OH + H`
and its reverse `O2r`, the three water photolysis branches `O3  H2O + hv -> OH + H`,
`O4  H2O + hv -> H2 + O(1D)` and `O5  H2O + hv -> O + H + H`, the hydroxyl
photolysis `O7  OH + hv -> O + H`, and `O6  O(1D) + H2 -> OH + H`, which is the
single sink that lets O(1D) be eliminated (T1.7). Of these, only the photolysis
channels deposit energy today: `heat_fuv` sums the excess of each absorbed band
photon over the bond energy for H2O and OH
(`ionization_equilibrium.f90:2123-2134`, with the O(1D) statement in the comment at `:2113-2122`). The four collisional channels O1, O1r,
O2, O2r deposit nothing, and the code states at the same site that the 1.96 eV
of electronic excitation carried out of O4 by O(1D) is deliberately kept out of
the photolysis heat and that its release in O6 "is NOT deposited by this
network". `molecular_reaction_heat.f90:69-73` records the third gap: the
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
`molecular_reaction_heat.f90:43-48` already states for the H2/He network
("Not as a list of reaction enthalpies, a list can disagree with itself, but
from ONE table of species enthalpies"). The requirement is therefore not a new
term but the extension of the one table to the oxygen and carbon species, and
the deletion of any channel-specific oxygen heat that would then be a second
count.

Three things the extension must state, because each is a place where the
present code loses energy:

1. **An elimination transfers the reservoir with the nuclei.** O(1D) is
   eliminated by the fast closure of T1.7, so the flux through O6 carries both
   the oxygen nucleus and `eps(O(1D))` to the products at the same rate. A
   closure that moves the nuclei and drops the excitation energy is exactly the
   `C25` gap, and it is the general rule for every fast-process elimination,
   not a statement about oxygen.
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

Replaces **D0 C2** (`ionization_equilibrium.f90:406` declares the density
`intent(inout)` and `:2152-2157` overwrites it with `calc_rho` of the
equilibrium composition; `EXHALE_main.f90:1236` pushes that into the
conserved mass row; the code's own `EXHALE_UPDATE_MAP` reports it as a
`chem` column absent from the steady residual,
`EXHALE_main.f90:2620-2634`). Answers D0 open question 3 in the direction
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

with `m_s` from the one species mass table and `tol_mass` derived from the
declared electron-mass convention of that table (whether the electron mass
is carried with the ion or as a separate term). The convention is stated at
exactly one site and `tol_mass` follows from it; it is not tuned until an
existing state passes. A violation is a failure of the composition solve
(reported through B6), never a correction applied to `rho`.

**T2.3.** Reconstruction of `rho` from the composition is a diagnostic. The
steady equations gain **no** chemical mass source: the fix is to remove the
term from the production map, not to add a matching fictitious term to the
residual, so that the marching map and the stationary map have the same
fixed point.

**Interface consequence.** The `n_io` argument of `ioniz_eq` becomes
`intent(in)`. That is a public-interface change inside the source tree and
is listed in section 10 as a point where the user's confirmation is wanted
before coding.

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
escape probability (`Cool_coeff.f90:3596-3601`, the Cen-type fit
`7.5e-19/(1+sqrt(T/1e5)) exp(-118348/T)`, summed into `coex` at
`util_ion_eq.f90:1719-1721`), in the same cells where `lya_rt` computes an
escape probability far below one, while the collisional de-excitation heat
returns `E_21` per de-excitation (`excited_hydrogen.f90:212-213`,
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

READ, `lya_rt.f90:251-340, 645-691`. The code builds

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
`Jbar = S (1 - beta_tot)` (`lya_rt.f90:715-716`). The two FACES of the slab are
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
(`lya_rt.f90:228-241`, `excited_hydrogen.f90`):

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
the same denominator (`lya_rt.f90:238-241`), so the target reuses one
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
background with no counter-flux, `binary_element_diffusion.f90:1958-2071`),
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
`binary_element_diffusion.f90:1456-1466` builds
`eEf = -kb_erg * TK * dln(n_e T)/dr` by central difference, one-sided at the
ends, gated by `he_ambipolar`; the routine's own header
(`:1417-1418`) states the same expression. The target therefore **adopts the
field the code already has** and supplies the derivation it was missing,
together with the conditions under which it holds.

**Omitted order of (F).** Three terms are dropped and each is named:

1. electron inertia, `O(m_e/m_i)` relative to the retained pressure term;
2. the thermal force, `alpha_T dlnT/dr`, which the code carries as a separate
   optional term with `he_alphaT` defaulting to zero
   (`binary_element_diffusion.f90:1493-1499`), and which Koskinen et al.
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
  can omit it: `diffusive_photochemistry.f90:80-85` states that "All four
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
`solve_trace_element_in_hydrogen` (`binary_element_diffusion.f90:1958-2071`)
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
(B9) is satisfied for those components too. Either the trace elements are
solved jointly with the H/He binary in one matrix inversion, which satisfies
(B10) by construction, or a sequential solve is retained and its residual
`| sum_s rho_s w_s | / ( sum_s rho_s |w_s| )` is measured and bounded at a
stated order. A trace approximation is retained only with that measured
residual and a stated validity domain.

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
`Pe_s` to be evaluated per cell for every species held at `w_s = 0` and
reported through B6 where it falls below a stated threshold.

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
routine's own header names (`binary_element_diffusion.f90:1417-1421`):
`m_He/m_H - 1` in neutral gas, `m_He/m_H - 3/2` in an H+ plasma,
`2 m_He/3 m_H - 1` in a He++ plasma, that is 2.9715, 2.4715 and 1.6477 with
`m_He/m_H = 3.9715`.

**AT-4d (element conservation in a transient).** A diffusion-only column run
through a transient conserves each element's total nucleus count to
round-off, which the conservative form of T4.5 gives and the advective form
of `C28` gives only to hydro truncation.

---

## 5. CO: constraint model or exclusion

Replaces **D0 C26** (CO has no balance row and no kinetic destruction rate:
an inert transported reservoir until an algebraic equilibrium ceiling removes
it). Follows rev 3 section 4.2 item 9 and D0 open question 5.

### 5.1 What the code does now

READ. With the oxygen chemistry on, the CO density is either the local
equilibrium value `co_equilibrium_density(n_C, n_O, T)` or, when carrier
transport is on, the transported value, and it is then clipped to the carbon
and oxygen the cell has:

```text
nCO = min( nox_eq(:,3), n_el(C), n_el(O) )         (ionization_equilibrium.f90:723-725)
n_tot(O) = max( n_el(O) - nCO, 0 ) ,  n_tot(C) = max( n_el(C) - nCO, 0 )   (:732-734)
```

with a second ceiling in `limit_to_element_budget`
(`diffusive_photochemistry.f90`). The cumulative record already exists:
`co_ceiling_CO_removed` (`diffusive_photochemistry.f90:617`, molecules
removed) and `co_ceiling_seen` (`:618`, per cell).

### 5.2 Alternative A: a one-sided destruction model

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
ranges and `J_CO` with its band data. This draft does **not** supply them:
no rate compilation for CO destruction in this temperature and density
regime was consulted here, and a specification that names rates without
reading their source would violate the citation rule. This is a blocked
part, listed as such in the report.

### 5.3 Alternative B: exclusion, with the record as the observable

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

**Recommended position**, following section 5.3 of
`To_be_determined_by_user_recommend_20260906.md`: adopt B now, hold A open.

**Owner:** (unassigned).

**Acceptance test AT-5.** (i) A run whose ceiling never fires reports an
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
piecewise fits join discontinuously; `h3p_cooling.f90:134-135, 46-61`).
Answers D0 open question 6. This section states **what the target requires of
the published fit definitions**; it does not re-derive the fits.

READ, `h3p_cooling.f90:1-40`: the code uses Miller, Stallard, Tennyson and
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
column at 6 (`h3p_cooling.f90:135`), because it leaves a finite emission
where there are no colliders. In the low-collider limit, where radiative
decay empties the emitting levels faster than collisions populate them, the
emission per molecule is set by the collisional excitation rate and is
therefore linear in `n(H2)`. Below the table's lowest column the target
either extrapolates on that limit, `s ~ s(T, 10^6) * n(H2)/10^6`, or refuses
to evaluate and reports an out-of-domain closure through B6. Which of the two
is chosen is stated with the choice.

**T6.3 (temperature limits and extrapolation policy).** The four published
segments do not join, and the mismatch is the paper's: the code's own header
records that continued across 300 K the 300-800 K polynomial stands 15 to 43
per cent above the 30-300 K one over 200-300 K, and the 1800-5000 K
polynomial 2.4 per cent below the 800-1800 K one at 1800 K
(`h3p_cooling.f90:14-20`, a value the code site records; not re-measured
here). The target states the policy rather than smoothing it away: each
segment is used only on its published range; below 30 K and above 5000 K the
evaluation is clamped and the cell is marked as evaluated outside the fit's
domain through B6; the discontinuity at a segment boundary is a known
property of the published fits and is reported, not hidden by an
interpolation that no publication supports.

**Acceptance test AT-6a.** (i) Low-collider limit: `Lambda / n(H3+) -> 0`
proportional to `n(H2)` as `n(H2) -> 0`. RED today by the clamp. (ii) LTE
join: for `n(H2) >= 1e14 cm^-3`, `Lambda = n(H3+) 4 pi E(T)` within the
tolerance the specification states. (iii) Domain: a cell at 20 K or at
6000 K produces a B6 out-of-domain record.

### 6.2 H2 thermodynamics

Replaces **D0 C22** (chemistry and the caloric equation of state use
different H2 internal-state models, so one thermodynamic potential does not
generate both).

READ. The EOS evaluates `u_rv` and `c_rv` from the complete observed bound
rovibrational ladder of H2 X^1 Sigma_g^+, 302 levels to 51966 K, from Roueff
et al. (2019, A&A 630, A58) table 2, with the nuclear-spin weights inside the
sum, hence the ortho/para equilibrium mixture
(`caloric_eos.f90:49-68`). The chemistry evaluates its partition sum
`q_rovib_H2` from Huber and Herzberg spectroscopic constants
(`mol_rates.f90:102-106`: `D0 = 36118.11 cm^-1`, `omega_e = 4401.213`,
`omega_e x_e = 121.336`, `B_e = 60.853`, `alpha_e = 3.062`), tabulated on
100-20000 K and interpolated (`mol_rates.f90:108-130`).

**T6.4 (one potential).** One partition function `Q_H2(T)`, from one level
set, with one zero of energy at `v = 0, J = 0`, generates both:

```text
u_rv(T) = k_B T^2 d ln Q_H2 / dT              (the EOS)
ln K_eq(T) = ... ln Q_H2(T) ... - D0 h c / (k_B T)   (the chemistry)
```

so that the forward and reverse rates, the caloric equation of state and the
reaction enthalpy of section 1's reservoir are thermodynamically consistent
by construction rather than by coincidence. The zero of the ladder is the
same zero the reservoir of `T1.3` uses for `eps(H2) = -D0`, which is what
makes section 1's non-overlap rule verifiable for H2.

**Acceptance test AT-6b.** (i) Van 't Hoff: `-d ln K_eq / d(1/T)` computed
from the shared `Q_H2` reproduces the reaction enthalpy assembled from the
section 1 reservoir plus the EOS heat capacities, to a stated tolerance,
over 100-5000 K. (ii) `c_v(H2)` from the shared `Q_H2` is reported against
both of today's models, and the difference is stated rather than assumed
zero: the code site already records that the observed ladder exceeds a rigid
rotor plus harmonic oscillator by 1.4 per cent at 200-900 K, 2.3 per cent at
1300 K and 7.3 per cent at 3000 K (`caloric_eos.f90:56-61`, code site
records). (iii) A run with no molecules reproduces the constant `gamma_ad`
arithmetic exactly, the structural branch the EOS already carries
(`caloric_eos.f90:76-84`).

**Owner:** (unassigned).

---

## 7. Rollback contract

Rev 2 section 5.1: the B1 checkpoint must cover **every** state-changing
operation before the final adoption boundary, viscosity, conduction and the
optional Shapiro filter included, or state a narrower interim contract.

### 7.1 The enumerated step

READ from `src/EXHALE_main.f90`, in the order the main program executes
them. The **adoption boundary** is line 1291 (`update_map_end_step`): from
line 1294 onward the loop reads the state it has adopted and computes
diagnostics and control flags. Everything above 1291 is inside the attempted
step and must be restorable.

| # | operation | lines | state it changes |
|---|---|---|---|
| 0 | `eval_dt(W,dt,dt_loc)`, update-map scaling of `dt`, `dt_loc` | 992-1001 | `dt`, `dt_loc` (main) |
| 1 | `u_old = u` | 1015 | the checkpoint anchor; the retry loop restarts every attempt from it (`:1008-1011`) |
| 2 | `retry_step` loop: three RK stages, `Apply_BC`, `positivity_limited_fluxes` | 1057-1170 | `u`, `u1`, `u2`, `WL`, `WR`, `dF`, `S` (main); `n_faces_flux_positivity_limited` (`RK_rhs.f90:29`), `n_faces_flux_positivity_limited_accepted` (`RK_rhs.f90:45`), `n_dt_halve`, `n_steps_dt_halved` (`EXHALE_main.f90:288`), and `dt`, `dt_loc` on a halving |
| 3 | `U_to_W`, `get_species_densities`, `comp_T_from_p` | 1187-1197 | `W`, `rho`, `v`, `p`, `nhi`...`nheiTR`, `nm`, `ne`, `n_tot`, `T` (main) |
| 4 | `element_diffusion_step` (`he_diffusion`) | 1204 | `f_sp` (main); `he_fraction_over_one`, `he_fraction_under_zero` (`binary_element_diffusion.f90:313-314`), `he_fraction_newton_steps`, `he_fraction_newton_resid` (`:317-318`), `trace_ratio_under_zero` (`:323`) |
| 5 | `photochemical_transport_step` | 1214 | `f_sp` carrier rows (main); `diffusive_photochemistry.f90`: `n_carrier` (`:455`), `carrier_solved` (`:465`), `cbg_*` saved backgrounds (`:520-528`), `cph_klw`, `cph_jh2o`, `cph_joh` photolysis arrays (`:532-533`), `adv_corr` (`:569`), `fc_base` (`:570`), `pct_newton_resid` (`:574`), `pct_newton_resid_before_limit` (`:581`), `pct_cell_constrained` (`:603`), `pct_verdict` (`:607`), `co_ceiling_CO_removed` (`:617`), `co_ceiling_seen` (`:618`), `pct_worst_limit` (`:624`), `pct_Dco` (`:628`), `row_terms` (`:633`), `col_scale_H2`, `row_scale_H2` (`:652`), `headroom_H2` (`:655`) |
| 6 | `excited_H_update` (`use_excited_H`) | 1219 | globals `gph_balmer_HI`, `heat_balmer` (`parameters.f90:1235-1236`), `Jlya_arr` (`:1238`), `n2s_arr`, `n2p_arr` (`:1239-1240`), `Sproton_arr` (`:1241`), `Hpe_arr`, `Hdx_arr` (`:1242-1243`); `excited_hydrogen.f90`: `jlya_rt_loaded` (`:58`), `jlya_rt_grid` (`:59`), `Tdiag`, `nhidiag`, `nediag` (`:62`), `nhiidiag` (`:63`), `taulya` (`:64`); `lya_rt.f90`: `jint_arr`, `jstar_arr` (`:60`) |
| 7 | `ioniz_eq(T,rho,f_sp,heat,cool,eta)` | 1223 | `f_sp`, `heat`, `cool`, `eta` **and `rho`** (the `C2` write, `ionization_equilibrium.f90:406` `intent(inout)`, overwritten from `calc_rho` at `:2152-2157`); module state: `nmol_eq` (`:53`), `NH2_col_lw` (`:61`), `f_shield_lw` (`:62`), `k_lw_diss` (`:63`), `p_lw_single`, `p_lw_absorbed` (`:69`), `tr_lines_lw` (`:78`), `P_H2_eq` (`:82`), `nox_eq` (`:110`), `n_o1d_eq` (`:111`), `NH2O_col`, `NOH_col` (`:112`), `j_h2o_fuv`, `j_oh_fuv`, `tau_fuv` (`:113`), `heat_fuv` (`:114`), `heat_chem` (`:117`), `bg_cell` (`:127`), `bg_cell_adopted` (`:149`), `bg_cell_best` (`:150`), `ieq_sweep_state_kind` (`:176`), `ieq_marching_ledger` (`:238`), `ieq_steady_iterate_ledger` (`:239`), `ieq_steady_candidate_ledger` (`:240`), `ieq_nonroot_streak` (`:288`), `ieq_acc_nprint` (`:309`) |
| 8 | `get_species_densities`, `comp_p_from_T`, `W` assembly, `W_to_U` | 1226-1238 | `ne`, `n_tot`, `p`, `W`, `u` (main). **This is the `C1` projection** and section 1 removes it as an independent update |
| 9 | `solve_energy_semi_implicit` (or the explicit `u(3,:)` update at `:1256`) | 1253 | `u` (main); `energy_semi_implicit.f90`: `n_energy_floor_hits` (`:30`), `energy_floor_first_step` (`:31`), `energy_floor_last_step` (`:32`), `energy_floor_cell_hits` (`:33`) |
| 10 | `Apply_BC(u)` | 1259 | ghost cells of `u` |
| 11 | `viscous_conduction_step` with its `U_to_W`, `comp_T_from_p`, `Apply_BC` (`transport_active()`) | 1272-1276 | `W`, `rho`, `v`, `p`, `T`, `u` (main); `viscous_conduction.f90`: `n_conduction_floor_hits` (`:165`), `conduction_floor_first_step` (`:166`), `conduction_floor_last_step` (`:167`), `conduction_floor_cell_hits` (`:168`) |
| 12 | `shapiro_filter(u)`, `Apply_BC(u)` (`shapiro_eps > 0` and on the `shapiro_every` cycle) | 1283-1287 | `u` only; the filter holds no module state (`Apply_BC.f90:271-288`) |
| -- | **adoption boundary** | 1291 | |
| 13 | post-adoption reads and control state: `U_to_W`, `check_base_inflow_is_subsonic`, `get_species_densities`, `comp_T_from_p`, `du` and arming flags, residual diagnostics, staged secondary-ionization flips at `:1649` and `:1698` | 1294-1441, 1649, 1698 | control and diagnostic state, not the physical state; `sec_ion_active` (`parameters.f90:179`) and `sec_ion_armed_step` are mode state and belong to A0, not to the rollback |

The update-map snapshots `u_umA` (`:1179`), `u_umB` (`:1241`), `u_umC`
(`:1261`), `u_umD` (`:1279`) mark the same operator-split seams and are the
existing instrumentation of this boundary. `EXHALE_main.f90:1171-1178`
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
`co_ceiling_CO_removed`, `bg_cell_adopted`, `ieq_nonroot_streak` and the
floor-hit counters of rows 9 and 11.

**T7.2.** Counters are partitioned into three classes, and the checkpoint
treats them differently:

- **physical accumulations** (the CO record, reaction and heating budgets,
  elapsed physical time): restored on rejection, so a rejected trial leaves
  no contribution, which is the B6 semantics of review 2 section 5.3;
- **attempt statistics** (attempted steps, retries, halvings): deliberately
  **not** restored, since they count attempts;
- **diagnostic extrema** (`he_fraction_over_one`, `trace_ratio_under_zero`,
  floor-hit counters): the specification states one policy per counter and
  the reason, since an extremum reached only in a rejected trial is
  informational and must not invalidate the accepted result.

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
| 2 | RK stages, `Apply_BC`, positivity-limited fluxes | (4) reconstruction faces scaled toward the cell average and ghost extrapolation dropped to zero gradient, in an accepted step; (3) the same repairs inside an attempt the `dt` bisection discards | no for the state records; the attempt counts stay. The code already separates the two (`calls`, `faces`, `faces in accepted steps`, `EXHALE_main.f90:2306-2321`) and that separation is the model the rest of the table follows |
| 3 | `U_to_W`, `get_species_densities`, `comp_T_from_p` | none of its own; an inadmissible state here is a rejection, not a validity state | not applicable |
| 4 | `element_diffusion_step` | (4) the helium fraction range clip and `trace_ratio_under_zero`; (2) Blanc's law and the Chapman-Enskog first approximation, both stated at the code site and neither checked; (5) the inner Dirichlet reservoir of the operator | no. The clip counters are diagnostic extrema and follow the T7.2 policy, which is open question 10; the (5) record is a property of the configuration and is not produced by the step at all |
| 5 | `photochemical_transport_step` | (1) the CO ceiling, if it is read as unvalidated physics; (4) the CO ceiling as a correction with no rate behind it, and the carrier element rescale of `limit_to_element_budget`; (3) a carrier trial above the headroom, which is refused rather than clamped; (2) the H2 self-shielding table clamps | no. The CO record is a physical accumulation under T7.2 and is restored, which is exactly what `AT-5` (iii) tests; the refusal count is an attempt statistic and stays. Which of (1) and (4) the ceiling is remains open question 12 |
| 6 | `excited_H_update` | (2) the damping-wing slab closure evaluated where `(a tau)^(1/3) > 10` fails, and the Sobolev channel outside its own domain; (1) the H(n=2) lag, which is a closure evaluated one outer pass behind the state it closes and whose lag error is not measured | no; both are properties of the attempted state |
| 7 | `ioniz_eq` | (1) the oxygen ledger gap (section 1.6), the H3+ model domain, the H2 caloric and chemistry state mismatch, the metal Voronov transcription note; (2) every clamped closure of B1a section 4.2 that the sweep evaluates, the CHIANTI coronal cutoff, the Ca II and C I range clamps, the infrared band clamps, the Penning and Oklopcic-Hirata ranges, the secondary-ionization extrapolation, the Gaunt table clamp; (3) steady candidate probes, which are already tagged `ieq_state_steady_candidate`; (4) the simplex clamp accepted as class 3, the class-4 relaxation amnesty, and the `rho` rewrite until section 2 lands | the (1) states are static properties of the configuration and are not produced by a step, so a restore does not touch them; every (2) and (4) record produced inside a rejected attempt is discarded with the attempt; the probe tags are attempt statistics and stay |
| 8 | the composition projection, `comp_p_from_T` and `W_to_U` | (4) the `C1` projection itself, in every configuration, with no run record today | no. Section 1 removes the producer; until then its activity refuses certification of a physical step (A2 section 5) |
| 9 | `solve_energy_semi_implicit` or the explicit update | (4) the energy floor clamp; until B2, the absence of a residual test on the returned state is itself a (1) state of the marching energy update | no |
| 10 | `Apply_BC(u)` | (5) the base `(p, s)` reservoir, the base composition handoff and the base infrared field; (2) a supersonic base inflow, which over-specifies the characteristic boundary and is diagnosed only today; (4) the outer ghost extrapolation dropped to zero gradient | the (5) records are configuration properties and always survive; the (2) and (4) records of a discarded attempt do not |
| 11 | `viscous_conduction_step` | (4) the conduction floor clamp | no |
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
no CO ceiling record and no reaction budget, contains a contribution from a
rejected attempt. (iv) The attempt statistics do contain it.

---

## 8. Derived products: the `_adv` rule

Replaces **D0 C8** (the post-process temperature is solved with a different
channel set from the hydrodynamic temperature: no H3+, no molecular infrared
bands, molecular ions absent from `n_e`, and on the heating side photoheating,
the He recombination coupling and the He(2^3S)+H Penning heat only, so no
Balmer, no Lyman-Werner, no molecular chemical and no FUV heat;
`T_equation.f90:81-90, 111-206`, `post_process_adv.f90:820-851`) and **D0 C20**
(the `_adv` profiles the transit tools read are reconstructed on a
molecule-free background, and their molecular columns are the equilibrium
sweep's values; `post_process_adv.f90:5-19`, the limitation printed into the
file by `write_output.f90:139-150`). D0 open questions 11 and 12 are the two
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
channel was missing, and the uncorrected profile keeps its own name. Writing a
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
the whole missing channel. `AT-8h`: an `_adv` product made from an `init`
snapshot carries `mode=init` and no `t_phys`; made from a physical state it
carries `mode=phys` with the state's `t_phys`; a product whose source state
carries a static exclusion (oxygen ledger, CO ceiling) carries the same
exclusion, and a tool reading the file can tell.

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
| AT-5 | section 5 | ceiling inactive, active, active-in-a-rejected-trial, active-in-initialization | exclusion only for accepted physical activity; record auditable | rev 3 item 9; review 2 section 5.3 |
| AT-6a | section 6.1 | `n(H2) -> 0`; `n(H2) >= 1e14`; out-of-range `T` | linear vanishing; LTE join; B6 out-of-domain record | Miller et al. (2013) Tables 5, 6 |
| AT-6b | section 6.2 | van 't Hoff on the shared `Q_H2`; `c_v(H2)`; molecule-free run | consistency to the stated tolerance; difference reported; constant `gamma_ad` exact | Roueff et al. (2019) table 2 |
| AT-7 | section 7 | injected failure after two different operations | complete restoration; correct elapsed time; no rejected contribution in a physical accumulation; attempt statistics keep it | rev 2 section 5.1 |
| AT-8a | B3, physical mode | temporal refinement: successively halved steps to one fixed end time | observed temporal order matches the claimed one; splitting error controlled | rev 3; D0 section 6.6 records that no such measurement exists |
| AT-8b | B3c | splitting-order check: operator order reversed at fixed end time | difference at the order the splitting claims | |
| AT-8c | B2 | returned-temperature residual on a falling cooling branch, a strongly temperature-dependent heating, an unavailable bracket, a floor encounter | final residual tested; a floor is a specified reservoir with budgeted heat or a defined failure | |
| AT-8d | B5 | species stationary residual for every active transported balance | evaluated on the independent space of T1.7 for that configuration, same source assembly as B3c, same face fluxes as B4 | B1a section 2.2 |
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
| 5 | oxygen | the rows of 4, plus AT-1d (section 1.6), AT-5, and AT-8d for the four transported balances | AT-4a to AT-4d | the states of 4, plus the static oxygen ledger exclusion (T1.10) and the CO ceiling record, which is the `AT-5` observable |
| 6 | diffusion | the rows of 3, plus AT-4b, AT-4c, AT-4d, and AT-8d for the elemental balance | AT-4a (no trace-metal rows: `He_metal_diffusion` is not set), AT-5 | the states of 3, plus Blanc's law and the Chapman-Enskog first approximation, both stated and unchecked, and the advective form of the elemental transport (`C28`) |
| 7 | profile | AT-1a, AT-1b, AT-1c, AT-2, AT-3, AT-4a, AT-4b, AT-4c, AT-4d, AT-6b (iii), AT-7, AT-8a to AT-8h | AT-5, AT-6a, AT-6b (i) and (ii) | the states of 2, plus `C28` and `C29`, the trace-metal rows with no counter-flux, which this configuration owns |
| 8 | H+ transport | the rows of 4, plus AT-8d for the H2 and H+ balances | AT-4a to AT-4d, AT-5 unless oxygen is added | the states of 4, plus `C18`, the two base ghosts re-solving H and H+ locally, and `C19`, no steady solver for the configuration, so certification is only through a marching solution |

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
| base reservoir ghosts | `1-Ng .. 0` | the characteristic face condition at `r_edg(0)`: the LODI `C^-` compatibility relation plus the `(p, s)` reservoir, with the number of conditions set by the eigenvalue signs, so `0 < v < c` gives two reservoir conditions and one interior relation (`base_boundary.f90:25-33, 281-391`); the ghosts are volume averages of the hydrostatic isentrope through the face state, not copies of it (`Apply_BC.f90:157-166`) | that the base inflow is subsonic, since a supersonic inflow over-specifies the boundary (`check_base_inflow_is_subsonic`, `Apply_BC.f90:76-100`) | the subsonic condition becomes a stated domain of the boundary closure: outside it the state carries a category (2) validity record and certification is refused, in place of the present diagnosis only (`C33`) |
| base reservoir ghosts, composition | same | the H2 Dirichlet inflow composition, applied as an advective inflow term when the base face flux is inflowing and a handoff states the partition (`diffusive_photochemistry.f90:1243-1244, 1878-1886`); no diffusive flux at face 0; the diffusion operator's own inner Dirichlet reservoir (`binary_element_diffusion.f90:511-537`) | the element and charge constraints of the ghost composition | the proton deliberately gets no Dirichlet value and the imposed partition applies only for `j >= 1`, so the two ghosts re-solve their ionization locally (`C18`). That local equilibrium is itself a closure and the target requires it to be stated as one, with its residual evaluated by A2 like any other closure, instead of being an unstated difference between the ghosts and the first physical cell |
| physical cells | `1 .. N` | the hydrodynamic rows, the local rows of T1.7, and the transported balances of the configuration; in the stationary solve `set_base_fix` may additionally anchor the first `nfix` cells by replacing their rows with `Y - Yfix` (`steady_newton.f90:284-291, 363-373`) | the constraints of section 2 and of T1.7: element nucleus totals, charge neutrality, the species mass sum against the supplied density | an anchored row is an imposed condition and is declared as one in the certification record, so that a certified state says how many of its cells carried an equation of the physics and how many carried an anchor |
| outer ghosts | `N+1 .. N+Ng` | free outflow as a zero-gradient copy, with linear extrapolation added under WENO3 as an accuracy device (`Apply_BC.f90:190-217`) | that the extrapolation left `rho > 0, p > 0`; where it did not, the drop back to zero gradient is counted | the drop is a category (4) record of the state it was applied to, as it already is by count; the target adds that a state whose outer boundary is not outflowing is out of the closure's domain, since the free-outflow condition imposes no incoming characteristic |

---

## 10. Open questions for the user

Only where two readings are defensible.

1. **Reservoir zero for the metals.** Should `eps_s` for a metal ion be the
   sum of its ionization potentials from the neutral ground state, or from
   the element's dominant state at the base, which would make the reservoir
   smaller but tie the zero to a configuration?
2. **Electron mass convention (section 2).** Should the species mass table
   carry the electron mass with its ion, or as a separate species, given
   that the choice sets `tol_mass` and changes the meaning of `calc_rho`?
3. **`ioniz_eq` interface.** Making the density `intent(in)` changes a
   routine called from the marching loop, the steady residual and the
   post-process; confirm before it is coded.
4. **The residual H I collisional-excitation cooling (section 3).** After the
   `1s -> 2p` share moves into the 2p ledger, should the remaining Cen-type
   term be kept for the higher levels as it stands, or replaced by a fit that
   excludes Lyman-alpha explicitly?
5. **Trace metals in the diffusion matrix (section 4.4).** Solve the metals
   jointly with the H/He binary in one inversion, or keep the sequential
   solve with a measured and bounded (B10) residual?
6. **H+ diffusion (section 4.6).** Should H+ enter the (B9) set with its
   charge now, which changes every ionized-region profile, or stay at
   `w_s = 0` with the `Pe` diagnostic reported until B4 lands?
7. **CO (section 5).** Adopt exclusion now and hold the destruction model
   open, as recommended, or wait for the rate compilation before deciding?
8. **H3+ below the table (T6.2).** Extrapolate on the linear low-collider
   limit, or refuse to evaluate and mark the cell out of domain?
9. **The shared H2 partition function (T6.4).** Take the Roueff observed
   ladder as the single source, accepting that `K_eq` moves by the amount
   the two models differ, or keep both and pin them with a test, which
   section 5.3 calls useful in transition but not the permanent ownership
   model?
10. **Diagnostic extrema across a rollback (T7.2).** Should an extremum
    reached only inside a rejected trial be kept as informational, or
    discarded with the trial?
11. **The identity row of a transported species (T1.7).** Should a
    transported carrier keep its identity row `x_i - x_i^fix` inside the local
    Jacobian, so that the local system keeps one dimension for each species the
    cell contains, or should the row and its column be removed and the imposed
    value substituted into the remaining rows? The first keeps one layout for
    every configuration; the second gives a smaller and better conditioned
    local system but a dimension that changes with the transport flags.
12. **Which validity category the CO ceiling belongs to (T7.4).** It removes
    carbon and oxygen with no rate behind it, which reads as active
    unvalidated physics, and it moves chemical energy with no source, which
    reads as an unbudgeted accepted correction. The counters serve both
    readings; the category decides whether a ceiling-active run is labeled or
    excluded. B1a open question 3 is the same question from the inventory
    side.
13. **The reservoir of an eliminated excited species (T1.9).** Should O(1D)
    stay eliminated by its fast closure, with `eps(O(1D))` transferred to the
    products at the closure's own rate, or should it become an explicit species
    with its own `eps` and its own row, which costs an unknown and removes the
    need for the transfer rule? The first is the smaller change and the second
    is the one that cannot silently lose the excitation energy again.
14. **What a refused corrected product suppresses (T8.1).** When a channel
    active in the run cannot be evaluated in the post-process, should the
    `_adv` file suppress only the affected columns, keeping the rest of the
    corrected profile, or should it not write a corrected profile at all for
    that configuration? D0 open question 11 asks the same thing about the
    molecular columns that were never corrected.

---

## 11. What this draft does not establish

- Nothing here was measured. Every code number is READ with its line; every
  published number is quoted from the publisher PDF.
- The independent species space (T1.7), the validity states of the enumerated
  step (T7.4) and the boundary equations by cell class (section 9.2) are READ
  from B1a's inventory and carry its provenance, not an independent reading of
  the tree.
- Section 5.2 names the form of a CO destruction model but supplies no rate
  coefficients: no compilation for CO destruction in this regime was
  consulted, and none is quoted from memory.
- Section 1.6 states what the oxygen ledger must contain and which three gaps
  it closes; it does not supply the species enthalpies themselves. `eps` for
  OH, H2O, CO, O I and O(1D) must be taken from a published table with the
  reference state of T1.2 stated, and none is quoted here.
- Section 8 states the rule for the derived products; it does not enumerate
  which channels of the present post-process would have to be added to satisfy
  it beyond those `C8` already names.

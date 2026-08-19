# Oxygen chemistry in EXHALE: what exists, and the directions open to us

2026-08-19. Written to make the deferred item "(H) oxygen chemistry scope" a
choice between stated alternatives rather than an open-ended one. Nothing here
is implemented; the point is to separate what oxygen chemistry would buy from
what it would cost, so the scope can be fixed before code is written.

"Oxygen chemistry" covers two problems that share an element and almost nothing
else:

- **Problem A — the neutral base.** The H2-to-H partition at ~1 microbar is set
  by catalytic cycles on O, OH and H2O that EXHALE has no species for. This is
  the blocker recorded as item (H) in `TO_BE_DONE.md`.
- **Problem B — atomic oxygen in the ionized wind.** O I/II/III are already
  solved, cooled and charge-exchanged. What is missing here is observational
  (no O I line in the transit tool) and a few completeness questions.

They can be done in either order, or one without the other.

---

## 1. What the code does today

Verified against the source, not the documents.

| Piece | Where | State |
|---|---|---|
| Ionization stages O I / O II / O III | `species_table.f90` (`im_OI = 4 … im_OIII = 6`) | solved inside the coupled system with the other nine metals |
| Photoionization cross sections | `cross_sec.f90` `sigma_OI` (13.62 eV), `sigma_OII` (35.12 eV) | Verner et al. (1996) VFKY96 fits |
| Recombination | `Cool_coeff.f90` `bad_OI`, `bad_OII` | Badnell RR + adf48 DR |
| Collisional ionization | Voronov (1997) | as for every metal |
| Charge exchange | `charge_exchange.f90` cases 14, 33, 48 | A14 `O+ + H`, C6 `O+ + He`, D17 `Fe + O+` (Huang et al. 2023 Table 4) |
| Line cooling | `cool_OI_chianti`, `cool_OII_chianti` | closed-form CHIANTI v11 fits |
| Ground-term fine structure | slots `ifs_OI63`, `ifs_OI145`, `ifs_OI44` | O I $^3$P solved in exact statistical equilibrium, the three lines line-trapped |
| Diffusive settling | `species_diffusion.f90` | O settles independently when `He_metal_diffusion` is on (never exercised for O specifically) |

Not present:

- **No O-bearing molecule anywhere.** The molecular network is H2, H2+, H3+ and
  HeH+ only (`mol_rates.f90`, reactions R1–R23). There is no O, OH, H2O, CO or
  O2 in it.
- **No oxygen line in the transit tool.** `EXHALE_transit.py` carries He I
  10830, Ly-alpha, H-alpha, H-beta, Mg II h&k, Ca II H&K and Na I D.
- **No oxygen in the lower-atmosphere handoff.** `base.inp` carries exactly
  `T_base`, `r_base`, `HeH_base`, `Kzz_base`, `q_H2_base`, `p_base`.

One asymmetry worth keeping in mind: the O I/H I resonant charge exchange (case
14) is fast and nearly resonant, because IP(O I) = 13.6181 eV sits just above
IP(H I) = 13.5984 eV. The O I fraction therefore tracks the H I fraction closely
wherever that reaction dominates. This is what makes an O I line a probe of the
hydrogen structure, and it is also why an error in the O ionization balance
shows up as a hydrogen-like profile rather than an obviously wrong one.

---

## 2. Problem A — the base partition (item (H))

### 2.1 The gap, restated

EXHALE's molecular network settles at `q_H2 = 0.861` at a 1 microbar base where
VULCAN gives 0.75. H2 Lyman-Werner photodissociation was implemented to test
exactly this and moves `q_H2` by 0.5 % of the gap: the base partition is a
formation-destruction balance in which `n_H/n_H2` goes as the square root of the
destruction rate, so closing the gap needs a rate about 4100x larger, while the
H2 column of 4.3e21 cm^-2 suppresses the band by 1e6. The band physically
cannot reach 1 microbar.

### 2.2 Why oxygen is the answer, and the mechanism to implement

The destruction that photochemical codes find at that level is catalytic, with
water as the carrier. The cycle is

```
H2O + h.nu  ->  OH + H          (water photolysis, ~1150-1900 A)
OH  + H2    ->  H2O + H         (barrier ~2100 K)
-----------------------------------------------
net:  H2 + h.nu -> H + H,  with H2O and OH unchanged
```

The reason this works where direct photodissociation does not is that H2
self-shields its own Lyman-Werner bands, while water absorbs at longer
wavelengths that H2 does not block. Oxygen is thus not a new coolant or a new
opacity here -- it is a catalyst that lets ultraviolet photons destroy H2 far
below the LW front.

This is not a new reading. `docs/lower_atmosphere_coupling.md` already records
it as the thing a photochemical lower model uniquely provides (Lavvas &
Arfaux 2014/2021): OH-catalyzed H2 destruction dissociates H2 at 500-1000 K
where thermal equilibrium would not, and the H2-to-H transition can sit exactly
at the 1 microbar handoff level. What is new here is only that the wind code
would carry the cycle itself instead of importing its result.

**Before any of this is coded, the cycle and its branching must be checked
against the literature** (Lavvas & Arfaux for the hot-Jupiter context, Moses et
al. 2011 for the network, a kinetics compilation such as Baulch et al. 2005 or
the KIDA/UMIST databases for the OH + H2 rate and its barrier, and a photolysis
cross-section source such as the Leiden database or PHIDRATES for H2O). The
stoichiometry above is the standard one, but the branching ratios and the
temperature dependence of the barrier are what set the answer, and they are not
something to take from memory.

### 2.3 The options

**A0 — Status quo: pin `q_H2_base` from the handoff.**
Cost: none, it works today. `vulcan_to_base.py` interpolates VULCAN's H2/H/He
mixing ratios at the base pressure and writes `q_H2_base`; the molecular base
then uses that instead of the chemical-equilibrium fit.
Limit: exactly one number crosses the boundary, the composition is imposed
rather than computed, and there is no feedback from the wind onto it. This is a
deliberate modeling choice, not a patch, and it is defensible in a paper as long
as it is stated.

**A1 — Widen the handoff vector.** *(Already sketched: this is the `base.inp`
of the Tier-3 plan in `docs/lower_atmosphere_coupling.md`, whose example file
carries atomic-metal release fractions below the six keys that were actually
implemented. A1 is finishing that plan, not a new proposal.)*
`vulcan_to_base.py` already loads VULCAN's full species list and mixing-ratio
array; it simply does not export anything but H2/H/He. Export the O-bearing
species (O, OH, H2O, and CO if the carbon budget matters) and the metal element
abundances at the base level, as new `base.inp` keys. The metal side is the part
with the larger scientific return -- Lavvas finds Na, K, Mg, Ca and Fe returning
to atomic form below ~1e-3 bar, which is a physically motivated base abundance
for `metals.inp` in place of a scaled-solar total.
Buys: a base composition that is self-consistent with a photochemical model
rather than one number plus an equilibrium fit; a defensible O abundance at the
boundary; the possibility of a metals-from-below abundance rather than a scaled
solar one.
Costs: a few keys in the reader, the composition bookkeeping to decide what a
non-equilibrium O abundance means for `ntot_bc` and the electron budget, and a
decision about which species EXHALE is allowed to ignore.
Does **not** buy: any chemistry in the wind. The new species would be boundary
values, not transported quantities. In particular it does not close item (H) --
`q_H2` is still imposed.
Effort: small, mostly Python and input parsing. Default off, so old goldens
reproduce.

**A2 — A minimal catalytic set inside the coupled molecular system.**
Add H2O and OH (and neutral O as a *molecular-network* member, see the trap
below) to `System_HeH_mol`, with the handful of reactions that carry the cycle
above plus water photolysis with its own shielding.
Buys: the base partition becomes computed rather than imposed; item (H) closes;
the molecular base stops depending on an external code.
Costs: the coupled system grows by 2-3 unknowns in every cell where molecules
are active; new rate coefficients and cross sections, each of which has to be
sourced and documented like the existing R1-R23; a real risk of stiffness,
because water chemistry spans rates the H2 network does not; and the validation
burden of showing the extended system reduces exactly to the present one when
the oxygen abundance is set to zero.
**The design trap:** O I already exists as a trace-metal ionization stage. If a
neutral-O molecular species is added, the code would carry oxygen in two
solvers, and the two would drift apart -- exactly the failure mode the "one
kind of thing, one rule" convention exists to prevent. The right shape is
probably to let the metal system own all oxygen and have the molecular network
consume and return it through the shared `n_e`/particle budget, not to duplicate
it. This has to be decided before the first line is written.

**A3 — A reduced C/H/O photochemical network in-code.**
Fifteen to twenty-five species, i.e. a small VULCAN inside EXHALE.
Buys: full independence, and the ability to say the lower atmosphere and the
wind are solved with the same code.
Costs: large, and mostly in validation rather than in coding -- a network of
that size is only trustworthy after it reproduces a published benchmark. It also
changes the character of the code from a wind solver with a base condition into
an atmosphere solver.
Worth it only if a scientific claim depends on the base composition being
computed, not assumed.

**A4 — Two-way offline iteration, keeping the one-way handoff.**
Run VULCAN, hand off, converge the wind, feed the resulting upper boundary back
to VULCAN, repeat until the base state stops moving.
Buys: self-consistency between the two codes without new in-code chemistry.
Costs: a driver, a convergence criterion, and the wall-clock of several VULCAN
runs per planet; the result is only as good as the assumption that the two codes
overlap in a region where both are valid.
Does not close item (H) either -- it makes the imposed number self-consistent
instead of computed.

### 2.4 Which photochemistry code, if the answer comes from outside

A0, A1 and A4 all import the base composition from a photochemistry code. That
choice is not neutral, and it has already been measured:
`docs/vulcan_photochem_comparison.md` (2026-08-09, HD 189733 b) ran VULCAN and
Photochem with the same T(p), the same stellar spectrum, the same gravity and
the same elemental abundances.

| run | q_H2 | q_H at 1 microbar |
|---|---|---|
| VULCAN, `NCHO_photo_network` | 0.618 | **0.238** |
| Photochem, Zahnle set | 0.806 | **0.033** |
| Photochem, `NCHO_photo_network` (via `vulcan2yaml`) | 0.717 | **0.130** |

The spread in the neutral-hydrogen fraction is **7.2x**, and it splits into a
factor 4.0 from the *reaction network* and 1.8 from the *code* -- the latter an
upper bound, because Photochem's domain top had to be truncated. The network
dominates.

Two consequences for the options above.

- **The choice now reaches EXHALE.** When that memo was written the handoff
  carried four keys and the H2 dissociation state was not one of them, so the
  comparison ended "none of it reaches EXHALE". `q_H2_base` became a read key on
  2026-08-10 and `p_base` on 2026-08-15, so that conclusion no longer holds: the
  7.2x now propagates into the molecular-base particle count and therefore into
  the wind solve. The memo carries a dated note saying so.
- **For Problem A the network matters more than the code.** Photochem's own
  Zahnle set is built for rocky-planet atmospheres and gives the most molecular
  base of the three; the reading that it carries fewer H2 dissociation paths at
  ~860 K is inferred from the result, not from a reaction-level audit. Whichever
  code runs, it has to run a network that carries the H2O/OH catalytic paths of
  Section 2.2, or Problem A is answered by the network's omissions.

Reasons to move to Photochem anyway, none of them about the chemistry: it ran
15-20x faster on this comparison (Fortran over CVODE against single-threaded
Python), it ships a gas-giant extension (`photochem/extensions/gasgiants.py`,
`EvoAtmosphereGasGiant`, quench-based initial conditions and automatic boundary
management), and it ships a climate model (`clima`). The last is the real
argument: it is the only path to replacing the prescribed Guillot T(p) in
`run_lower.py` with a self-consistent one, and it makes A4 cheap enough to
iterate.

Practical state: the clone at `../photochem` carries both extensions, but the
package is **not importable from any Python on this machine** (checked
2026-08-19), so a switch starts with reinstalling it. The conversion path
(`photochem.utils.vulcan2yaml`) and its three traps -- the degenerate
`He <=> He` reaction, the grid-top requirement, and VULCAN's `atom_list` having
to match the network's atoms -- are already recorded in the comparison memo.

Neither code releases atomic metals (Na, Mg, Ca, Fe), so that part of the base
stays user-supplied either way.

**Recommended order if this is taken up:** fix the network before touching the
code. Decide the reference network, run both codes on it to separate the
residual 1.8x from the domain-truncation artifact, and only then move to
Photochem for the speed and for `clima`. Independently of which code wins,
`base.inp` should record *which code and which network* produced `q_H2_base` --
today that provenance exists only as a comment the reader ignores.

---

## 3. Problem B — oxygen in the wind

These are independent of everything above and can be done first.

**B1 — The O I 1302/1304/1306 triplet in the transit tool.**
Oxygen escaping from HD 209458 b was detected in exactly this multiplet
(Vidal-Madjar et al. 2004, ApJ 604, L69, "Detection of Oxygen and Carbon in the
Hydrodynamically Escaping Atmosphere of the Extrasolar Planet HD 209458b"; see
also Koskinen et al. 2010, ApJ 723, 116 on the same UV transits). EXHALE already
writes the O I number density in `Ion_species.txt`/`*_adv.txt`, and the transit
tool's metal-line table is a list of `(label, lambda0, f, A21, mass, n_lower,
R_instr)` tuples -- adding a species is adding rows plus a density accessor.
Buys: the one oxygen observable this code can actually compare against, on a
planet already in the paper set, using data already produced.
Costs: small in code. The physics caveats are real, though: the triplet is
optically thick, it sits in the FUV where interstellar absorption and
geocoronal emission contaminate the line, and the three components blend, so the
comparison has to be made the way the observers made it (a band-integrated
depth, not a line-center one -- the same trap already recorded for Mg II in the
Huang comparison).
Open physics question to settle first: O I 1025.76 A is nearly coincident with
H Ly-beta, so solar-system aeronomy treats Ly-beta pumping of O I as a real
source. Whether it matters for the 1302 multiplet in this regime needs a
literature check before the line is quoted.

**B2 — Charge-exchange completeness for oxygen.**
The code carries `O+ + H`, `O+ + He` and `Fe + O+`. Since O I is locked to H I
through the first of these, the O ionization balance is only as good as that
list. An audit against Huang et al. (2023) Table 4 groups B/C/D for any O2+
reaction, and a check of whether the reverse (endothermic) direction is treated
with the right Boltzmann factor, is cheap and would either close the question or
find something.
Note the code already gets one subtlety right, with the reasoning written at the
site: IP(O I) > IP(H I), so the barrier belongs on `O + H+`, not on `O+ + H`.

**B3 — Cooling: probably complete, one question open.**
O I and O II have CHIANTI fits, the O I ground term is in statistical
equilibrium, and the three fine-structure lines are trapped. The measured effect
of swapping the CNO coefficients to the CHIANTI set was 0.004 % on Mdot and
0.8 K on the temperature profile, so the oxygen coolant is not where the leverage
is. The open question is the same Ly-beta pumping as in B1: a pumped O I
population would radiate, and none of the current channels knows about it.

**B4 — Oxygen diffusion.**
`He_metal_diffusion` already lets each metal settle independently, but this has
never been exercised for oxygen specifically. A single numerical experiment
would show whether O separates enough from H to matter for B1. Cost: a run, not
a code change.

---

## 4. What the choice depends on

Three questions decide almost everything above.

1. **Is the target the paper or the code?** If the paper: B1 is the only option
   with an observable attached, and A0 remains a defensible stated assumption.
   If the code: A2 is the item that removes an external dependency.
2. **Should the base composition be an input or a result?** A0/A1/A4 keep it an
   input and differ only in how much crosses the boundary and how consistent it
   is. A2/A3 make it a result, at a cost that rises steeply between them.
3. **Is a VULCAN dependency acceptable in the production path?** If yes, A1 is
   the efficient answer and A2 loses much of its value. If no, A2 is the minimum
   that works.

## 5. Suggested order, if the answers are "paper first, VULCAN acceptable"

1. **B1** — O I 1302 in the transit tool. Highest return per line of code, uses
   data already on disk, and gives the paper a second element.
2. **B2** — the charge-exchange audit. Cheap, and it protects B1.
3. **A1** — widen the handoff. Removes the "one number crosses the boundary"
   limitation without new physics in the wind.
4. **A2** — only if a claim has to rest on a computed base partition. Settle the
   ownership question (Section 2.3) before starting.
5. **A3 / A4** — not now.

## 6. Validation gates, so the scope stays honest

- **B1**: reproduce the published HD 209458 b O I depth within its stated
  uncertainty, band-integrated the way the observation was reduced, with the
  interstellar treatment stated. Metals-off runs must skip the line, as the
  other metal lines already do.
- **B2**: every reaction in the audit either present with its source, or
  recorded as absent with a reason.
- **A1**: a base composition round-trip against the VULCAN profile it came from;
  the new keys default to absent, and a `base.inp` without them reproduces the
  current goldens byte-for-byte.
- **A2**: `q_H2` at 1 microbar reproduced against VULCAN for at least the
  hot-Uranus gate and the HD 209458 b run; the extended system reduces exactly to
  the present one at zero oxygen abundance; the whole regression matrix still
  passes with the new physics default off.

## 7. What not to do

- Do not add O-bearing molecules as new trace-metal ionization stages. They are
  not ionization stages, and the metal solver's assumptions (coronal-like rates,
  no three-body reactions) do not hold for them.
- Do not let two solvers own oxygen. Whichever way Section 2.3 is decided, there
  must be exactly one place that says how much oxygen exists in a cell.
- Do not take rates or cross sections from memory or from a review's summary
  table. Every reaction added here needs the same treatment R1-R23 got: a named
  source, its validity range, and the range written at the code site.

---

## 8. Related documents, and what each one owns

This memo is the decision document. It does not restate what the others
establish, and it should not be duplicated into them.

| Document | Owns |
|---|---|
| `docs/lower_atmosphere_coupling.md` | the three-tier coupling design; why a photochemical lower model matters (Lavvas 2014/2021); the Tier-3 `base.inp` sketch that A1 would finish |
| `docs/base_composition_handoff_plan.md` | the design of the implemented handoff: `q_H2_base`, `p_base`, and how the photochemical value replaces the equilibrium fit in the base particle count |
| `docs/vulcan_photochem_comparison.md` | the measured VULCAN-vs-Photochem comparison, its traps, and the runtime numbers quoted in Section 2.4 |
| `TO_BE_DONE.md` item (H) | the open blocker itself: why Lyman-Werner cannot close the gap and why oxygen is required |
| this file | the directions, their costs, the ownership trap, and the validation gates |

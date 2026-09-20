# D7c: the helium state the He <-> H charge-exchange pair reacts from

Item D7c of `docs/PLAN_20260918_rev2.md`, raised by the D7a report as the
first of the things it noticed outside its own scope.

## Verdict

The two He 2^3S systems handed the He <-> H pair of Huang et al. (2023),
ApJ 951, 123, Table 4 group B the SUMMED He I density, singlet plus
metastable. That is physically wrong: the tabulated rate is for ground-state
helium, and He(2^3S) + H+ is a different reaction with a different rate that
this code carries nowhere. Both systems now pass the ground singlet
`n_heiSI`, so the four ionization systems that carry helium and the
constrained chemical equilibrium all react the pair from one reservoir.

The correction removes the term k1 n(He 2^3S) n(H+) from the He + H+ rate,
so it moves every state in which the metastable is tracked and the pair is
on, which is the default. The movement measured is small because the
metastable is a part in 1e-5 of the neutral helium and because the forward
rate is barrier-suppressed at wind temperatures; the size is in section 4.

## 1. Which helium state the rate is for

READ from the journal PDF, `references/Huang_2023_ApJ_951_123.pdf`, Table 4
and section 2.4.2. The two group B rows are

| reactants | rate [cm^3 s^-1] | source |
|---|---|---|
| He + H+ | 1.75e-11 (T/300)^-0.75 exp(-12.75/T4) | Glover & Jappsen (2007) |
| He+ + H | 1.25e-15 (T/300)^0.25 | Glover & Jappsen (2007) |

with T4 = T/1e4 K. Section 2.4.2 says the reverse rates of the table are the
paper's own microscopic-balance estimates where only the exothermic direction
is measured; it says nothing that would make either row a metastable one, and
the source is a primordial-chemistry network, which carries no excited helium.

The barrier settles it. He + H+ -> He+ + H is endothermic by the difference
of the two ionization potentials, and for GROUND helium that is
24.587 - 13.598 = 10.989 eV, which is 10.989/8.6173e-5 = 1.2753e5 K, the
12.75/T4 of the table (READ from the table, the ionization potentials READ
from the standard values the code's own thresholds use). He(2^3S) lies
19.82 eV above the singlet, so its ionization potential is 4.77 eV and the
same collision is exothermic by 8.83 eV: no barrier, and a rate of a wholly
different size. Charging the metastable at the ground-state rate is therefore
not an approximation of that reaction, it is the wrong reaction.

## 2. What each call site passed

READ from the source. The line numbers of the two corrected files are their
entry-text ones; the corrected calls sit at `System_HeH_TR.f90:82` and
`System_HeH_TR_metals.f90:148`.

| call site | the helium density passed to the pair | is it the ground singlet |
|---|---|---|
| `System_HeH.f90:82` (and its Jacobian at `:135`) | `n_hei = (1-x2-x3) n_He` | yes; this system has no metastable, so its He I is the singlet |
| `System_HeH_metals.f90:171` (and its Jacobian at `:267`) | `n_hei = (1-x2-x3) n_He` | yes, same reason |
| `System_HeH_TR.f90:74` | `n_hei = (1-x2-x3) n_He`, the TOTAL He I | **no, this was the defect** |
| `System_HeH_TR_metals.f90:139` | `n_hei`, the total He I, its own line 95 saying so | **no, this was the defect** |
| `System_HeH_mol.f90:647` | `n_heiSI = n_hei - n_heiTR` | yes |
| `System_HeH_mol_metals.f90:268` | `n_heiSI` | yes |
| `constrained_chemical_equilibrium.f90:1210` | `sden(is_HeI_SI)` | yes |

The comment at `constrained_chemical_equilibrium.f90:1150` says "The He <-> H
pair below takes the ground singlet alone, exactly as in the fraction
systems". With the two triplet systems corrected that sentence is true; at
the entry text it was not.

Neither triplet system carries an analytic Jacobian: `ion_system_HeH_TR` and
`ion_system_HeH_TR_metals` are solved by `hybrd1` with a numerical one
(READ, `ionization_equilibrium.f90` lines 2590 and 2596), and the only
callers of `he_h_cx_jac` are `System_HeH` and `System_HeH_metals`, which the
dispatch reaches only when `thereis_HeITR` is false (READ, the same block).
There was therefore no Jacobian counterpart to correct; `he_h_cx_jac`'s
`n_hei` is the singlet in every configuration that reaches it.

## 3. The metastable's own charge exchange is absent, not double counted

`ion_residual_core::tr_triplet_row` (lines 88 to 102) carries the metastable's
photoionization, radiative recombination into the level, the collisional
excitation and de-excitation pair, the 2^3S -> 1^1S decay, its electron-impact
ionization and the total He(2^3S)+H ionization Q31 (Penning plus associative).
It carries no charge exchange with H+. So the entry text was not counting the
metastable's own reaction twice: it was applying the GROUND-state rate to the
metastable population on the summed He I row, while leaving the metastable
balance itself untouched by it, which does not correspond to any reaction.
After the correction the metastable takes no charge-exchange channel at all,
which is the state of the code's reaction set and is now said in the comment.
Adding the real He(2^3S) + H+ reaction is a separate physics item and needs a
published rate.

## 4. Impact

Control: a private build of the tree with only these two files at their entry
text. Measured: a private build of the delivered tree. Verified that no other
source moved between them (MEASURED: md5 over every `src/**/*.f90`, only the
two files differ). Both run with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`
on scratch copies.

### `wasp_he23off`, the triplet off

MEASURED: every data line of all ten output files is byte-identical; the only
difference in any file is the `run=` field of the provenance header. The
system this case solves is `System_HeH_metals`, which was not touched.

### `wasp_full`, the triplet system with metals, 400 steps

MEASURED, over the physical cells. log10(Mdot) = 7.20 g/s in both.

| column | max relative | at r [Rp] | max absolute [cm^-3] | at r [Rp] |
|---|---:|---:|---:|---:|
| `HeI` | 3.98e-6 | 1.2705 | 5.61e-2 | 1.0159 |
| `HeII` | 2.24e-10 | 1.2262 | 1.17e-2 | 1.0642 |
| `HeIII` | 1.58e-11 | 1.0670 | 9.42e-4 | 1.0653 |
| `HeITR` | 1.06e-7 | 1.2725 | 7.12e-8 | 1.0642 |
| `HI` | 9.87e-10 | 1.2648 | 6.30e-1 | 1.0159 |
| `HII` | 1.14e-11 | 1.2764 | 4.08e-2 | 1.0665 |
| `T` | 3.17e-11 | 1.2764 | | |
| `rho` | 1.14e-11 | 1.2764 | | |
| `CaI` | 2.14e-4 | 1.2764 | | |
| `FeI`, `NI` | 4.0e-5, 2.5e-5 | 1.458 | | |

The largest relative movements sit on columns whose density there is
negligible: `HeI` moves by 4e-6 at 1.2705 Rp, where n(He I) is 3.4e-7 cm^-3,
and `CaI` by 2e-4 at 1.2764 Rp, where n(Ca I) is 4.4e-14 cm^-3. In the dense
launch region the movement is at the 1e-11 level.

At the He 2^3S peak of this snapshot, r = 1.0588 Rp, the metastable is
7.22e-6 of the neutral helium (MEASURED) and the signed movements are
He I +1.12e-11, He II -9.43e-12, He III -8.47e-12, He(2^3S) -9.23e-12,
H II +2.05e-12 relative. The direction is the expected one: the pair no
longer ionizes the metastable, so a little more helium stays neutral.

The size is set by two factors, both MEASURED or READ: the term removed is
k1 n(2^3S) n(H+) against k1 n(1^1S) n(H+), a ratio of 7e-6 there, and k1
itself is 3.66e-18 cm^3 s^-1 at 1e4 K, barrier-suppressed by exp(-12.75)
against its 1.75e-11 prefactor.

### `.L26/fid_resolve`, the certified atomic fiducial, evaluate route

A scratch copy of the run directory with `Restart intent: stationary
evaluate` and the state's own `_IC` pair, run with both binaries.

MEASURED: every line of the certification block is identical, entry by entry
(active equations 6 of 30; hydrodynamic mass 1.202e-9 at cell 1, momentum
3.495e-12 at cell 484, energy 1.612e-8 at cell 24, elemental transport He/H
4.051e-1 at cell 2 with 6.679e-6 at cell 217 in the gated window, level
balance He 2^3S 1.624e-19, eliminated-species closure 7.449e-17; CERTIFIED
in the wind, `certified=T cert_reason=certified_in_wind`, original claim
REPRODUCED). The only difference anywhere in the log is the cooling-breakdown
arithmetic check, `max |sum(channels)/cool - 1|`, 3.22e-16 against 3.53e-16.

The composition and state files move at 1e-15 relative and below (MEASURED:
largest relative movement of any column of `Ion_species.txt` is 3.46e-15,
`HeII` at 1.3974 Rp; of `Hydro_ioniz.txt`, 1.96e-15 on the cooling column,
with `rho`, `v` and `p` byte-identical). The metastable is a larger fraction
of the neutral helium here than on WASP-121b, 1.25e-4 at its peak (MEASURED),
but the wind of this planet is at 3e3 to 5.5e3 K, where exp(-12.75/T4) puts
the forward rate ten orders of magnitude below its prefactor, so the term
removed is nothing at all.

## 5. Tests

`src/tests/charge_exchange_rows/he_h_pair_reservoir.f90`, a new row of the
D7a suite. It runs the residual of each of the four systems twice at one
state, with the pair on and off, and takes the difference, which is the
pair's contribution to each row and nothing else. See the suite `README.md`
for the state and the assertion list.

RED against the entry text, MEASURED: 16 of 30 assertions fail. The two
triplet systems' pair sources sit 1.22e-2 cm^-3 s^-1 away from
k2 n(He+) n(H0) - k1 n(He 1^1S) n(H+), which is exactly the
k1 n(He 2^3S) n(H+) the diagnostic line of the test prints; they disagree
with the molecular systems by the same amount; and the x(He 2^3S) column of
the pair's derivative is 0 where the singlet reservoir requires 0.1220403.

GREEN with the delivered text, MEASURED: the whole suite is 132 assertions,
0 failures, the four systems agreeing on the pair's sources to 0.0 exactly
and the derivative column at 1.2204032e-1 against the analytic
1.2204025e-1. The two other suites whose drivers reach these modules,
`ionization_imposed_fractions` (48 assertions) and `carrier_helium_inventory`
(15), also pass.

## 6. Noticed and not fixed

- **`System_implicit_adv_HeH_TR.f90:155`** passes `xhei = xheiS + xheiTR`,
  the summed He I fraction, to `he_h_cx_fvec_adv`. The comment just above it
  says the pair is charged to the singlet "the same one the
  equilibrium systems make"; with this item that sentence is now false in the
  other direction. The advection-correction system needs the same one-token
  change and its comment rewritten. Not in this item's files.
- **The metal + He reactions of Table 4 group C take the summed He I** in
  every metal system (`System_HeH_TR_metals.f90:131`,
  `System_HeH_mol_metals.f90:261`, `constrained_chemical_equilibrium.f90:1203`
  through `n_hei`), and `constrained_chemical_equilibrium.f90:1148` states
  that as a deliberate choice, "the He I reservoir the metal charge exchange
  reacts with (the metastable is a level inside the neutral stage)". Those
  rates are Si+ + He, C+ + He and O+ + He of Table 4 and are ground-state
  rates like group B, so the same argument applies to them; the difference is
  a part in 1e-5 of the rate and the group is off unless `cx_full 1` is set.
  Correcting it means one coordinated change across three files, not one, and
  is left for the item that owns them.

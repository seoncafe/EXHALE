# LHS 1140 b, item L36e: the ionization stages carried in an atomic gas

Fourth increment of item L36 of `docs/PLAN_20260917.md` section 10, which
carries stages B to F of `docs/lhs1140b_stationary_L12a_design_20260913.md`
as `docs/lhs1140b_stationary_L12b_derivation_20260916.md` sequenced them.
It follows `docs/lhs1140b_stationary_L36_20260918.md` (the stage flux as a
module), `_L36c_` (the H II row solved in a molecular run) and `_L36d_` (the
two helium stages beside it), and it closes the one thing those three left:
the ATOMIC configuration, which is the target of the whole L12 item.

Every number is MEASURED (run or computed here) unless marked READ.

---

## 0. Verdict, before the detail

1. **The atomic run carries the three stages, and the state it converges is
   CERTIFIED.** On the certified atomic fiducial
   `LHS1140b/models/.L26/fid_resolve` with `Ionization transport: True`
   through the partitioned stationary route, the outer iteration accepted at
   pass 56: every active equation within its own tolerance, the three stage
   rows at 8.620e-7 (H+), 7.147e-7 (He+) and 2.360e-6 (He++) against their
   1e-5 gate, and the two stage-sum entries at 2.161e-16 and 3.078e-16.
   Section 4.

2. **Twenty-five passes are not enough, and the reason is a relaxation
   front, not a defect of a row.** Over passes 1 to 25 the worst gated
   species row falls only from 3.51e-2 to 1.76e-2 while the cell that binds
   it marches OUTWARD from 217 to 491, about 3.3 cells a pass; the state is
   an ionization front relaxing from the local root of the restart towards
   the transported partition. The binding cell reaches the outer boundary at
   pass 26, and from there the measure falls monotonically: 1.75e-2 (26),
   1.27e-2 (40), 8.73e-3 (45), 2.03e-3 (50), 8.68e-6 (56, accepted).
   Section 4.

3. **The two closures of the He 10830 line AGREE once the stages are
   carried, and stage D is closed on the negative side.** On the certified
   state, whose three ionization fractions are within their own gate for the
   first time, the post-process moves the He 2^3S radial column by
   **0.127 per cent** (8.0846e10 against 8.0949e10 cm^-2) and the He I 10830
   equivalent width by 2.7 per cent (1.3678 against 1.4053 per cent A).
   The rule of L12a section 3 promotes the level when that column moves by
   more than 5 per cent with the three fractions at their gate; it does not,
   so **the He 2^3S level does NOT join the carried set**, which is the side
   its own Damkohler measurement supports (smallest Da 15.7 anywhere in the
   physical domain, READ, L12a section 1.1). Section 5.

4. **L36's factor 13.12 was a statement about the LOCAL-EQUILIBRIUM
   closure, and the carried stages remove it.** The He I 10830 equivalent
   width of the same planet: 20.3950 per cent A from the local-equilibrium
   solution, 1.5544 from its `_adv` profile (both reproduced here to four
   figures, READ from L36), and **1.3678 from the transported composition**,
   whose own `_adv` reads 1.4053. The transported solution and its
   post-process differ by a factor 1.027 where the local-equilibrium
   solution and its post-process differ by 13.12. Section 5.

5. **The gate is one predicate, named for what it says.**
   `transported_rows_exist()` (`ionization_equilibrium.f90`) is true where
   the transport-chemistry operator has a row to solve: the molecular
   carriers where the network is on and its transport selected, the three
   ionization stages on their own key and in any gas. It replaces
   `thereis_mol .and. carrier_transport` at thirteen sites of
   `src/EXHALE_main.f90`, at the operator's own three entry conditions, at
   the frozen background and at four sites of the certification. Section 1.

6. **The row source of an atomic gas is the sweep's own rows.**
   `carrier_source` takes the three stage sources from `heh_tr_rows`
   (`ion_residual_core`), the routine `System_HeH_TR` and `System_HeH` solve
   in that configuration, plus the He <-> H charge exchange of Huang et al.
   (2023) Table 4 group B. There is one spelling of each rate: a state
   stationary for the operator is stationary for the sweep. Section 2.

7. **The He <-> H charge exchange was missing from the stage rows and is a
   defect this item fixes.** Every ionization system that carries helium
   adds `he_h_cx_fvec` to rows (1) and (2) after the balance rows; the
   transport operator did not, so the transported partition was relaxed
   towards a root the sweep does not have. It is on by default. The one
   key-on regression fixture moves by 1.2479e-4 (n(He II)) because of it.
   Section 3.

8. **Key-off byte identity: MEASURED** on `wasp_full` (400 steps),
   `mol_diffusion` and `mol_carrier` (12000 each): the data rows of
   `Hydro_ioniz.txt` and `Ion_species.txt` are byte-identical to the control
   build on all three. No golden was refreshed, `EXHALE.x` was never touched
   and nothing was written under `backup/regression/`. Section 3.

9. **Four defects outside the stage rows were found and fixed**, and they
   are reported as their own item: an out-of-bounds write in the update-map
   diagnostic, a gas particle density of zero in every atomic run, a row
   record whose scale was not the scale the acceptance divides by, and a
   restart classification that made two of its own notes unreachable.
   Section 7.

---

## 1. The predicate, and the sites it replaced

The transport-chemistry operator solves TWO independent sets of rows. The
molecular carriers -- H2, and OH, H2O and CO with the oxygen cycle -- exist
where the molecular network is on and `Molecular carrier transport` selects
them. The three ionization stages x(H II) per hydrogen nucleus and
x(He II), x(He III) per helium nucleus are stages of an ELEMENT: their flux
is that element's own nucleus face flux and a hydrogen and helium mixture
has them whether or not it has molecules.

```
    transported_rows_exist() = (thereis_mol .and. carrier_transport)
                               .or. ionization_transport
```

in `ionization_equilibrium`, which is the lowest module every consumer
already reads. What it now conditions:

| where | what it opens |
|---|---|
| `EXHALE_main.f90`, 13 sites | the operator-split step of the marching loop (line 2499 through `photochemical_transport_step`), the fixed-wind relaxation of the partitioned stationary route, the pass cap and `species_alternated` of the outer iteration, and its ten diagnostics and progress-control sites |
| `diffusive_photochemistry.f90` | the three entry conditions of the operator itself (`carrier_transport_interval`, `carrier_steady_residual`, `relax_photochemical_composition`) and the advected-carrier registry |
| `ionization_equilibrium.f90` | `bg_cell(j) = ieq_cell` and `bg_ready`, and the four keep/adopt/install routines of the stationary solver's background |
| `certification.f90` | the carrier rows of a state, their entries, the carrier history and the step verdict |

Two sites keep `thereis_mol .and. carrier_transport` and say why: the
update-map record of `EXHALE_main.f90` line 6023, which follows the H2 row
alone, and `carrier_h2_chemical_root`, which is a statement about H2.

`carrier_set_init` now states the same split: the molecular carriers are
solved rows only where the network exists AND its transport is selected,
the stages on their own key. Before this item `carrier_solved(ic_H2)` was
true in every configuration, which in an atomic gas would have put an H2
row into a gas with no H2.

**What the frozen background needed.** `bg_cell` is the coefficient state
each cell was left in by the last sweep, and it was filled inside
`if (thereis_mol)`. Three things moved out of that block, unchanged in
value: the He 2^3S coefficients (which are now zeroed in an `else` branch
rather than inside the molecular one, so an atomic cell state carries them
either way), the cell temperature and gas particle density, and the imposed
ionization-stage fractions with their three flags. The molecular part of
the cell state stays where it was.

## 2. The row source of an atomic gas

`carrier_source` evaluates the chemistry of every row on the frozen
background. Which balance writes the H and He ionization rows is now the
gas, and in each case it is the balance the LOCAL SWEEP of that same gas
solves:

| gas | rows (1), (2), (3) from | charge exchange |
|---|---|---|
| molecular | `mol_heh_rows` (`System_HeH_mol`), taken whole, with its molecular sinks R17/R23, the He+ + CO channel and the H2 photoionization branches | `he_h_cx_fvec(..., n_heiSI, ..., +1)`, the orientation and the reservoir `System_HeH_mol` passes |
| atomic | `heh_tr_rows` (`ion_residual_core`): photoionization with its secondaries, collisional ionization, radiative and dielectronic recombination, the He 2^3S channels and the Penning proton source | `he_h_cx_fvec(..., n_hei, ..., -1)`, the orientation and the reservoir `System_HeH_TR` passes |

`heh_tr_rows` is the general atomic form: with the metastable untracked its
coefficients are zero and its four rows reduce term by term to the rows of
`System_HeH`, so the atomic branch is ONE branch and not two. Its row (2)
is the summed He I balance written He-I-gain positive and its row (3) the
He++ balance, so the stage sources are

```
    src(H II)  =  f(1) ,
    src(He II) = -f(2) - f(3) ,
    src(He III) = f(3) ,
```

which is the statement that a helium nucleus leaving He I enters He+ and
one leaving He+ enters He++; the three stage sources of the element sum to
zero, which is what makes the transport of the partition conserve helium
nuclei. Row (4), the He 2^3S balance, is not a stage of the partition --
the metastable is a sublevel inside He I -- and is not read.

**The He 2^3S level stays with the sweep**, for the Damkohler reason of
L12a section 3 and now with the rule of section 5 evaluated.

## 3. Key-off byte identity, and the one key-on case that moves

Control build `EXHALE_L36e_ctl.x`, md5 `cf3d2941d1ebe7864f14bd411d5bdfb9`,
a build of the ENTRY TEXT of every file of this item, which is the binary
the L36d increment delivered. Measured build `EXHALE_L36e.x`, md5
`6a36529386102b37a38cb40fde0a611c`. Scratch case copies,
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, identical bounds on the two
builds.

| case | key | steps | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|---|
| `wasp_full` | off | 400 | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_diffusion` | off | 12000 (its pinned count) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_carrier` | off | 12000 (its pinned count) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `carrier_model_a_newton` | **on** | 3 outer passes, `EXHALE_JFNK_MAXIT=5` | 4.2725e-05 | 1.2479e-04 |

`mol_carrier` and `mol_diffusion` run the transport-chemistry operator
itself with the ionization key off, so they are the cases that reach every
changed gate with `carrier_solved(ic_Hp)` false; `wasp_full` is atomic with
helium, the He 2^3S level and metals, so it reaches the changed cell state
and the changed `n_tot` of section 7 with no molecules at all. Their
printed step lines are identical over all 400 and all 12000 steps.

**The key-on fixture moves, and the movement is the charge exchange.**
Column by column on `carrier_model_a_newton`: n(He II) 1.2479e-4 (cell 76,
1.256321e2 to 1.256164e2), n(H I) 2.5916e-5, n(He III) 1.6337e-5, n(HeH+)
1.1289e-5, and 4.2725e-5 on the velocity. Every gate this item rewrote
evaluates identically in that configuration (`thereis_mol` and
`carrier_transport` are both true there), and the key-off molecular cases
above show the gates and the cell-state restructure move nothing at all, so
what is left is the one term the stage rows gained: the He <-> H charge
exchange of Huang et al. (2023) Table 4 group B, which enters rows (1) and
(2) alone and whose largest movers are exactly helium and neutral hydrogen.

That fixture's golden joins the three `hp_*` goldens L36c and L36d reported
as stale. **No golden was refreshed here.**

## 4. The atomic fiducial with the key on

`LHS1140b/models/.L26/fid_resolve` (LHS 1140 b, He/H = 2.13, `He_Kzz = 1e9`,
`He_diffusion: True`, He 2^3S on, metal-free, 500 cells to 29.03 R_p), on a
scratch copy, with `Ionization transport: True` and
`Restart option change: iontrans` added, `Restart intent: stationary`,
`Coupled carrier solve: False`, `EXHALE_PTC_DTAU0=1.0`, `OMP_NUM_THREADS=8`.

**Twenty-five passes, the bound the item was briefed with: NOT certified.**

| pass | worst gated species row | cell | row |
|---|---|---|---|
| 1 | 3.51e-02 | 217 | `carrier balance He++` |
| 2 | 2.05e-02 | 416 | `carrier balance He++` |
| 10 | 1.91e-02 | 444 | `carrier balance He++` |
| 20 | 1.80e-02 | 475 | `carrier balance He++` |
| 25 | 1.76e-02 | 491 | `carrier balance He++` |

every pass at `hydro info=0`. The final certification of that state:

| entry | max | cell | gated at r >= 1.20 | verdict |
|---|---|---|---|---|
| `carrier balance H+` | 7.960e-02 | 1 | 9.209e-03 at 445 | ABOVE |
| `carrier balance He+` | 6.622e-03 | 464 | 6.622e-03 at 464 | ABOVE |
| `carrier balance He++` | 1.760e-02 | 491 | 1.760e-02 at 491 | ABOVE |
| `ionization stage nucleus sum H` | 2.196e-16 | 233 | reported | within |
| `ionization stage nucleus sum He` | 4.022e-16 | 226 | reported | within |
| `elemental transport He/H partition` | 2.256e-05 | 339 | 2.256e-05 at 339 | ABOVE |
| hydrodynamic mass / momentum / energy | 5.719e-10 / 2.197e-13 / 1.318e-08 | | | all within |

**What the refusal is.** The binding cell marches outward at about 3.3
cells a pass, which is an ionization relaxation FRONT and not a stalled
row. Read term by term from the row record
(`EXHALE_CARRIER_ROW_TERMS=1`) of that state, every rate a volumetric rate
in cm^-3 s^-1:

| row | cell | r [R_p] | stage flux divergence | net chemical source | residual | row scale | measure |
|---|---|---|---|---|---|---|---|
| H+ | 445 | 11.747 | 8.402e-02 | 2.880e-03 | 8.114e-02 | 8.811e+00 | 9.209e-03 |
| He+ | 464 | 15.967 | 3.119e-02 | 2.412e-03 | 2.878e-02 | 4.347e+00 | 6.622e-03 |
| He++ | 491 | 24.963 | 9.341e-04 | 5.767e-05 | 8.764e-04 | 4.980e-02 | 1.760e-02 |

The transport divergence stands 13 to 29 times the net chemical source and
the residual IS that divergence: the front has not yet reached those cells,
so nothing there balances. This is the opposite of the regime L36d measured
on `carrier_model_a_newton`, a hot Uranus whose gated cell sits at 1.20 R_p
and whose He+ chemistry stands 68 times its stage flux divergence.

**Ninety passes: ACCEPTED at pass 56, CERTIFIED.** The binding cell reaches
the outer boundary at pass 26 and the measure then falls monotonically:

| pass | worst gated species row | row |
|---|---|---|
| 26 | 1.75e-02 | `carrier balance He++` |
| 40 | 1.27e-02 | `carrier balance He++` |
| 45 | 8.73e-03 | `carrier balance He++` |
| 50 | 2.03e-03 | `carrier balance H+` |
| 53 | 7.11e-05 | `elemental transport He/H partition` |
| 56 | 8.68e-06 | `elemental transport He/H partition` -- ACCEPTED |

| entry | max | cell | verdict |
|---|---|---|---|
| `carrier balance H+` | 8.620e-07 | 150 | within |
| `carrier balance He+` | 7.147e-07 | 177 | within |
| `carrier balance He++` | 2.360e-06 | 183 | within |
| `ionization stage nucleus sum H` | 2.161e-16 | 294 | within (reported) |
| `ionization stage nucleus sum He` | 3.078e-16 | 166 | within (reported) |
| `elemental transport He/H partition` | 8.679e-06 | 500 | within |
| `level balance He 2^3S` | 7.263e-21 | 442 | within |
| `eliminated-species closure System_HeH_TR` | 5.043e-15 | 495 | within |
| hydrodynamic mass / momentum / energy | 5.670e-10 / 9.788e-12 / 1.045e-08 | | all within |

`CERTIFIED: every active equation was evaluated and is within its
tolerance`, and the stationary solve returned `info = 0`.

**The rows of the certified state, term by term**, at the cell that binds
each of them:

| row | cell | r [R_p] | stage flux divergence | net chemical source | residual | row scale | measure |
|---|---|---|---|---|---|---|---|
| H+ | 500 | 29.031 | 2.792e-04 | 2.792e-04 | 1.281e-08 | 2.269e-01 | 5.646e-08 |
| He+ | 217 | 1.201 | 1.538e+01 | 1.538e+01 | -2.403e-05 | 7.048e+02 | 3.409e-08 |
| He++ | 498 | 28.073 | 1.729e-05 | 1.729e-05 | 3.076e-10 | 4.158e-03 | 7.398e-08 |

Transport against reaction to eight digits, which is what a stationary
stage row is.

**The row scale, which is the question L36d raised and which is the user's
and not this item's.** The scale a stage row is judged on is
`row_terms_phys`: the two face fluxes by their own MAGNITUDES, the
advective magnitudes, the reaction terms and the 1e-20 free-element floor.
At the He+ row's binding cell of the certified state it is 704.8 cm^-3
s^-1 against a residual of 2.4e-5, so the row certifies with four decades
to spare; on this configuration the scale does not bite. Where it bit is
L36d's hot Uranus, and what it would take to change it -- measuring a fast
stage against its production PLUS its loss instead of against the net --
needs the production/loss split at the one place the terms are written and
is a change to a shared acceptance interface with its own anchoring. **No
tolerance was set or changed by this item.**

## 5. The line, and stage D

`EXHALE_transit.py` with `EXHALE_TRANSIT_STATE` selecting the input pair
(L36 built the selector). The equivalent width below is
`int (1 - T) dlambda` over the He I 10830 window on the rotation- and
instrument-convolved curve, in per cent A; it reproduces L36's two numbers
to four figures, which is what fixes the definition.

| state | EW [per cent A] | max depth |
|---|---|---|
| local-equilibrium solution (key off) | 20.3950 | 66.1194 |
| its `_adv` profile | 1.5544 | 5.1959 |
| transported solution, 25 passes (not certified) | 2.9823 | 10.4172 |
| **transported solution, certified** | **1.3678** | **4.5745** |
| its own `_adv` profile | 1.4053 | 4.7338 |

**The factor 13.12 is a property of the local-equilibrium closure.** With
the stages carried the solved state and its post-process agree to a factor
1.027, and the line the transported composition gives is 1.3678 per cent A
against 20.3950 from the same run's equations with the stages left local.

The composition behind it, MEASURED on the certified state against the
key-off one: x(H II) reaches 0.0583 against 0.6657, and
x(He II) + x(He III) reaches 0.1888 against 0.9117. That is the
over-ionization of a local root in a wind whose ionization time is longer
than its flow time (`docs/k22_electron_density_excess.md` sec. 7, READ),
removed.

**Stage D, evaluated for the first time and CLOSED.** The rule of L12a
section 3: the He 2^3S level joins the carried set if the post-process
moves its column by more than 5 per cent WHEN the three ionization
fractions already agree to their own gate. The certified state is the first
state on which that precondition holds. MEASURED on it, He 2^3S radial
column over 1 to 29.03 R_p:

| closure | N(He 2^3S) [cm^-2] |
|---|---|
| the solved state | 8.0846e+10 |
| its `_adv` profile | 8.0949e+10 |

a difference of **0.127 per cent**, twenty times below the 5 per cent
threshold. **The level is not promoted.** For the record, the same column
of the local-equilibrium closure is 4.9257e+11 against its `_adv`
8.0109e+10, a factor 6.15 (L36 measured 6.15 and this item reproduces it),
which is what the rule would have read had it been evaluated on a state
whose fractions do not meet their gate -- and is why the precondition is
part of the rule.

## 6. Where the approximation is valid, restated beside the result

The stage rows carry ONE velocity: an ion and its neutral move together.
MEASURED by L36 section 3 on this same fiducial (READ):

- the ambipolar drift of an ion against its own neutral stays below 0.1 of
  the bulk velocity only inside **1.220 R_p** (hydrogen) and **1.243 R_p**
  (helium) without the resonant charge-exchange channel, 1.231 and 1.273
  R_p with it;
- over `r >= 1.2 R_p` it reaches **0.625** (H) and **0.430** (He) without
  that channel, 0.518 and 0.336 with it;
- the momentum transfer itself is fast: `t_coll/t_flow < 0.1` inside
  14.51 R_p (H) and 8.25 R_p (He), so the stages are collisionally locked
  in TIME SCALE and still drift, the drift being set by the ambipolar field
  against the friction and not by the friction alone;
- above about 13 R_p the continuum equation is itself unvalidated by the
  run's own Knudsen measure on this state (max Kn 0.88 in the heating and
  acceleration region).

So the certified state of section 4 is a solution of an equation whose
one-velocity closure is exact only inside about 1.2 R_p, and whose
neglected drift is of order 0.02 to 0.6 of the advective stage flux over
the range the stages are carried. An ambipolar stage drift is a separate
physical term and is not in this item; the approximation is carried at the
site with its measured size, in the module header of
`ionization_stage_transport.f90`.

## 7. Defects found and fixed, outside the stage rows

Each was confirmed in the source, fixed, and its fix verified on the path
it touches.

1. **An out-of-bounds write in the update-map diagnostic**
   (`EXHALE_main.f90`). `res_um_all` and `terms_um_all` were allocated
   `(1:N,4)` and handed to `carrier_steady_residual`, whose dummies are
   `(1:N,n_carrier_max)` with `n_carrier_max = 7`: the routine zeroes the
   whole dummy, which wrote 3N doubles past the end of both arrays. It
   became wrong when the carrier set grew past four. FIXED: both are
   allocated by `n_carrier_max`.

2. **The gas particle density was zero in every atomic run**
   (`ionization_equilibrium.f90`). `calc_ntot` was called only under
   `thereis_mol` and `n_tot` was set to zero otherwise, so an atomic run's
   frozen background carried `ntot = 0`, the pressure the sweep's retry
   forms was the electron pressure alone, and the transport operator's
   `wfac = rho n0 / max(n_tot, 1e-99)` would have been of order 1e99 for
   any row that read it. The stage rows do not read it -- they cycle past
   the mixing-ratio face flux and past `carrier_face_coefficients`, which
   is why the numbers of section 4 do not move -- but the record's own
   fraction column read 1e307 and said so. FIXED: the particle density is a
   property of the gas and is formed in every configuration from the one
   stoichiometric sum. `wasp_full` (atomic, helium, He 2^3S, metals) is
   BYTE-IDENTICAL across the change, so no existing run moves.

3. **The carrier row record did not print the scale the acceptance divides
   by** (`diffusive_photochemistry.f90`). The writer rebuilt a scale from
   the flux DIVERGENCE, while `carrier_residual` builds it from the two
   face fluxes by their magnitudes; in a smooth wind the two faces cancel
   and the difference is two decades. MEASURED at the H+ binding cell of
   the 25-pass state: the record read 0.934 where the certification read
   9.209e-3 for the same row of the same state. FIXED: the record carries
   `row_terms_phys` itself and the measure formed from it, and now
   reproduces the certification's three gated numbers exactly
   (9.2093e-03, 6.6219e-03, 1.7600e-02 at cells 445, 464, 491).

4. **`heh_tr_rows` returns no production/loss split, and the proton record
   invented none.** In an atomic gas the record's H+ production and loss
   were both zero while its source was not. FIXED: the proton carries its
   split where `mol_heh_rows` hands one out and its NET on the side its own
   sign falls where no split exists, which is what the two helium rows
   already did, so `production - loss` is the row's own source for every
   row in every gas.

5. **`iontrans` was classified as a token that changes how many unknowns a
   state has** (`load_IC.f90`), which refused every restart across the key
   and made the two notes beside it unreachable. It is not one: n(H II),
   n(He II) and n(He III) have a column in every state file this code
   writes, and what the key changes is whether those columns are each
   cell's local root or a partition the flow carried. FIXED: it is a token
   a `Restart option change:` line may name, which is what the option
   ladder is for, and what the certified run of section 4 uses. The four
   tokens that DO change the columns (`metals`, `mol`, `oxychem`,
   `carrier`) are unchanged.

6. **Three statements that named the proton a molecular carrier** are
   corrected in place: the `Molecular carrier transport` key documentation
   and the `Ionization transport` key documentation of `input_read.f90`,
   and the header of `carrier_transport_inert_report.sh`.

## 8. Tests

`src/tests/grid_and_gates/ionization_transport_atomic.sh` is new, five
rows, registered as `iontrans_atomic`:

| row | measured | RED on the entry text |
|---|---|---|
| `ionization_transport_atomic_runs` | rc=0, no refusal | rc=1, refused in `input_read` |
| `ionization_transport_atomic_setup_names_rows` | yes | no |
| `ionization_transport_atomic_stage_sum` | 2 entries evaluated | 0 |
| `ionization_transport_atomic_rows_measured` | 3 rows evaluated | 0 |
| `ionization_transport_molecular_local_refused` | rc=1, refused by name | rc=1 (green then as now) |

The four RED rows are RED on the control build for the reason the item
exists, and the fifth guards the refusal that is KEPT.

`grid_and_gates`' four `route_change_*` rows, which L36d left failing,
are GREEN. Their block tests the ROUTE token and nothing about the
ionization stages, and it built its rungs from
`backup/regression/carrier_model_a_newton`, which carries
`Ionization transport: True` and whose rung B forces
`Coupled carrier solve: True` -- a configuration the refusal L36c added
aborts. The block now runs on `backup/regression/carrier_elem_newton`,
which admits both routes and carries no ionization key. **The refusal was
not lifted**: the four rows are GREEN on the control build too with the
corrected fixture, which is how one can tell that the fixture was the
defect and not the binary.

Every suite whose driver links a changed object, on the measured build:

| suite | result |
|---|---|
| `ionization_stage_flux` | 38 passed, 0 failed |
| `ionization_imposed_fractions` | 48 passed, 0 failed |
| `carrier_reference_scales` | 16 passed, 0 failed |
| `steady_species_rows` (with the reloads, `EXHALE_SPECIES_EXE`) | 201 passed, 0 failed |
| `element_operator` | 40 passed, 0 failed |
| `species_face_flux` (with the whole-binary rows) | 15 passed, 0 failed |
| `carrier_boundary_jacobian` | 12 passed, 0 failed |
| `carrier_outer_boundary` | 16 passed, 0 failed |
| `carrier_constraint_attribution` | 10 passed, 0 failed |
| `carrier_returned_state_acceptance` | 36 passed, 0 failed |
| `certification` | 84 passed, 0 failed |
| `attempted_step` | 70 passed, 0 failed |
| `adv_static_limit` | 53 passed, 0 failed |
| `coupled_source_step` | 30 passed, 0 failed |
| `molecular_seed` | 25 passed, 0 failed |
| `run_mode` | 31 passed, 0 failed |
| `state_mapper` | 76 passed, 0 failed |
| `acceptance_classes` | 22 passed, 0 failed |
| `constrained_network_layout` | 8 passed, 0 failed |
| `steady_completion_flag` | 3 passed, 0 failed (on the certified run's log) |
| `steady_selfconsistent_residual` | 7 passed, 0 failed (same log) |
| `grid_and_gates`: `restart_option_change`, `iontrans_atomic`, `carrier_transport_inert` | 29 passed, 0 failed |
| `carrier_retry` | 151 passed, 1 RED, section 9 item 1 |
| `coupled_block_jacobian` | 2 RED on BOTH builds, section 9 item 2 |

**Two fixtures had to state the configuration they mean.**
`carrier_returned_state_acceptance` said "the oxygen cycle on, so H2, OH,
H2O and CO all carry an equation" and set `thereis_oxychem` alone; a
carrier is a row of the operator where the gas HAS it and the run
TRANSPORTS it, so the fixture now sets `thereis_mol` and
`carrier_transport` too. `carrier_retry` has three rows about a
configuration with nothing to transport; with the ionization key left on
by its own setup those configurations now DO have something to transport,
so the three rows hold the key off and put it back.

## 9. Noticed outside this item

1. **`carrier_retry`'s knife-edge row is unchanged by this item.**
   `one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
   reads 1.100521341448987e-06 against `ieq_res_tol` = 1e-6 on the measured
   build AND on the control build with the same fixture, so the code change
   does not move it; what moved it from L36d's 1.044879e-06 is the fixture
   correction of section 8. The binding species is still HeH+ at 1.3e-12 of
   the gas in the base cell, the chain still settles at a third sweep
   (L36d section 8, READ), and the tolerance and the 1e-20 presence
   threshold were NOT touched.

2. **`coupled_block_jacobian` has two RED rows on BOTH builds**, MEASURED
   here with the control objects:
   `the_assembled_action_reproduces_the_central_difference_carrier`
   8.9540842e-01 against 1e-2 and `the_species_rows_scale_is_the_certifications`
   3_of_13 against 0_of_13. Pre-existing; that suite is about the coupled
   route, which the ionization key refuses.

3. **The metal charge exchange is still absent from the stage row
   sources.** `cx_add_to_fvec` writes group A (metal + H/H+) into row 1 and
   groups C and D into rows 2 and 3 of every metals system, and
   `carrier_source` does not call it, so in a run with metals the stage
   sources and the sweep's rows still differ by those terms. The atomic
   fiducial is metal-free, so it does not enter anything measured here.
   Adding it needs the frozen metal densities per stage and the layout
   `cx_metal_base` of the system being judged, neither of which the
   transport operator holds today; REPORTED and not done.

4. **Four goldens of key-on cases are stale**: `hp_front`, `hp_zero_seed`
   and `hp_trace_seed` (L36c and L36d, READ) and now
   `carrier_model_a_newton` by 1.2e-4. The refresh they need is one refresh
   at the close of the series and is not taken here.

5. **Twenty-five outer passes is not a bound at which this configuration
   can be judged.** The ionization relaxation front crosses the domain at
   about 3.3 cells a pass from the restart's local root, so a 500-cell
   column needs of order 30 passes before any gated row can begin to fall
   and 56 before it certifies. A campaign that judges a transported
   ionization state at 25 passes will read a front and not a solution.

## 10. Where the L12 stages stand after this increment

| stage (L12b's sequence) | state |
|---|---|
| B, the stage flux on the shared element flux | landed with L36, carried for both elements since L36d, and now carried in an atomic gas |
| C, the rows solved in the run | CLOSED. Hydrogen with L36c, helium with L36d, the atomic configuration here, and the atomic state certifies |
| D, the He 2^3S promotion rule | CLOSED, on the negative side: on a state whose three fractions meet their gate the post-process moves the He 2^3S column by 0.127 per cent against the rule's 5 per cent, so the level is not promoted (section 5) |
| E, the certification entries and the tolerance anchoring | CLOSED with L36d; this item adds no entry and changes no tolerance |
| F, the `_adv` retirement as the source of the carried fractions | CLOSED. The selector landed with L36; the result it was for is here: with the stages carried the solved state and its post-process give He I 10830 equivalent widths of 1.3678 and 1.4053 per cent A, a factor 1.027, against the 13.12 of the local-equilibrium closure |

# Item D3b: the stage sum entries made gates at their floating-point bound

The acceptance-interface half of item D3 of `docs/PLAN_20260918_rev2.md`,
on the measurement of `docs/lhs1140b_stationary_D3a_20260918.md`. The user's
decision (the D3b row of that plan's "Decisions made" table, 2026-09-18) is
carried out here and nothing else in the interface moves: the He II gate and
its denominator stay exactly as they were, no other tolerance changes, and
the only entries whose verdict is new are the two named below.

Every number is MEASURED (run or computed here) unless marked READ.

| binary | md5 | what it is |
|---|---|---|
| `EXHALE_D3b_ctl.x` | `1abf66c7024f0c2a066814e9749f2f21` | the tree with this item's three production files at their entry text |
| `EXHALE_D3b.x` | `147b0fbc6015c5fd181eaf81d23c77a4` | the same tree with them as delivered |

A source manifest of all `src/**/*.f90` was taken before the first build and
again after the second: between them no source but this item's three
production files and its one test file changed, so the pair differs by this
item alone. Both are private builds (`OBJDIR=build_D3b`, `build_D3b_ctl`) of
the conda-forge gfortran on PATH, deleted at the end; the tree's `build/` and
`EXHALE.x` were never touched and nothing was written under
`backup/regression/`.

---

## 0. Verdict

1. **The two stage sum entries gate**, at the floating-point bound of the
   identity they measure. They were `evaluated` with `tol = 0` and
   `within_tol = .true.` written unconditionally, which is a row that could
   not refuse anything whatever it measured.
2. **The tolerance is `2 (5 nk + 2) eps g`** of D3a section 6, with `g` the
   ratio of magnitude sums MEASURED at the face the measure was taken at, so
   each row is held to the rounding of its own face. It is formed in one
   place, `ionization_stage_sum_rounding_bound` beside the face flux it
   bounds, and the entry reads it from there. On the states this code
   produces `g` is 1 to five decimals and the tolerance is **3.109e-15 for
   hydrogen and 5.329e-15 for helium**. Where the operator reports no `g` the
   bound falls back on the algebraic CEILING of `g` (3/2 for one carried
   stage, 2 for two), which is 4.663e-15 and 1.066e-14.
3. **The three states D3a measured all certify on these rows, at 0.06 to
   0.08 of their own tolerance**, and every other entry of their reports is
   unchanged to every printed digit against the control binary. The
   constructions that break the identity on purpose stand ten decades above
   (1.5e-4 against 5.329e-15 is 2.8e10).
4. A non-finite measure now refuses the row. Under the entry text it
   certified, because `within_tol` was written `.true.` without reading the
   measure at all.
5. The entry's `scale:` line now names the `g` the tolerance was taken at, so
   a reader of a report can see which face's rounding the row was held to.

---

## 1. What the entries were, and what they are

The identity is

```
   sum_k F_k(f) + F_close(f) = N_el(f)        at every face f,
```

the stage fluxes of one element summing over ALL of its stages, the closing
neutral one included, to that element's own nucleus flux
(`ionization_stage_transport`, equation 2). It holds by construction, so
what stands in it is the rounding of the sums.

The entry text (READ, `certification.f90` as this item entered it) wrote

```
   rep%e(rep%n)%tol         = 0.0d0
   rep%e(rep%n)%within_tol  = .true.
```

with the reason `reported, not gating: no anchored tolerance`. D3a took that
anchor. The entry now reads

```
   rep%e(rep%n)%tol         = ionization_stage_sum_rounding_bound(nk, g)
   rep%e(rep%n)%finite      = finite_real(dmax)
   rep%e(rep%n)%within_tol  = rep%e(rep%n)%finite .and. (dmax .le. tol)
```

with the reason `gating at the floating-point bound of an algebraic
identity (anchor 8)`, which names anchor (8) of
`docs/certification_tolerance_anchoring_20260910.md`. What is unchanged: one
entry per element, the measure read from the stationary evaluation that
formed the stage fluxes rather than re-formed, the `unavailable` status for
an element none of whose stages is solved, and the entry's absence from a
run that transports no stage at all.

---

## 2. The bound, and where it is formed

`ionization_stage_sum_rounding_bound(nk [, g])` in
`src/modules/functions/ionization_stage_transport.f90` returns
`2 (5 nk + 2) eps g`, next to `ionization_stage_face_flux`, whose rounding it
counts. `g` is passed in, measured at the face the identity was measured at;
without it the function returns the bound at
`ionization_stage_sum_scale_ratio_limit(nk)`, the algebraic CEILING of `S/S'`
derived in the same file and in anchor (8):

| element | nk | tolerance at the measured `g = 1` | ceiling of `g` | tolerance at the ceiling |
|---|---|---|---|---|
| hydrogen | 1 | **3.1086244689504383e-15** | 1.5 | 4.6629367034256575e-15 |
| helium | 2 | **5.3290705182007514e-15** | 2.0 | 1.0658141036401503e-14 |

The `g = 1` column is D3a section 6's table, reproduced by the function, so
the two documents describe one number.

**Where `g` comes from.** `ionization_stage_face_flux` returns the stage eddy
terms `E_k` through an optional output, so the eddy term keeps the one
spelling the module header asks of the stage flux;
`ionization_stage_nucleus_sum` forms `S` and `S'` from them at each face and
returns `S/S'` at the face that carries the maximum, beside the maximum
itself; `ionization_stage_sum_measure` carries it to the certification with
the measure and the face. A negative value means no face of the column
carried a nonzero scale, and the entry then takes the ceiling.

---

## 3. The three production states

Each is a state D3a measured, re-entered on a scratch copy with
`Restart intent: stationary evaluate`, single threaded, and run twice: once
with the control binary and once with the measured one, from the same
restart pair.

| state | run directory of D3a | `... sum H` | `g` | of tol | `... sum He` | `g` | of tol |
|---|---|---|---|---|---|---|---|
| hot Uranus, transported: `carrier_model_a_newton`, `Ionization transport: True`, 25 outer passes | `D3a/runs/cma_eval_meas` | 2.203e-16 at cell 221 | 1.00000 | 0.071 | 4.195e-16 at cell 264 | 1.00000 | 0.079 |
| hot Uranus, its local root: the same fixture with the key off, solved to `info = 0`, evaluated with the key on | `D3a/runs/cma_eval_ctl` | 2.144e-16 at cell 249 | 1.00000 | 0.069 | 3.049e-16 at cell 330 | 1.00000 | 0.057 |
| atomic fiducial, certified: the L36e 90-pass state of `LHS1140b/models/.L26/fid_resolve`, key on | `D3a/runs/fid_eval_ch` | 2.161e-16 at cell 294 | 1.00000 | 0.069 | 4.106e-16 at cell 229 | 1.00009 | 0.077 |

All six readings are `within`, at tolerances 3.109e-15 (hydrogen) and
5.329e-15 (helium) times their own `g`. Each measure reproduces D3a's own
table for the same state to every digit it quotes, and so does each ratio,
because D3a's table is written at `g = 1` and these faces carry `g = 1`.

**Every one of the six worst faces carries `g = 1` to five decimals.** The
identity's largest rounding sits where the advective nucleus flux is largest,
and there the stage eddy term is negligible beside it, so `S` and `S'` differ
in the fifth decimal at most. The ceiling of `g` is therefore never the
operative number on a state this code produces.

**The rest of each report does not move.** The whole run log of the control
and of the measured build differs in exactly the two `tol=` fields of these
two entries and in nothing else: 4 differing lines on the transported state,
4 on the fiducial, 8 on the local root (which prints its report twice). The
four output files of each pair are identical row for row, the only
difference in any of them being the run timestamp of a comment header.

The third state is the one the brief asked for if it was on disk from L36e
with the key on, and it is: `fid_eval_ch` is that state evaluated with
`Ionization transport: True`, so it is a reading of these entries and not a
substitute for one.

---

## 4. What `g` costs, and what its ceiling is for

`g = S/S'` is a property of the face the identity was measured at, and it is
measured there. The ceiling remains in the code for the one case the operator
cannot supply it: a column no face of which carried a nonzero scale, where the
measure is zero and there is no face to take a ratio at.

The two numbers differ by 1.5 (hydrogen) and 2 (helium), so the choice is
worth that factor on the tolerance and no more. What the measured form buys:

- **a state is held to the rounding of its own face**, not to the largest
  rounding the construction admits anywhere. A measure at 0.9 of the ceiling
  bound and at 1.35 of its own face's bound is refused, which is the row
  `a_measure_the_ceiling_admits_is_refused_at_its_own_face`;
- **the separation from a broken construction grows** from 1.4e10 to 2.8e10
  (1.5e-4 against the helium tolerance);
- **the report says which face's rounding was used**: the `scale:` line of the
  entry carries `g` to five decimals.

What it cannot do is change a verdict on any state measured here: all six
readings sit at 0.06 to 0.08 of the tolerance either way.

---

## 5. Tests

`src/tests/certification/run.sh`: **103 assertions, 0 failed** (84 before this
item, MEASURED by compiling the entry text of the suite against the same
objects). The nineteen new rows state the verdict of the entry on measures
chosen here rather than on one a solve happens to produce, which is what the
public `ionization_stage_sum_entry` is for:

| row | measure | verdict |
|---|---|---|
| `the_hydrogen_bound_is_the_rounding_of_its_own_sums` | the formula transcribed | 4.6629367034256575e-15, equal bit for bit |
| `the_helium_bound_carries_its_second_stage` | the formula transcribed | 1.0658141036401503e-14, equal bit for bit |
| `an_exact_identity_is_evaluated` | 0 | `evaluated` |
| `and_it_is_gated_against_the_bound` | 0 | `tol` is the bound, not 0 |
| `an_exact_identity_is_within_the_bound` | 0 | within |
| `half_the_bound_is_within_it` | 0.5 x bound | within |
| `one_and_a_half_times_the_bound_is_above_it` | 1.5 x bound | ABOVE |
| `the_helium_entry_is_gated_against_its_own_bound` | | the helium bound, not the hydrogen one |
| `a_broken_construction_is_above_the_bound` | 1.5e-4 | ABOVE |
| `the_entry_text_accepts_the_broken_construction` | 1.5e-4 | within, under the transcribed entry text |
| `a_non_finite_measure_is_not_within_the_bound` | NaN | ABOVE |
| `an_unmeasured_identity_is_unavailable` | no measurement | `unavailable` |
| `no_transported_stage_leaves_no_hydrogen_entry` | key off | no entry |
| `nor_a_helium_one` | key off | no entry |
| `the_tolerance_is_the_bound_at_the_measured_g` | `g = 1` | 3.1086244689504383e-15, bit for bit |
| `and_not_the_bound_at_the_ceiling` | `g = 1` | strictly below the ceiling bound |
| `a_measure_the_ceiling_admits_is_refused_at_its_own_face` | 0.9 x ceiling bound, `g = 1` | ABOVE |
| `while_the_same_measure_without_a_g_is_admitted` | the same measure, no `g` | within |
| `an_unmeasured_g_falls_back_on_the_ceiling` | `g < 0` | the ceiling bound |

1.5e-4 is the size the constructions that break the identity on purpose
reach (`src/tests/ionization_stage_flux/`: 1.48e-4 to 3.62e-4, READ from D3a
section 6).

**RED, before.** Two readings of it, because the routine the new rows call
does not exist in the entry text and a suite that does not compile states
nothing:

- the transcribed entry text is a row of the suite itself,
  `the_entry_text_accepts_the_broken_construction`: with `tol = 0` and
  `within_tol = .true.` written unconditionally, a measure of 1.5e-4
  certifies;
- on all three production states the control binary prints
  `tol= 0.0E+00 within` for both entries whatever the measure is, and the
  measured binary prints `tol= 4.7E-15` and `tol= 1.1E-14` (section 3).

**GREEN, after**: the table above, and the same three states within their
tolerances.

**The suites that link either changed module**, all with
`EXHALE_OBJDIR=build_D3b`, single threaded:

| suite | passed | failed |
|---|---|---|
| `certification` | 103 | 0 |
| `stage_row_balance` | 10 | 0 |
| `ionization_stage_flux` | 45 | 0 |
| `element_operator` | 40 | 0 |
| `adv_static_limit` | 53 | 0 |
| `krylov_and_dogleg` | 334 | 0 |
| `steady_species_rows` | 197 | 0 |
| `carrier_helium_inventory` | 15 | 0 |
| `charge_exchange_rows` | 250 | 0 |

The last two link the carrier module the `g` measurement was added to.

`the_stage_sum_identity_is_within_its_rounding_bound_hydrogen` and
`..._helium` of `stage_row_balance` are the executed validation of the bound
itself: largest ratio to the face's own bound 0.0703 and 0.0580, `g` 1.00 to
1.50 (MEASURED here, reproducing D3a section 6).

---

## 6. Key-off byte identity

The gate exists only where `Ionization transport: True`, so a run with the
key off carries no such entry at all and nothing in it can move. MEASURED,
control against measured build, single threaded, on scratch copies of the
regression cases:

| case | steps | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|
| `wasp_full` | 400 | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_carrier` | 12000 (its pinned `maxsteps`) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |

**Three cases of the matrix DO transport the stages**, and the gate is live in
them. `hp_zero_seed`, `hp_trace_seed` and `hp_front` carry
`Ionization transport: True`; each was run at its pinned 100 steps with both
binaries:

| case | `... sum H` | `g` | of tol | `... sum He` | `g` | of tol | `Hydro_ioniz.txt`, `Ion_species.txt` |
|---|---|---|---|---|---|---|---|
| `hp_zero_seed` | 2.210e-16 | 1.00000 | 0.071 | 3.682e-16 | 1.00000 | 0.069 | data rows BYTE-IDENTICAL |
| `hp_trace_seed` | 2.216e-16 | 1.00000 | 0.071 | 3.496e-16 | 1.00000 | 0.066 | data rows BYTE-IDENTICAL |
| `hp_front` | 2.181e-16 | 1.00000 | 0.070 | 3.857e-16 | 1.00000 | 0.072 | data rows BYTE-IDENTICAL |

All six readings are `within`, so no case of the matrix changes its verdict.

No golden was refreshed. NOTED, and not this item's: the three cases above
differ from `backup/regression/golden/` in both compared files with the
CONTROL binary as well, so the live tree has moved past the goldens through
another item's change this increment, not through this one.

---

## 7. What changed, file by file

- `src/modules/functions/ionization_stage_transport.f90`:
  `ionization_stage_sum_rounding_bound(nk [, g])` and
  `ionization_stage_sum_scale_ratio_limit(nk)`, both public, with the
  derivation and the executed validation in their comments; and an optional
  `Ek_out` on `ionization_stage_face_flux`, which returns the stage eddy
  terms the bound is written against so that they keep one spelling. Nothing
  else in the module was touched and no existing argument moved.
- `src/modules/lower_atmosphere/diffusive_photochemistry.f90`:
  `ionization_stage_nucleus_sum` forms `S` and `S'` from those eddy terms and
  returns `S/S'` at the face that carries the measure through an optional
  `g_out`; `stage_sum_g(2)` records it beside `stage_sum_max`; and
  `ionization_stage_sum_measure` returns it through an optional `g`. `dmax`
  and `jworst` are formed exactly as they were.
- `src/modules/time_step/certification.f90`: the entry construction moved
  out of `certify_state` into the public `ionization_stage_sum_entry`, which
  takes the measured `g`, sets the tolerance from the bound at it, sets
  `finite` from the measure, takes the verdict from both and names `g` in the
  entry's `scale:` line. The loop over the two elements is otherwise as it
  was.
- `src/tests/certification/certification_contexts.f90`: the two new cases of
  section 5, and the local `entry_index` replaced by the module's own
  accessor (section 8).
- `docs/certification_tolerance_anchoring_20260910.md`: anchor (8).
- `docs/lhs1140b_stationary_D3b_20260918.md`: this file.

---

## 8. Noticed outside this item's scope

1. **Fixed here, in this item's own test file.** The local `entry_index` of
   `src/tests/certification/certification_contexts.f90` returned 1, the
   index of the first entry, when the inventory held no entry of that name,
   so an assertion about a row that had disappeared would have read the
   hydrodynamic mass row's status instead of failing. The module already
   exports `certification_entry_index`, which returns 0 for an absent entry;
   the test now uses it through an `entry_status` that returns -1 when the
   entry is absent, and the five call sites go through that. No assertion
   changed its verdict: the entry text of the suite passes 84 rows against
   these objects, and those 84 are 84 of the 103.

2. **Reported by this item and corrected in the tree since**: the header of
   `stage_sum_rounding_bound` in
   `src/tests/stage_row_balance/stage_row_balance_tests.f90` gave the bound
   at `g = 1` as 1.55e-15 and 2.66e-15, which is `(5 nk + 2) eps` without
   the factor 2 the same comment's own derivation carries. It now reads
   3.109e-15 and 5.329e-15, which is what the routine's own DIAGNOSTIC line
   prints (`2.0d0*7.0d0*eps_d` and `2.0d0*12.0d0*eps_d`) and what D3a
   section 6 and anchor (8) state. The executed rows always used the correct
   expression.

3. **Reported.** The `reason` text of an entry is printed only when the
   entry is `unavailable` or when it refuses, so the new
   `gating at the floating-point bound of an algebraic identity (anchor 8)`
   is carried by the report and appears in a log only if one of these rows
   ever refuses. The same was true of the text it replaces. What a reader
   sees change on an entry that passes is the `tol=` field and the `scale:`
   line, which now names the `g` the tolerance was taken at.

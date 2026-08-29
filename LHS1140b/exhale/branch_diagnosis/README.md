# Are the two converged states of the K_zz scan physics?

Measurements made 2026-08-28 on the two solutions
`kzz_profile_scan/kzz1e7/heh3p4/k02` (cool: `log10 Mdot` 7.54, He 10830
red-pair EW 0.4287 %A) and `kzz_profile_scan/kzz1e7/heh3p399/k03` (hot: 7.87,
1.3248 %A), which the scan reports at the same K_zz and the same reservoir
He/H to 0.03 %.

## Verdict

**The hot state is a numerical artifact, not a second physical state.** It
carries a contiguous band of cells, 1.0333-1.6835 R_p, in which the total
number density of carbon, nitrogen and oxygen -- every ionization stage
summed -- is exactly zero, while the same elements sit at the handoff
reservoir ratio immediately below and immediately above it. There is no sink
for elemental carbon in these equations, so a shell emptied of carbon between
two shells at C/H = 2.639e-4 is not a state of the atmosphere.

That band is also what makes the state hot. On the cool solution the C I
metal lines carry **95.6 %** of the volume-integrated cooling over
1.0-1.2 R_p and 91.2 % over 1.2-2.0 R_p. Deleting the carbon removes that
cooling: the integrated cooling of 1.0-1.2 R_p falls by a factor 23
(1.94e-8 -> 8.29e-10 erg/s) at 0.91x the heating, and `T max` rises from
3899 K to 5879 K.

**Causality is measured, not inferred.** Refilling the deleted C, N and O at
the handoff's own reservoir ratio and re-solving the wind with the *same* hot
handoff collapses the hot state onto the cool branch:

| | `log10 Mdot` | `T max` | red-pair EW |
|---|---|---|---|
| cool solution | 7.540 | 3899 K at 1.340 R_p | 0.4287 %A |
| hot solution | 7.870 | 5879 K at 1.377 R_p | 1.3248 %A |
| **hot, C/N/O refilled, re-solved** | **7.530** | **3890 K at 1.334 R_p** | **0.4065 %A** |

## What this changes in `kzz_profile_scan/results.txt`

- Result 1, "the wind has two converged states at every K_zz that was
  reached": the second state is the defect. There is one branch, the cool
  one.
- The reservoir He/H = 11.73 recorded at K_zz = 1e9 is on the surviving
  branch, and both arms that bracket it are free of the band, so as far as
  this defect goes it stands as *the* crossing rather than "the crossing on
  the cool branch". Every other caveat already attached to it -- above all
  that the wind is subsonic through the whole domain and the continuum result
  is unvalidated -- is untouched.
- The hot-branch column of the answer table (1.64 at 1e9, 2.26 at 1e8, 2.889
  at 1e7) is not a result.
- Result 2, "at K_zz = 1e7 the cool branch ends at He/H = 3.8 ... the cool
  state ceases to exist, it does not merely fail to converge", is **refuted**.
  The 4.2 continuation acquired the band. Refilling it and re-solving
  continues the cool branch smoothly (`cool_branch_continuation.txt`):

  | He/H | 3.0 | 3.4 | 3.8 | 4.2 (this work) |
  |---|---|---|---|---|
  | `log10 Mdot` | 7.560 | 7.540 | 7.520 | 7.520 |
  | red-pair EW [%A] | 0.3931 | 0.4287 | 0.4661 | 0.5287 |

  So there is no fold and no turning point at 3.8 -- the branch was
  interrupted, not terminated. Whether the cool branch reaches the measured
  EW of 1.108 %A at K_zz = 1e7 is now an open question that the scan did not
  answer; it would need the ladder continued past 4.2 with the band watched
  for.
- The four `kzz3p16e7` arms at `log10 Mdot` 7.76-7.77 that fill the gap
  between the two clusters all carry the band and all have window spreads of
  0.13-0.18. They are arms in the act of acquiring the defect, not
  intermediate physical states.

### Which rows survive and which do not

Measured arm by arm in `affected_results_rows.txt` and
`branch_vs_dropout.txt`:

- **Clean, and the crossings rest on them**: every `flux_closure/` arm
  (`heh3`, `heh5`, `heh8`, `heh9p7`, `heh10p3`, `heh10p7`, `heh11p1`,
  `heh12`, `heh12_ctl`, `hi`, `lo`, `ref`), so the K_zz = 1e9 crossing
  11.717-11.731 is unaffected; and `kzz1e8/heh11p1`, `heh13p620`,
  `heh14p643`, `heh15p090`, `heh18p0`, so the K_zz = 1e8 crossing
  14.798-14.808 is unaffected. Also `kzz3p16e7/heh11p1`, `kzz1e8/heh3`,
  `kzz1e7/heh3`, `heh3p4`, `heh3p8`, `kzz1e6/heh3`, `kzz1e3/heh11p1`.
- **Contaminated**: the `heh11p1` EW ladder loses its 1.78e7 rung -- that arm
  carries the band (`T max` 8266 K against 5464 K at 3.16e7 and 5644 K at
  1e8), so the ladder 1.0853 -> 1.0251 -> 0.9704 -> 0.9788 should stop at
  3.16e7. The `T max` comparison quoted for K_zz = 1e7 at He/H = 11.1
  compares one band arm (1e7, 7706 K measured here) with another (1.78e7,
  8266 K), not a hot state with a cool one. The four `kzz3p16e7` arms above
  13.7 that "leave the cool branch" all carry the band.

### How much margin the surviving arms have

Three arms that carry no band were re-solved from their own converged state
with the composition diagnostic on (`clamp_margin.txt`). A settled arm barely
moves and stays well below the clamp:

| arm | diffusion steps | max X | 1 - max X | `X = 1` hits |
|---|---|---|---|---|
| K_zz 1e9, He/H 11.1 (the 11.73 crossing) | 122 | 0.9779 | 2.2e-2 | 0 |
| K_zz 1e8, He/H 14.643 | 273 | 0.9832 | 1.7e-2 | 0 |
| K_zz 1e7, He/H 3.8 | 173 | 0.9381 | 6.2e-2 | 0 |

against 7916 steps and 142 hits in the solve that acquired the band. The
defect appears in the arms that make a long excursion, not in the arms that
sit still -- but the margin is only 2 % in helium mass fraction at
He/H = 11.1, so it is not a large one.

### Blast radius

Every `output*/Ion_species.txt` in the EXHALE tree was scanned
(`dropout_blast_radius.txt`): 550 files, 130 carry a partial zero band. 128
are inside `LHS1140b/exhale/kzz_profile_scan`; the other two are
`WASP-121b/output_nan_20260810` and `output_pre_newtonfix_20260810`, one cell
each at 1.6098 R_p, both already on record as failed runs.

No regression golden and no `backup/regression/` case output carries one
(`regression_golden_dropout.txt`), so a fix to the operator would not have to
move a golden on that account.

## The measurements

| file | what it holds |
|---|---|
| `compare_branches.py`, `T1_profile_diff_kzz1e7_heh3p4.txt`, `T1_summary_kzz1e7_heh3p4.txt` | T1: the two states overlaid, relative difference against radius |
| `T1_divergence_onset.txt` | T1: the smallest radius at which each quantity separates |
| `T1_collisional_validity.txt` | `src/utils/collisional_validity.py` on both states |
| `T2_fixed_composition.py`, `T2_fixedBC.log`, `T2_fixedBC/`, `T2_profile_identity.txt` | T2: the 2x2 of handoff x seed at a fixed boundary condition |
| `T3_time_marching.py`, `T3_march/`, `T3_march_30k/`, `T3_30k_*.log`, `T3_march_endpoints.txt`, `T3_march_1000_endpoints.txt` | T3: the RK march started from each state |
| `T4_channel_budget.py`, `T4_channel_budget.txt` | T4: heating and cooling by channel on both states |
| `T4_carbon_state.txt`, `T4_metal_dropout.txt` | T4: the carbon and oxygen partition, and the extent of the zero band |
| `T4_dropout_census.txt` | every closure iteration of the scan, with or without the band |
| `branch_vs_dropout.py`, `branch_vs_dropout.txt` | the band against `log10 Mdot` and the equivalent width, arm by arm |
| `metal_dropout_reproduce.py`, `metal_dropout_reproduce.txt`, `metal_dropout_runs/` | the band created again from a clean seed with the current binary |
| `metal_restore_causal_test.py`, `metal_restore_causal_test.txt`, `metal_restore_EW.txt`, `metal_restore/` | the band refilled and the wind re-solved |
| `cool_branch_continuation.py`, `cool_branch_continuation.txt`, `cool_continuation/` | the cool branch carried past He/H = 3.8 |
| `T_restored_vs_cool.txt`, `T1_profile_diff_restored_vs_cool.txt` | the refilled state against the cool state, profile by profile |
| `clamp_margin.py`, `clamp_margin.txt`, `clamp_margin/` | how near the surviving arms come to the `Xhe = 1` clamp |
| `affected_results_rows.txt` | `T max` and the band, arm by arm, for the rows `results.txt` quotes |
| `dropout_blast_radius.txt` | every `output*/Ion_species.txt` in the tree scanned for the band |
| `regression_golden_dropout.txt` | `backup/regression/` goldens and case outputs checked |
| `hd209_precedent_dropout.txt` | the HD 209458 b baseline runs checked for the same band |

Nothing in `src/` was changed and no existing output was overwritten. Every
run made here is inside this directory.

## T1 -- the two states differ from the base upward, not beyond the exobase

The grids are identical cell for cell. The relative difference exceeds 2 % by
`r = 1.0004 R_p` in T, v, rho, x(HII), x(HeII) and n(2^3S) alike, and 30 % by
1.0013 in T. At 1.5 R_p, T is 3731 K against 5365 K (`_adv`); at 2 R_p, rho
differs by a factor 2.0 and n(2^3S) by 1.8; at 3 R_p, n(2^3S) by 3.0.

So this is **not** the pathology of `docs/Update_EXHALE.md` section 83. There
the two solutions differed only beyond the exobase and agreed below it; here
they differ everywhere, and most of the He 10830 difference is made in the
collisional part of the atmosphere (`Kn_bulk` about 1e-3 at `T max`).

`collisional_validity.py` on the two:

| | cool | hot |
|---|---|---|
| peak heating | 1.108 R_p | 1.066 R_p |
| `T max` | 3876 K at 1.346 R_p | 5391 K at 1.442 R_p |
| critical point | none in domain, max Mach 0.562 | none in domain, max Mach 0.577 |
| exobase (`Kn_bulk = 1`) | 20.823 R_p | 27.321 R_p |
| `Kn_bulk = 0.1` at | 5.743 R_p | 6.860 R_p |
| `log10 Mdot` | 7.515 | 7.847 |

(The `T max` here is from the `_adv` profiles the tool reads; the equilibrium
profiles give 3899 K and 5879 K.)

Both are subsonic through the whole domain, so both carry the standing caveat
of `docs/collisional_validity.md`: the mass flux is set at the outer boundary
and the continuum result is unvalidated. That caveat is untouched by anything
here and applies to the surviving branch too.

## T2 -- the handoff is not what distinguishes them

The two arms' handoffs state the same base: T0 to 4.4e-7, R0 to 4.6e-7, He/H
to 1.4e-4, C/H to 2.6e-3 relative. The only large difference between the two
files is the trial escape flux in the header (`F_He` 1.456e7 against
5.095e7 g/s), and EXHALE does not impose it -- `input_read.f90`
(`apply_lower_atmosphere_profile`) takes T0, R0, p_base, q_H2, He/H, the
metal ratios and `K_zz(p)` from the profile and nothing else.

The 2x2, each cell one `solve_escape_wind` call with no composition pass:

| | seed = cool | seed = hot |
|---|---|---|
| handoff = cool arm | info 0, `log10 Mdot` 7.540 | info 0, **7.870** |
| handoff = hot arm | info 0, **7.540** | info 0, 7.870 |

The result depends on the seed alone. The composition Picard loop is not what
holds the two states apart -- the seed's own C/N/O column is.

## T3 -- the march does not settle the question, and here is why

Each state was handed to the RK integrator with the steady solver off, the
same handoff, PLM (the scheme both were converged under) and `du_th` set below
anything the run can reach, so the stop is the step cap.

| state | `T max` | `log10 Mdot` | red-pair EW | zero-C cells |
|---|---|---|---|---|
| cool reference | 3899 K | 7.54 | 0.4287 | 0 |
| hot reference | 5879 K | 7.87 | 1.3248 | 172 |
| 1 000 steps from cool | 3885 K | 7.54 | 0.2159 | 0 |
| 1 000 steps from hot | 5874 K | 7.87 | 0.7924 | 172 |
| 30 000 steps from cool | 3845 K | 7.55 | 0.3976 | 0 |
| 30 000 steps from hot | 5782 K | 7.88 | 1.2648 | 172 |

**Neither JFNK root is a fixed point of the marching discretization.** Both
leave immediately -- `du` grows from about 0.005 at step 2 -- and the flux
spread rises to 16 (from cool) and 10 (from hot) around step 5000 before
falling back to 2.0 and 1.7 by step 30 000, still far from flat. The steady
residual the marching path reports on the loaded state, 0.5 in its own
`||R||(ref)`, is the same quantity the `Resid tol` gate uses, against the
1.3e-4 and 1.1e-4 the JFNK reported for the same states -- so the two solvers
do not agree on what a steady state is here. That is a separate finding and it
is not resolved here.

So T3 cannot certify or reject either state. What it does show is that the
separation is carried by the composition and not by the flow: through a
30 000-step excursion and back, each march stays on the branch it started on,
because the metal column cannot change. That is the same conclusion T2 reaches
from the other side.

## T4 -- the mechanism

Cooling, share of the volume-integrated total:

| band | channel | cool | hot |
|---|---|---|---|
| 1.0-1.2 R_p | C I | 95.6 % | 0.9 % |
| | recombination | 2.7 % | 63.6 % |
| | bremsstrahlung | 1.1 % | 26.7 % |
| 1.2-2.0 R_p | C I | 91.2 % | 52.6 % |
| | He I collisional excitation | 1.5 % | 16.8 % |

Integrated cooling over 1.0-1.2 R_p: 1.94e-8 (cool) against 8.29e-10 (hot)
erg/s, a factor 23. Integrated heating over the same band: 3.64e-8 against
3.30e-8, a factor 0.91, and its channel shares barely move (He I 80.0 % ->
80.8 %). The heating is the same; the cooling is gone.

The hypothesis this was set up to test -- that higher ionization destroys the
neutral carbon and so removes its cooling, a positive feedback -- is **not**
what happens. Carbon in the hot state is not C II or C III over that band; it
is absent:

| r [R_p] | C I cool | C II cool | C/H cool | C I hot | C II hot | C/H hot |
|---|---|---|---|---|---|---|
| 1.010 | 1.428e8 | 1.266e6 | 2.632e-4 | 1.324e8 | 1.154e6 | 2.639e-4 |
| 1.050 | 8.749e6 | 1.469e5 | 2.632e-4 | 0 | 0 | **0** |
| 1.499 | 7.216e4 | 5.435e3 | 2.632e-4 | 0 | 0 | **0** |
| 1.800 | 2.422e4 | 3.612e3 | 2.632e-4 | 2.633e4 | 2.571e3 | 2.639e-4 |

Cells 117-288 of 504, r = 1.0333-1.6835 R_p, hold exactly zero C, N and O in
every stage, in the equilibrium profiles and in the `_adv` profiles alike.
Hydrogen and helium are present there (nH = 2.1e9 cm^-3, He/H = 1.25 at
1.154 R_p). Mg is zero on both states everywhere, which is correct: the
handoff carries `--atoms H,He,N,O,C` and no magnesium.

The output is written list-directed, so `0.0000000000000000` is an exact
zero and not an underflowed small number.

## Where the band comes from

`T4_dropout_census.txt` lists every closure iteration of the scan. The band
is present in every arm the scan calls hot and absent from every arm it calls
cool. It is **created inside a wind solve**, not carried in from the
chemistry: `metal_dropout_reproduce.txt` re-runs the step
`kzz1e7/heh3p8` (clean) -> `kzz1e7/heh4p2` with the current binary and
reproduces it exactly -- 0 zero cells in, 157 out, r = 1.0429-1.6716 R_p, the
same band the recorded arm has, with `info = 0`.

Once created it is permanent. It survives every later JFNK solve, every
closure iteration, the restart path of `load_IC` (whose rebuild-from-abundance
branch fires only when an element is absent from the *whole* column, `tmp <= 0`
summed over all cells, so a partial band passes and is merely rescaled), and
the RK march.

**The code path that does exactly this** is in
`src/modules/functions/binary_element_diffusion.f90`. `element_diffusion_step`
clamps the helium mass fraction at line 439,

```fortran
where (Xhe .gt. 1.0d0) Xhe = 1.0d0
```

and `project_elements` (line 1477) then computes, cell by cell,

```fortran
nucH_new  = (1.0d0 - Xhe(j))*msum(j)/m_1
if (nucH_old(j) .gt. 1.0d-30) then
   rH = nucH_new/nucH_old(j)
else
   rH = 0.0d0
endif
...
do im = 1, n_mion                       ! metals slaved to hydrogen
   f_sp(j,mion_fsp(im)) = f_sp(j,mion_fsp(im))*rH
enddo
...
if (nucH_new  .gt. gotH)                                         &
   f_sp(j,isp_HI)  = f_sp(j,isp_HI)  + (nucH_new  - gotH)
```

At the clamp `Xhe = 1` gives `nucH_new = 0`, hence `rH = 0`: every hydrogen
species and every metal ion in that cell is multiplied by zero. Hydrogen is
then restored -- by the shortfall deposit, either at that step or at the next
one, where `nucH_old <= 1e-30` sends `rH` to zero again while `nucH_new > 0`
puts the hydrogen back as neutral H. **The metals have no such deposit.**
Their zero is final: every later call multiplies it by something.

Measured in the reproduction run: the helium mass fraction reaches exactly
1.0 in **142 of 7916** diffusion steps (`X min` never reaches 0), and that
same solve turns 0 zero-carbon cells into 157. The helium mass fraction is
already 0.93-0.98 at these reservoir ratios (He/H = 3.4 gives X = 0.931,
He/H = 11.1 gives 0.978), so the clamp is close at hand and a transient
excursion reaches it.

That causal link from `Xhe = 1` to the deleted cell was read from the source
and is consistent with every measurement above. It was **not** instrumented
cell by cell, because that would mean changing `src/`.

Once a cell's metal columns are zero, nothing short of a cold start or a
whole-column rebuild can lift them. The complete list of writers to
`f_sp(:,mion_fsp(...))` is `ionization_equilibrium.f90:971` (repartition
within an element, `nm_tot` preserved exactly),
`binary_element_diffusion.f90:554/558/560` (the trace-metal branch, gated on
`he_metal_diffusion`, off in these runs), `binary_element_diffusion.f90:1525`
(the `rH` multiplication above, ungated), `set_IC.f90:191` (cold start) and
`load_IC.f90:330/390/447` (restart). Every one of them either multiplies or
rebuilds the whole column.

**No convergence test in this pipeline can see it.** The elemental
conservation check of the handoff (`lower_profile_schema.py:509`) is stated on
the *chemistry* profile at its deepest level, not on the wind. The closure
driver tests the radial spread of the elemental *fluxes* over a window at
2 R_p and above, outside the band. The steady residual is a residual of the
hydro and energy equations, whose solution with the metals removed is a
perfectly good solution of the equations as posed -- which is why both states
report `info = 0` and `||R||` of a few times 1e-4.

## Precedent -- the HD 209458 b `K_zz = 0` case is a different failure

`docs/Update_EXHALE.md` section 70 and `docs/binary_diffusion_design.md`
section 9.2 record two steady HD 209458 b states at `K_zz = 0` a factor 2.8
apart in `Mdot`, both `info = 0` with `||R||` inside 1e-3, on compositions 1 %
apart, and leave it open as "a property of the steady solver at this
configuration".

Every `backup/phase_d_baseline/` run was checked here
(`hd209_precedent_dropout.txt`): **none of them has a partial zero band.** Si
and K are absent from all of them, which is correct -- they are not in that
`metals.inp`. So the HD 209458 b bistability is not this defect and remains
open. Caveat: only the final state of each run is on disk; the odd and even
passes of the limit cycle in `new_kzz0_undamped` were not saved with species
densities, so the alternating states themselves were not checked.

`docs/Update_EXHALE.md` section 83 is a third, separate phenomenon: there the
two solutions agreed below the exobase and differed only above it. The three
should not be merged.

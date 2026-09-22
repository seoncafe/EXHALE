# LHS 1140 b, M2: the atomic-to-molecular conversion audited on its own checkpoint

Phase 8 of `docs/PLAN_20260920_rev9.md` section 15, the audit of section 11.3,
carried out on 2026-09-21.

**Mode R throughout** (rev9 section 3): every number below is an evaluation of
the CURRENT operator on an immutable stored state, or a conversion run that
reads a stored atomic state and writes only into a scratch directory. No state
was published, no golden was refreshed, nothing was written into
`LHS1140b/models/`, and the conversion was not changed. `load_IC` rebuilds the
boundary from the physical column and this run's reservoir, so none of this is
a reproduction of a historical certificate.

Every number is labeled MEASURED (this document ran it), READ (from a named
file) or INSPECTED (source, with `file:line`). One number is labeled DERIVED and
says from what.

---

## 1. Verdict, first

**The conversion is not defective, and the state it produces is not the reason
the rows refuse.** The two possibilities the plan names do not exhaust this
case, so the answer is stated in three parts:

1. **No invariant of the conversion failed.** Every invariant it promises holds
   at the arithmetic floor on the state it actually runs on today: the H and He
   nucleus budgets, the mass budget, composition admissibility, the chosen
   thermodynamic invariant, and the equation-of-state round trip
   (section 4). The charge invariant, which the conversion computes and never
   compares, HOLDS: measured, the ion charge density per unit mass is bitwise
   unchanged in 503 of 504 rows and differs in the remaining one by 1.5 times
   the rounding floor of the census sum itself (section 5).

2. **The refusing rows are an off-equilibrium state, so in the plan's terms
   this is an initialization problem and not a conversion defect.** But the
   decisive measurement is stronger than that and it disqualifies the seed
   partition as the cause: **with the partition set to zero, that is with
   nothing converted at all and the interior identical to the certified atomic
   wind to 1.3e-15, the molecular configuration still puts the base mass row at
   1.001 of its scale and the base face mass flux at -5997 times the wind flux**
   (section 8). The same column, in the ATOMIC configuration, certifies, with
   the base face carrying +1.000 of the wind flux (MEASURED). What separates the
   two is the lower boundary the molecular run builds, not the H2 the conversion
   writes.

3. **What the measurement does not decide**, and what would: whether that base
   boundary state is itself wrong or merely far from the fixed point of the
   molecular equations. Nothing here establishes that either way. The
   measurement that would decide it is named in section 9.

The work this points at is therefore neither repairing the conversion nor
choosing another seed partition. It is the molecular base boundary of this
column: the reservoir the molecular base particle count implies stands 5.07 per
cent denser than the atomic one at the same level, and the characteristic face
condition turns that into a face velocity 12.3 times larger and of the opposite
sign to the wind (MEASURED, section 8).

---

## 2. Identity

### 2.1 The checkpoint

| | |
|---|---|
| case | `molecular_photochem_gj1132_kzzprofile/HeH9` |
| generation | `g0005_20260920T053700Z_f59de41d` (`latest_complete`; the index names no certified generation) |
| record | `docs/lhs1140b_identity_records_20260921.md` section 3.5 |
| `Hydro_ioniz.txt` md5 | `348a6ab88a84b1b967d7582f14f9a723` (MEASURED before the first run and again after the last, unchanged; equal to the manifest) |
| `Ion_species.txt` md5 | `e39c8734708d77894f578fe1df25283c` (MEASURED, same) |
| boundary model of the state | `characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3` (READ, state header) |
| options | `mol=T molbase=T carrier=T carrier_newton=F iontrans=F oxychem=F he_diff=T sec_ion=T visc=F cond=F wellbal=T` (READ, state header) |
| historical binary | `EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c` (READ, manifest) |

Rows of its certificate (READ, `states/<gen>/certification.txt`): mass
1.041E+00 against 3.4E-10 at cell 2, momentum 1.450E-03 at cell 1, energy
1.021E+00 at cell 1, carrier H2 gated 1.000E+00 at cell 432, elemental He/H
gated 1.494E-04 at cell 217.

**The two decade counts, kept apart.** The CURRENT generation: MEASURED with
the unrounded values of this document, the mass row measure is 1.040949 and its
cell-2 tolerance is 10 times the cell's own rounding floor, 10 x
3.3560314651785441E-11 = 3.3560E-10, so the distance is 3.1017E+09, **9.4916
decades**, which is the `distance 3.102E+09` the certificate prints. The figure
**10.01 decades** (2.254E-01 against 2.2E-11) belongs to the 2026-09-18 memo
`not_solved.md` under binary `3146d11b`, is evidence about that occasion only,
and is not mixed with the above anywhere in this document.

### 2.2 The operator

| | |
|---|---|
| binary | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `a55d90086af7a18699ed86d2d39b31b2` (MEASURED, before and after every run, unchanged) |
| bytes, mtime | 4044168, 2026-09-21 09:16 (MEASURED) |
| compiler | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0` (MEASURED, `strings`) |
| linked libraries | `libopenblas.so.0`, `libgfortran.so.5`, `libgomp.so.1`, `libquadmath.so.0`, all from `/opt/miniconda3/lib` (MEASURED, `ldd`) |
| source manifest | none exists for this md5; working tree at `93eed86` with modifications |
| threads | `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and every run launched under `env -i` so no other variable of the calling shell reached it (MEASURED) |
| assembly selector | `EXHALE_RESID_QUAD` unset in every run, so `ROWS_PRODUCTION` in double (rev9 section 13.6) |

The current binary reproduces the checkpoint's certificate exactly to the
printed digits (MEASURED): mass 1.041E+00 at cell 2, momentum 1.450E-03 at cell
1, energy 1.021E+00 at cell 1, carrier H2 1.000E+00 at cell 432, elemental He/H
gated 1.494E-04 at cell 217, `NOT CERTIFIED: 5 entry/entries`. The conserved
columns round-trip bitwise in the radius, density and velocity columns and to
1.5e-14 in pressure.

### 2.3 The configuration

The case `input.inp` with the relative paths made absolute and
`Restart intent: stationary evaluate`, exactly the file the runner's own
post-processing pass uses (`runs/r20260920T053659Z_37074/eval/input.inp`, md5
of the scratch copy `22df71bd65a851ba5fd7a00da28219dd`, differing from the
runner's only in the profile path, which points at a byte-identical scratch
copy, md5 `8b3c9aea5834d51e836571bf082c5050`). Derived route keys, READ from
`EXHALE_resolved.out` of the case: `well_balanced T`, `carrier_transport T`,
`carrier_in_newton F`, `carrier_newton_on_stall F`, `ionization_transport F`,
`oxygen_chemistry F`. Numerical flux HLLC and reconstruction PLM as stated in
prose by `EXHALE_setup.out`, with the stationary route selecting WENO3 for the
residual, which is what the state header's `recon=WENO3` records.

### 2.4 The conversion's own input

The atomic state the conversion reads, in every run of this document:

| | |
|---|---|
| directory | `LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k04/output` |
| `Hydro_ioniz_IC.txt` md5 | `f55a75c2b1347177b217bcd8ffff1362` (MEASURED) |
| `Ion_species_IC.txt` md5 | `b6dd7275ed06caea52217baec63f8987` (MEASURED) |
| its own claim | `certified=T cert_reason=certified_in_wind mode=init` (READ, header) |
| its boundary model | `characteristic_face_ps_reservoir_C_minus_contact_upwind_v2` (READ, header), which is NOT the molecular case's `..._v3` |
| its options | `mol=F molbase=F carrier=F` (READ, header) |

---

## 3. What the conversion promises, as the source states it

INSPECTED, `src/modules/init/molecular_seed_from_atomic_state.f90`:

- Two thermodynamic invariants at `:638-648`. `invariant=p` is the default
  (`:169`, `:256-273`, key `EXHALE_MOLECULAR_SEED_INVARIANT`) and keeps rho, v
  and p, with T following the new particle count; `invariant=T` keeps T and
  reconstructs p.
- Budget checks at `:610-633` on the H nuclei, the He nuclei and the mass, with
  an `error stop` on a violation or on a negative H I or H2 fraction, against
  `molecular_seed_closure_tol = 1.0d-12` (`:205`).
- The energy is closed with the production equation of state at `:651-657`
  (`W_to_U` / `Apply_BC` / `U_to_W`), and the round-trip pressure difference is
  recorded and also held to the same tolerance.
- `transfer_h2` (`:704-748`) restores the atomic composition and then changes
  only H I and H2, both of which carry zero charge weight
  (`species_table.f90:152`, `bsp_charge`).
- A charge census is computed twice (`element_nuclei_and_charge` at `:450` and
  at `:610`, returning `nchg0` and `nchg1`) and **never compared**. That is a
  missing assertion, not a violation.
- The requested H2 fraction comes from `EXHALE_MOLECULAR_SEED_X2` (`:234-252`):
  absent means `handoff`, `local` takes the smaller of the thermochemical fit
  and `2 n_root / n_H` from `carrier_h2_chemical_root`, and a real in [0,1]
  states the fraction.

Neither invariant claims to preserve the source model's chemical or energy
equilibrium.

---

## 4. Measurement 1: the promised invariants, including the boundary

### 4.1 The conversion as it actually runs, today

Four conversions were run from the atomic state of section 2.4 under the binary
of section 2.2, each in its own scratch directory, each writing only
`output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt` and stopping.
MEASURED, from the `(molecular_seed)` report of each run:

| | `local`, `p` | `handoff`, `p` | `0`, `p` | `local`, `T` |
|---|---|---|---|---|
| `q_H2_base` | 5.1230983685145566E-02 | 5.2506039495835946E-02 | 0 | 5.1230983685145566E-02 |
| ceiling `0.5/(0.5+He/H)` | 5.2574751754452698E-02 | same | same | same |
| x2 requested (largest) | 9.7568641868732442E-01 | 9.9875825506558791E-01 | 0 | 9.7568641868732442E-01 |
| largest x2 transferred | 9.7568641868732442E-01 | 9.9875825506558802E-01 | 0 | 9.7568641868732442E-01 |
| cells capped by their own neutral H | 0 of 504 | **385 of 504** | 0 of 504 | 0 of 504 |
| cells clipped at the element-ratio ceiling | 191 | (not applicable) | (not applicable) | 191 |
| max rel. change of the H nuclei | 2.5720E-16 | 2.2036E-16 | 1.5249E-16 | 2.5720E-16 |
| max rel. change of the He nuclei | 0.0000E+00 | 0.0000E+00 | 0.0000E+00 | 0.0000E+00 |
| max rel. change of the mass | 2.3015E-16 | 4.3686E-16 | 0.0000E+00 | 2.3015E-16 |
| invariant kept | p | p | p | T |
| it moved T by (max, rel.) | 4.8732E-02 | 1.0809E-01 | 4.9884E-02 | 0 |
| it moved p by (max, rel.) | 0 | 0 | 0 | 4.8732E-02 |
| equation-of-state round trip | 3.7171E-16 | 5.2342E-16 | 2.4780E-16 | 3.1742E-16 |

Every budget stands twelve to sixteen decades below the `1.0d-12` the
conversion refuses at, and no `error stop` fired in any run, so no composition
was inadmissible and no H I or H2 fraction went negative.

**The `0`, `p` column's T displacement of 4.9884E-02 is entirely in the ghost
rows.** The maximum is taken over `1-Ng:N+Ng` (`max_rel_change`, `:757-767`),
and with x2 = 0 the interior of an atomic state is untouched while the two
lower ghosts, which `load_IC` fills with this run's MOLECULAR reservoir row and
which therefore already carry H2, have that H2 moved back into H I. The header's
statement at `:118-120` that x2 = 0 is "the conversion identity (nothing is
transferred, so every species column, the nuclei and the chosen primitive
invariants come back unchanged)" is therefore true of the interior and not of
the lower ghosts. MEASURED: over the 500 physical cells, the `0`, `p` state
agrees with the atomic source to 1.9e-16 in density, 2.2e-16 in velocity,
3.5e-16 in pressure and 1.3e-15 in temperature, so the identity claim holds
where it can.

### 4.2 The interior, mode against mode

MEASURED, each mode against the `0`, `p` conversion, so that everything except
the transfer is held fixed:

| over the 500 physical cells | `local`, `p` | `handoff`, `p` | `local`, `T` |
|---|---|---|---|
| max rel. change of rho | 0 (bitwise) | 0 (bitwise) | 0 (bitwise) |
| max rel. change of v | 0 (bitwise) | 0 (bitwise) | 0 (bitwise) |
| max rel. change of p | 4.20E-16 | 4.97E-16 | **4.87E-02** |
| max rel. change of T | **5.12E-02** | **1.21E-01** | 4.31E-16 |

This is the promise, met: `invariant=p` keeps rho, v and p to the rounding of
the equation-of-state round trip and moves T; `invariant=T` keeps T and moves p.

### 4.3 The base ghost and the base face state after `Apply_BC`

This is the part rev9 added after review1, and it has two answers that must not
be merged.

**(a) The ghost moved a great deal inside the conversion, and that movement is
the boundary replacing a placeholder.** INSPECTED, `load_IC.f90:1215-1221`: the
lower ghost rows of the restart pair are not read, and the loader sets
`rho(j) = rho(1)`, `v(j) = v(1)`, `p(j) = p(1)`, `T(j) = T(1)` for the two
lower ghosts together with the reservoir composition row. That copy of cell 1
is the ghost state the conversion's `Apply_BC` is handed. MEASURED, the ghost
`j = -1` of the written file against cell 1 of the same file:

| | cell 1 (before) | ghost `j=-1` (after `Apply_BC`) | ratio |
|---|---|---|---|
| rho [mH/cm3], `0`,`p` | 1.1779835E+14 | 1.5971311E+14 | 1.3557 |
| p [cgs], `0`,`p` | 9.1893542E-01 | 1.0970496E+00 | 1.1938 |
| v [cm/s], `0`,`p` | -1.5193948E+00 | -1.2748200E+02 | 83.90 |
| rho, `local`,`p` | 1.1779835E+14 | 1.5980833E+14 | 1.3565 |
| p, `local`,`p` | 9.1893542E-01 | 1.0970813E+00 | 1.1939 |
| v, `local`,`p` | -1.5193948E+00 | -1.2799245E+02 | 84.24 |

So yes, the ghost moved, by 36 per cent in density, 19 per cent in pressure and
a factor 84 in velocity. That is the hydrostatic isentrope continuation of the
boundary replacing the cell-1 copy, and it happens identically whether or not
anything was converted.

**(b) The conversion's own effect on the ghost is small, and rev9's mechanism
is real but three to four decades too small to be the cause.** MEASURED, the
written ghost rows of each mode against the `0`, `p` conversion:

| | ghost `j=-1` | ghost `j=0` | cells 1, 2, 3 |
|---|---|---|---|
| `local`,`p`: d rho | +5.962E-04 | +3.757E-06 | 0 (bitwise) |
| `local`,`p`: d p | **+2.890E-05** | +2.285E-06 | 0 / -1.2E-16 / -2.8E-16 |
| `local`,`p`: d T | +5.063E-02 | +5.123E-02 | +5.12E-02 |
| `local`,`p`: d v | -4.004E-03 | -4.599E-03 | 0 (bitwise) |
| `handoff`,`p`: d p | +2.961E-05 | +2.341E-06 | 0 |
| `local`,`T`: d v | **+1.849E+01** | +1.850E+01 | 0 (bitwise) |

The plan's sentence is confirmed as a mechanism: **the conversion preserves the
cell pressure exactly and still changes the ghost pressure**, by 2.9e-05
relative at `j=-1` and 2.3e-06 at `j=0`, and the ghost density by 6.0e-04. It
is not confirmed as the cause of anything: those displacements are three to
four decades below the displacement the target run's own boundary produces when
the same column is loaded (section 8), and in any case they never reach the
target run, because `load_IC` discards the file's lower ghost rows
(`load_IC.f90:1184-1213`). The channel by which a conversion can change the
base flux of the run that loads it is cell 1's composition, density and
pressure, not the ghost rows it wrote.

`invariant=T` is the exception worth naming: it moves the ghost velocity by a
factor 18.5 AND reverses its sign, from -128 cm/s to +2230 cm/s (MEASURED),
because the 4.9 per cent pressure drop it imposes on the interior changes what
the boundary continues from.

**(c) One composition, one boundary, one residual.** MEASURED in all six
evaluations of this document: the cached base face state against the boundary
of the installed composition is exactly `0.00000E+00`
(`report_base_face_state_consistency`, `EXHALE_BOUNDARY_TRACE=1`).

---

## 5. Measurement 2: charge

### 5.1 The invariant, stated

Define it as rev9 section 11.3 asks, as preservation of the **ion charge
density per unit mass** across the neutral transfer:

```
    Q(j) = sum_i q_i f_i(j) ,        f_i = n_i / (rho n0)
```

with `q_i` READ from `species_table.f90:152` (`bsp_charge`) and `:224`
(`mion_stage`), and `f_i` formed from the number densities the state file
carries (`write_output.f90:426-440` writes `n_i` in cm^-3) divided by the
density column. That is the space `transfer_h2` works in (`:704-748`), and in it
the statement is exact rather than approximate: the routine writes `f(H2)` and
`f(H I)` and nothing else, and both carry `q = 0`.

The tolerance is absolute plus relative and stays defined in a neutral cell:

```
    |dQ| <= a + b |Q| ,   with the meaningful absolute scale
    floor(j) = eps * sum_i |q_i| f_i(j) ,   eps = 2^-53 ,
```

the rounding of the census sum itself, which is nonzero wherever any ion is
present and which is the number an assertion should use.

### 5.2 The measurement

Comparing the `local`, `p` conversion with the `0`, `p` conversion isolates the
transfer exactly: both go through the same `load_IC`, the same reservoir, the
same `Apply_BC` and the same writer, and differ only in the fraction handed to
`transfer_h2`. MEASURED, over all 504 rows:

| | `local`,`p` vs `0`,`p` | `handoff`,`p` vs `0`,`p` | `local`,`T` vs `0`,`p` |
|---|---|---|---|
| rows where `dQ` is not exactly zero | **1 of 504** | 1 of 504 | 3 of 504 |
| max abs. `dQ` [charges per H mass] | 5.293956E-23 | 5.293956E-23 | 1.110223E-16 |
| max rel. `dQ/Q` | 1.692448E-16 | 1.692448E-16 | 3.388053E-16 |
| **max `dQ` in units of the floor** | **1.524** | **1.524** | **3.052** |
| where | row 1, r = 1.000000 (ghost `j=0`) | same | row 502, r = 29.5156 |
| H nuclei per unit mass, max rel. change | 4.46E-16 | 2.55E-16 | 4.46E-16 |
| He nuclei per unit mass, max rel. change | 2.27E-16 | 1.13E-16 | 2.27E-16 |

**The charge invariant holds at the arithmetic floor. The floor is
`eps * sum_i |q_i| f_i`, and the largest departure anywhere is 1.5 of them for
the default invariant and 3.1 for `invariant=T`.** An assertion added to the
conversion should assert exactly this, at a few floors, and it will pass on this
state. The residue that exists is not the transfer at all: it is the file write
and read, since the comparison is made through two written files.

For reference, the same quantity measured against the atomic source FILE, which
also carries `load_IC`'s own reprojection: over the 500 physical cells,
7.2240E-16 relative for every mode (MEASURED), and over all 504 rows 3.2227E-01,
the 32 per cent being the two lower ghosts, which `load_IC` replaces with the
reservoir composition row. That is the separation rev9 asks for: the boundary
answers to its own reservoir and charge closure, and preserving each old ghost
charge across `Apply_BC` is not the invariant of an interior neutral transfer.

### 5.3 Whether this is a test or an identity

**It is a test of the charge invariant and it would be an identity for
quasineutrality.** The charge density above was formed here, from the species
number densities the file carries and a charge table read from the source; the
code's own census (`element_nuclei_and_charge`) was not used, and the two
conversions compared are two separate processes. Quasineutrality is a different
matter: `calc_ne` builds `n_e` from the same weighted sum over the same
composition, so comparing `n_e` with `sum_i q_i n_i` would compare a number with
itself, and the state file does not carry an independently obtained electron
density. No such comparison is reported here.

---

## 6. Measurement 3: the refusing rows decomposed, in the code's own convention

### 6.1 The convention, and what is active

INSPECTED, `steady_residual.f90:294-306`:

```
    R_mass     = dF_mass     - S_mass
    R_momentum = dF_momentum - S_momentum  [ - Smom if transport is active ]
    R_energy   = dF_energy   - S_energy - (heat - cool)  [ - Sene if active ]
```

- `S_mass = 0` and `S_energy = 0` in every branch, and under the well-balanced
  option `S = 0` for all three rows (INSPECTED, `Source.f90:44-51`). The state
  carries `wellbal=T`, so the momentum branch in force is **WENO3, well
  balanced**: `dF_momentum = (A+ F+ - A- F-)/V + (q_up+ - q_dn-)/dr`, explicit
  source zero (rev9 section 13.4 table).
- **Transport is NOT active.** READ, the state header options: `visc=F cond=F`,
  so `transport_active()` is false (`viscous_conduction.f90:337-339`) and
  `Smom = Sene = 0` identically. They are named here because rev9 requires it,
  not because they contribute.
- The gravitational work of the energy row is INSIDE `dF_energy` as `dF3p/V`
  and enters with a positive sign (INSPECTED, `RK_rhs.f90:273-275`); no second
  gravitational term is subtracted anywhere below.
- **The energy zero** is `internal_energy_of_mixture = 1.5 T + x2 u_rv/T0`
  (INSPECTED, `caloric_eos.f90:486`), translational energy plus the H2
  rotational and vibrational excitation, NOT a binding energy measured from
  separated atoms. The dissociation energy lives in the chemical source
  accounting and is not added a second time here.

Code units, MEASURED by pairing the cgs columns of the written work state with
the code-unit ghost record of the same evaluation: `n0 = 4.11119947E+13 cm^-3`,
`T0 = 185.41987 K`, `p0 = n0 kB T0 = 1.05246612 erg/cm3`,
`v0 = 1.23680934E+05 cm/s`, `R0 = 1.1593574793765156E+09 cm` (READ, header),
`t_s = R0/v0 = 9.37377692E+03 s`, and the energy-row rate unit
`q0 = p0/t_s = 1.12277701E-04 erg cm^-3 s^-1`.

### 6.2 The mass row at cell 2, where it refuses

MEASURED, from `output/residual_profile.txt` and the base mass row trace
(`EXHALE_BOUNDARY_TRACE=1`), all in code units:

| cell | `A- F-/V` (lower face) | `A+ F+/V` (upper face) | `dF_mass = R_mass` | scale | `|R|/scale` |
|---|---|---|---|---|---|
| 1 | -2.241968E-02 | -4.721178E-02 | -2.479204E-02 | 4.721168E-02 | 0.5251 |
| 2 | -4.719329E-02 | +1.932519E-03 | **+4.912596E-02** | 4.719343E-02 | **1.0409** |

`S_mass` is zero, so the mass row IS the flux difference, which
`write_residual_breakdown` also prints for this cell as
`R = 4.912596E-02, flux divergence = 4.912596E-02, source = 0.000000E+00`
(MEASURED). Rebuilding the difference from the exported face fluxes, the
reconstructed face radii and the cell volume reproduces the printed row to
2.9e-06 relative, that residue being the nine printed digits of the face fluxes
and not a disagreement of the code.

**What the row says, physically.** The face mass fluxes `r^2 rho v` at the three
faces bounding these two cells are (MEASURED, code units)

```
    face 0 (the base face)  -4.334975E-06        inflow
    face 1                  -9.125140E-06        inflow, twice as large
    face 2                  +3.735211E-07        outflow
```

against a wind-window mean of `5.88765E-07` (READ, certificate). So the base
face draws mass IN at 7.36 times the wind flux, face 1 draws in at 15.5 times,
and face 2 carries out 0.63 times. The flow REVERSES between face 1 and face 2.
Cell 2 gains mass at 26 times the rate it loses it, and since the row scale is
by construction the larger of the two face terms (4.719343E-02, which is exactly
the lower-face term), `|R|/scale` slightly above one means the two terms do not
cancel at all: the row is the lower-face inflow, essentially unopposed. It is
not a cancellation that has gone wrong; there is no cancellation.

### 6.3 The energy row at cell 1, where it refuses

MEASURED, code units; `heat` and `cool` are the cgs columns of the work state
written by the same evaluation, divided by `q0`:

| cell | `heat` | `cool` | `heat - cool` | `R_energy` | `dF_energy = R + (heat-cool)` |
|---|---|---|---|---|---|
| 1 | +8.866575E-04 | +2.854822E-06 | +8.838026E-04 | **-4.351493E-02** | -4.263113E-02 |
| 2 | +9.254358E-04 | +2.387314E-06 | +9.230485E-04 | +8.304002E-02 | +8.396307E-02 |

The row scale the certificate divides by is the row's largest term, and it comes
out as `4.263114E-02` at cell 1 (from `|R| = 4.351493E-02` and the printed
`|R|/scale = 1.020731`), which is `|dF_energy|` to six digits. **So the energy
row at cell 1 is an unbalanced flux difference: the radiative term is 2.07 per
cent of it, and nothing else is there to balance it.**

Splitting the flux difference further, DERIVED (not exported by this binary):
the gravitational work is
`dF3p/V = [A+ F_mass+ (phi_i+ - phi_c) - A- F_mass- (phi_i- - phi_c)]/V` with
`phi(r) = -b0/r`, the spherical branch (INSPECTED, `grav_field.f90:16-20`, the
run has `Domain mode: Spherical`), and
`b0 = G Mp mu/(kB T0 R0) = 125.8868` (DERIVED from `Gc = 6.67430e-8`,
`MJ = 1.8982e30`, `Mp = 0.0176220 MJ`, `mu = 1.67353284e-24 g`, all INSPECTED
in `parameters.f90:842,858` and `input_read.f90:1832,2114`; it gives a base
gravity of 1661 cm/s^2, which is the planet's):

| cell | `dF3p/V` (gravitational work, inside `dF_E`) | enthalpy flux difference `dF_E - dF3p/V` | `dF3p / dF_E` | `(heat-cool)/dF_E` |
|---|---|---|---|---|
| 1 | -8.469276E-04 | -4.178420E-02 | +1.99E-02 | -2.07E-02 |
| 2 | -5.503718E-04 | +8.451344E-02 | -6.55E-03 | +1.10E-02 |

So at cell 1 the gravitational work and the net radiative source are each about
two per cent of the flux difference and happen to nearly cancel each other; the
row is the enthalpy flux difference `v(E+p)`, unbalanced. A consumer that
omitted `dF3p/V` would misstate the energy row of these cells by two per cent,
which is the size of the radiative term it is being compared against, so rev9's
warning is quantitatively justified here even though it is not what decides the
refusal.

### 6.4 The momentum row at cell 1

MEASURED: `R_momentum = -2.370512E-01` at cell 1 and `+5.128598E-02` at cell 2,
in code units, with `S_momentum = 0` (well balanced), so
`R_momentum = dF_momentum` including the WENO3 well-balanced departure term
`(q_up+ - q_dn-)/dr`. The certificate's `1.450E-03` is this divided by the row
scale `1.635292E+02`, which below the escape radius is the gravitational force
density. The momentum row is four decades quieter than the other two, relative
to its own scale, and is not where this state fails.

### 6.5 What could not be separated with this binary

`dF_energy` could not be split into the stored energy face fluxes and `dF3p`
BY THE OPERATOR: `write_residual_breakdown` prints signed terms for the single
worst cell of the whole table only (INSPECTED,
`steady_residual.f90:1046-1063`), and that cell is the mass row at cell 2. The
term export rev9 section 13.2 and phase 5 require now EXISTS in the tree as
`src/modules/time_step/conservation_budget.f90` with the key
`EXHALE_CONSERVATION_BUDGET`, and it exports `grav_work_over_volume` as a
separately named column with the sign convention stated (INSPECTED,
`conservation_budget.f90:44-56, 320-380`), but the file is untracked, was
modified at 13:14 on 2026-09-21 while this audit was running, and is NOT LINKED
into `EXHALE.x` md5 `a55d9008` (MEASURED: `strings` finds no
`EXHALE_CONSERVATION_BUDGET` in the binary; the run with the key set wrote no
file). That is phase 5's work and another worker's file, so no build was made
here. The split in section 6.3 is therefore DERIVED and should be re-measured
against that export when the phase-5 binary exists.

---

## 7. Measurement 4: the seed modes, historical and today

### 7.1 Which mode produced each historical seed of this case

Both seeds this case has ever had used **`local`** with **`invariant = p`**,
from the same atomic directory. READ:

| seed | evidence | mode | `q_H2_base` | x2 | T moved | round trip |
|---|---|---|---|---|---|---|
| 2026-09-16 (behind `output_pre_L34/`, the parent of `g0001`) | `output_pre_L34/case_files/seed.log` (md5 `872c8b8235b2592dce3143cadad47344`) | `local`, `p` | 5.1360529156327897E-02 | 9.7803306766890707E-01 | 4.8850E-02 | 3.2879E-16 |
| 2026-09-18 (behind `g0002` and `g0003`, and so behind `g0005`) | case `seed.log` (md5 `c64a0cfd75b32fa722ffa00cfef6410b`) | `local`, `p` | 5.1230983685747987E-02 | 9.7568641869823824E-01 | 4.8732E-02 | 3.7401E-16 |

Corroborated by `provenance/g0003_.../provenance_recovered.json`, which quotes
`REPRODUCE.md`: "molecular seed from
`.../atomic_photochem_gj1132_kzzprofile/HeH9/k04/output`, **x2 local**; then a
continuation from the first solve's written state at dtau0=1.0e8" (READ).
`state_index.json` gives the chain `g0001 -> ... -> g0003`, and `g0005` is an
`evaluate` of `g0003` with no step taken (READ, manifest).

Neither state file of this case carries the `# molecular-seed-from` and
`# molecular_partition` header lines the conversion writes
(`molecular_seed_from_atomic_state.f90:843-872`), because those lines are
written only by the process that performed the conversion and every stored
generation of this case was written by a later solve. The mode is therefore
recoverable only from the logs and the provenance records, which is worth
noting: had the logs been lost, the seed mode of this case would not be
recoverable from its states.

**The current binary reproduces the 2026-09-18 seed.** MEASURED, `local`, `p`
today against that log: `q_H2_base` 5.1230983685145566E-02 against
5.1230983685747987E-02 (6e-13 relative), x2 9.7568641868732442E-01 against
9.7568641869823824E-01 (1e-12), cells capped 0 of 504 in both, cells clipped at
the ceiling 191 in both, base layer cells 1 to 118 up to r = 1.0346 in both,
handoff over root 1.02 to 45.10 in both, T moved 4.8732E-02 in both.

### 7.2 What each mode gives on this state today

Each conversion was then evaluated by the same route and the same binary
(Mode R, `Restart intent: stationary evaluate`), so the four columns are the
same experiment on four seeds. MEASURED:

| | `local`,`p` | `handoff`,`p` | `0`,`p` | `local`,`T` |
|---|---|---|---|---|
| hydrodynamic mass row | 1.001E+00 at cell 1 | 1.949E+00 at cell 187 | 1.001E+00 at cell 1 | 1.454E+00 at cell 118 |
| hydrodynamic momentum row | 3.128E-02 at cell 1 | 1.464E+00 at cell 358 | 3.123E-02 at cell 1 | 5.812E-01 at cell 1 |
| hydrodynamic energy row | 1.424E+00 at cell 2 | 2.011E+00 at cell 274 | 1.058E+00 at cell 2 | 1.742E+00 at cell 12 |
| carrier balance H2 | 1.000E+00 at cell 1 | 9.989E-01 at cell 1 | UNAVAILABLE (no H2) | 1.000E+00 at cell 432 |
| elemental He/H, gated | within | 5.941E-01 at cell 358 | within | 3.155E-05 at cell 440 |
| cells without a chemical root | 0 | **98** | 0 | 0 |
| largest accepted reaction residual | 1.435E-10 | **2.784E+00** | 5.250E-11 | 1.298E-10 |
| composition mass closure | 8.719E-16 | 7.578E-16 | 8.395E-16 | 1.021E-15 |
| entries refusing | 4 | **7** | 4 | 5 |
| base face mass flux / wind mean | -6035 | -6036 | **-5997** | +126840 |
| base cell-1 mass row `|R|/scale` | 1.0005 | 1.0005 | 1.0006 | 1.0028 |

**`local` is the right historical choice and the measurement says why.**
`handoff` asks for more H2 than 385 of 504 cells have neutral hydrogen to give,
leaves 98 cells with no chemical root and a reaction residual of 2.784, and
refuses on seven entries instead of four. `invariant=T` reverses the base face
mass flux to a large OUTFLOW and is worse on every hydrodynamic row.

**And `local` buys almost nothing at the base.** The base face mass flux of the
`local` seed differs from the `0` seed, which converts nothing at all, by 0.6
per cent, and the cell-1 mass row by 0.01 per cent. Whatever the partition, the
base of this column enters the molecular run with a face inflow six thousand
times the wind flux.

**A chemistry relaxation at fixed hydrodynamics, if one is ever proposed here,
must state its thermodynamic invariant.** Section 4.2 measures why: on this
state, fixed pressure and fixed temperature differ by 4.9 per cent in the other
variable, because the heat capacity changes when 97.6 per cent of the hydrogen
nuclei become H2. Fixed conserved energy is a third experiment again and was not
run.

### 7.3 The metal-free scalar restart

Not re-measured here. `not_solved.md` records (READ) that `load_IC` refuses a
metal-free scalar state for this metal-bearing photochemical column with
"metadata field `reservoir`: this run carries a reservoir for C/H
( 2.778025321E-04) and the restart files do not". rev9 section 11.3 item 4 calls
that correct behavior and this audit found nothing against it.

---

## 8. The control that decides it: one column, two configurations

The `0`, `p` conversion writes a molecular state file whose 500 physical cells
are the certified atomic wind to 1.3e-15 and whose H2 column is zero everywhere.
Evaluating it under the MOLECULAR configuration, and evaluating the same atomic
column under its OWN configuration, is one controlled comparison of the two
equation sets and the two boundary models on one state. MEASURED:

| | atomic configuration, atomic column | molecular configuration, same column (`0`,`p`) |
|---|---|---|
| verdict | **CERTIFIED**, in the wind | NOT CERTIFIED, 4 entries |
| mass row | 3.004E-09 at cell 2, tol 8.3E-09, within | **1.001E+00 at cell 1** |
| momentum row | 3.705E-13 at cell 499, within | 3.123E-02 at cell 1 |
| energy row | 1.388E-08 at cell 20, within | 1.058E+00 at cell 2 |
| base cell-1 mass row, signed | -2.654E-12 | **+1.7057E+01** |
| base face mass flux | +5.783479E-07 | **-3.296142E-03** |
| the same, in wind-mean units | **+1.000** | **-5997** |
| interior velocity at the base face | -1.484002E+00 cm/s | -1.484002E+00 cm/s (bitwise the same) |
| reservoir density at the face | 1.397205E+14 mH/cm3 | **1.468023E+14 mH/cm3 (+5.07 %)** |
| face velocity the condition returns | -1.126628E+01 cm/s | **-1.381929E+02 cm/s (12.3 x)** |
| reservoir temperature | 182.104 K | 181.984 K |
| contact upwind `w_rev` | 0.0 (the reservoir states the entropy) | 0.0 (the same) |
| boundary model | `..._contact_upwind_v2` | `..._contact_upwind_ghost_fixed_point_seed_reservoir_row_v3` |

The interior is bitwise identical, the branch of the face condition is the same,
and the face velocity differs by a factor 12.3. What differs is the reservoir:
the molecular base particle count, `q_H2(base, photochemical) = 0.053 ->
ntot_bc = 0.950` (READ, the case's own `seed.log`), makes the gas at the same
base pressure and temperature 1/0.950 = 5.26 per cent denser, and the measured
reservoir density ratio is 1.0507.

Two further facts from the same records:

- With `0`, `p` the molecular boundary still imposes the handoff partition on
  the GHOST cells alone, so the ghost record shows `x_H2 = 0.99876` in both
  lower ghosts while cell 1 has no H2 at all (MEASURED). That composition step
  across the base face changes the base face mass flux by only 0.6 per cent
  relative to the `local` seed, in which the ghost and cell 1 are both about 98
  per cent H2, so the step is not what drives the inflow either.
- The two solves that produced `g0003`, and hence `g0005`, moved the base a long
  way and stopped: the base face mass flux went from -5997 wind means at the
  seed to -7.364 at `g0005`, and the ghost density from 1.5110E+14 to
  5.4584E+13 mH/cm3 (MEASURED). The base is still drawing in, at 7.4 times the
  wind flux, and the mass row at cell 2 is what measures that.

---

## 9. Verdict, in the plan's terms, and what would decide the rest

**Not a conversion defect.** No invariant the conversion promises failed, on the
state it actually runs on, under the current binary: nuclei, mass, admissibility
and the equation-of-state round trip all sit twelve to sixteen decades inside
the tolerance the conversion refuses at, and the charge invariant, which the
code computes and does not compare, holds at 1.5 rounding floors.

**An off-equilibrium state, so an initialization problem in the plan's
language, but not one the seed created.** The refusing rows are an unbalanced
base: the mass row at cell 2 is the lower-face inflow with nothing opposing it,
and the energy row at cell 1 is the enthalpy flux difference with a radiative
term two per cent of its size. That imbalance is present in full when NOTHING is
converted, so the seed partition, the seed mode and the thermodynamic invariant
are not its cause; they change the base face mass flux by less than one per cent
between them.

**What the measurement does not decide.** Whether the molecular base boundary
state of this column is itself wrong, or is right and merely far from the fixed
point that the solver has not reached. Nothing in this audit distinguishes those
two. What would:

1. **A boundary-only comparison at fixed interior.** Evaluate the same column
   under the molecular equation set with the base reservoir particle count
   forced to the atomic value, and again with the molecular value, everything
   else held. If the face inflow tracks `ntot_bc` as 1/0.950 predicts, the step
   is the reservoir definition and is a physical statement to be checked against
   the handoff, not a solver problem. This is a Mode R probe and needs no new
   state.
2. **Does a stationary solution exist for this base.** A Mode C continuation
   from `g0005` with a declared budget, watching the base face mass flux, which
   has already moved three decades toward the wind flux and stalled at -7.36 of
   it. rev9 section 16 already lists the existence question as open and this
   audit does not close it.
3. **The term export of phase 5**, so that the energy row of cells 1 and 2 can
   be read as stored energy face fluxes plus `dF3p` rather than as the derived
   split of section 6.3.

---

## 10. Things noticed that qualify the plan or the code

Reported, not fixed; none of them is in this item's scope.

1. **rev9 section 11.2's reasoning about charge is right and is now measured.**
   The transfer moves mass only between two neutral species, and the census
   charge is preserved algebraically. Section 5 turns that from an argument into
   a number and supplies the floor an assertion should use.
2. **The conversion's header overstates what x2 = 0 does.**
   `molecular_seed_from_atomic_state.f90:118-120` says a stated fraction of 0 is
   "the conversion identity (nothing is transferred, so every species column,
   the nuclei and the chosen primitive invariants come back unchanged)". It is
   the identity over the physical cells and not over the two lower ghosts, which
   `load_IC` fills with the molecular reservoir row carrying H2 and which the
   transfer then empties into H I. MEASURED: with x2 = 0 the reported
   temperature displacement is 4.9884E-02 and every part of it is in those two
   rows.
3. **The seed mode is not recoverable from a state file once a solve has been
   taken.** `write_molecular_seed_header` writes its two lines only in the
   process that performed the conversion (`:849`, `if (.not. seed_built)
   return`), so no stored generation of this case carries them and the mode had
   to be recovered from two logs and a provenance record.
4. **The `handoff` mode is inadmissible on this column and says so quietly.**
   It caps 385 of 504 cells at their own neutral hydrogen and the evaluation
   then finds 98 cells with no chemical root and a largest accepted reaction
   residual of 2.784. The conversion itself reports only the cap count and
   raises nothing.
5. **The phase-5 export exists as source and is not in the binary of record.**
   `src/modules/time_step/conservation_budget.f90` is untracked, was being
   edited while this ran, and `EXHALE.x` md5 `a55d9008` does not contain it. The
   `Makefile` already lists it at line 145. Nothing was built here.

---

## 11. Commands

All raw products are in the worker scratch; nothing was written into the
repository except this file, and nothing was written into
`LHS1140b/models/<case>/`.

```
SCR=/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/
      59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/m2
E=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
ATOM=$E/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k04/output

# one evaluation of the checkpoint, with the ghost record and the boundary trace
env -i PATH=/usr/bin:/bin HOME=$HOME OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1 $E/EXHALE.x
#   in $SCR/R0_eval_g0005, output/*_IC.txt = the stored g0005 pair       exit 0

# the same state through the residual diagnostic, for the rows of every cell
env -i ... EXHALE_RESIDUAL=1 EXHALE_BOUNDARY_TRACE=1 $E/EXHALE.x
#   in $SCR/R1_resid_g0005                                              exit 0

# the four conversions (Restart intent: stationary, Load IC? True)
env -i ... EXHALE_MOLECULAR_SEED=$ATOM EXHALE_MOLECULAR_SEED_INVARIANT=p \
    [EXHALE_MOLECULAR_SEED_X2=local|0]  $E/EXHALE.x
#   in $SCR/S_{local_p,handoff_p,zero_p,local_T}                    all exit 0

# one evaluation of each converted seed
env -i ... EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1 $E/EXHALE.x
#   in $SCR/E_{local_p,handoff_p,zero_p,local_T}                    all exit 0

# the base face state of the seeds and of the atomic control
env -i ... EXHALE_RESIDUAL=1 $E/EXHALE.x
#   in $SCR/F_{zero_p,local_p,atomic}                               all exit 0

# the control: the atomic column in its own configuration
env -i ... EXHALE_GHOST_RECORD=1 EXHALE_BOUNDARY_TRACE=1 $E/EXHALE.x
#   in $SCR/C_atomic_k04, input = the k04 input.inp with absolute paths
#   and 'Restart intent: stationary evaluate'                           exit 0
```

md5 of the four converted seed pairs (MEASURED):

```
local_p     Hydro b8fe4f3b557acce79b249d0e8507f50a  Ion e7101b2692174fea582ed5c710eba8c1
handoff_p   Hydro 9a5c53e83fee55b0a312c7038219fe2c  Ion ec3f94782d77bd7a0f56a71661b29ae0
zero_p      Hydro 97eae5c9e1204d9b4637aeff9c36bcf5  Ion 031cd879416bb1a05329fa72a20918f3
local_T     Hydro fc0d9553606530270006af0ef65798e6  Ion d1ecc03a853a160f3707d955138327b6
```

Products used: `output/residual_profile.txt`, `output/boundary_trace.txt`
(face-state gaps and the base mass rows of cells 1 to 12),
`output/ghost_record.txt` (the lower ghost and the base continuity row it
produces), `output/Hydro_ioniz.txt` and `output/Ion_species.txt` of each
evaluation, and the `(molecular_seed)` and `(certification)` blocks of each log.
The charge and budget arithmetic of section 5 was done from the written species
files alone, with `census.py` in the same scratch directory.

No background process of this work is left running.

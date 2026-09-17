# LHS 1140 b, item L10: the post-processing `_adv` profiles in the outer wind

Item L10 of `docs/PLAN_20260913_lhs_stationary.md`, opened by section 7 of
`docs/lhs1140b_stationary_L8_20260913.md`: on the certified 45 R_p wind the
advection-corrected temperature stands at 1199 K where the solution carries
339 K, and the 45 and 60 R_p `_adv` profiles disagree with each other above
30 R_p while the two solutions agree to 0.1 percent (READ from that
report; measured here as 0.28 percent over 2 to 44 R_p). This states which
correction moves the temperature, whether the correction is defined where it
is applied, what the He I 10830 line inherits from it, and what the tools
that read `_adv` should read instead.

Every number below is MEASURED unless marked READ. The states are the three
certified L8 solutions in `LHS1140b/models/.L8/{r30,r45,r60}`; the binary
that wrote them is `EXHALE.x`, md5 `97e10317a710b9ccc63addbedde3586a`, and
the control binary of section 4 -- the private build with only the entry
text of `post_process_adv.f90` in it -- reproduces every one of their `_adv`
files bit for bit.

---

## 1. Verdict

**The temperature split is not a post-processing defect and not an
undefined correction outside the exobase. It is the advection correction
doing exactly what it is built to do, and what it reports is that the
photoionization-equilibrium composition the wind was SOLVED with is wrong
above about 4 R_p, by up to a factor 4.5 in the neutral fraction and 5.9 in
the photoheating. The `_adv` temperature is the temperature that heating
implies; the solution's temperature is the one its own (unrelaxable)
equilibrium composition implies. Neither is a self-consistent state,
because the density and the velocity were never re-solved for the corrected
composition.**

The chain, each link measured:

1. The gas cannot relax. Over a flow time `r/v` the hydrogen recombines
   `t_flow/t_rec = 1.7` times at 2 R_p and 5e-3 times at 29 R_p, and is
   photoionized `t_flow/t_ion = 7e-2` times at 2 R_p and 1e-2 at 29 R_p
   (section 3.3). Above about 4 R_p both numbers are 1e-2 or smaller: the
   ionization state is frozen into the parcel and carried outward.
2. The correction says so. The equilibrium neutral fraction falls from
   0.963 at 2 R_p to 0.346 at 29 R_p; the advection-corrected one goes from
   0.974 to 0.941, that is, hardly at all.
3. Photoheating is proportional to the neutral density, so the corrected
   heating stands above the run's own by a factor that tracks the neutral
   one and runs a little ahead of it, the helium being more neutral too:
   1.0 at 4 R_p, 1.5 at 10, 2.7 at 20, 3.7 at 29 and 5.9 at 44 R_p, against
   a hydrogen neutral-fraction ratio of 1.1, 1.4, 2.0, 2.7 and 4.5 at the
   same radii.
4. The corrected temperature follows the corrected heating. Columns
   `T_adv/T` and `Q_adv/Q` of the table in section 2 track each other over
   the whole factor 5.9, agreeing with each other to 7 percent at every
   radius.

**One genuine defect was found on the way and is fixed.** The 45-against-60
disagreement of section 7 of the L8 report is ONE discarded root: MINPACK
returns `info = 4` ("no further progress") at three cells of the 45 R_p run
even though the iterate it returns satisfies the energy equation to
`|R|/s = 1.5e-17 ... 5.2e-17` of the largest term the equation holds, and
the code threw those roots away. Because the correction is an upwind
recursion, the discard restarted the temperature profile at r = 32.58 R_p
and every row above it. With the roots accepted by their residual, the 45
and 60 R_p `_adv` temperatures agree to 0.19 percent over 2 to 44 R_p -- the
0.28 percent at which the two SOLUTIONS agree. Sections 4 and 6.

**The He I 10830 equivalent width of the L8 report is not touched by the
temperature.** Feeding the transit tool the solution's temperature with the
`_adv` composition moves the 30 R_p equivalent width from 1.1702 to 1.1685
percent A (0.15 percent); feeding it the `_adv` temperature with the
equilibrium composition moves it to 16.1466 (a factor 13.8). The line is a
statement about the composition and about nothing else in these two files
(section 5).

---

## 2. Where the two profiles part, and by how much

`T` is the solution, `T_adv` the advection-corrected profile, `x_HI` the
neutral fraction of the two, `Q_adv/Q` the ratio of the photoheating the
post-process computes to the run's own, `Da_rec` and `Da_ion` the
recombination and photoionization Damkohler numbers of section 3.3. MEASURED
on `models/.L8/r45`.

| r [R_p] | T [K] | T_adv [K] | T_adv/T | x_HI | x_HI,adv | Q_adv/Q | Da_rec | Da_ion |
|---|---|---|---|---|---|---|---|---|
| 1.99 | 4111 | 4004 | 0.97 | 0.963 | 0.974 | 0.96 | 1.7e+00 | 6.7e-02 |
| 4.02 | 2038 | 2087 | 1.02 | 0.886 | 0.961 | 0.98 | 1.9e-01 | 2.4e-02 |
| 5.97 | 1368 | 1550 | 1.13 | 0.804 | 0.956 | 1.10 | 7.2e-02 | 1.8e-02 |
| 8.06 | 1018 | 1312 | 1.29 | 0.736 | 0.953 | 1.29 | 3.7e-02 | 1.3e-02 |
| 9.98 | 831 | 1208 | 1.45 | 0.688 | 0.951 | 1.51 | 2.4e-02 | 1.1e-02 |
| 15.00 | 574 | 1120 | 1.95 | 0.580 | 0.948 | 2.09 | 1.2e-02 | 8.6e-03 |
| 20.17 | 447 | 1126 | 2.52 | 0.478 | 0.945 | 2.69 | 7.7e-03 | 8.4e-03 |
| 24.93 | 377 | 1161 | 3.08 | 0.399 | 0.943 | 3.26 | 5.8e-03 | 8.8e-03 |
| 28.74 | 339 | 1199 | 3.54 | 0.346 | 0.941 | 3.74 | 4.9e-03 | 9.2e-03 |
| 44.22 | 252 | 838 | 3.33 | 0.208 | 0.937 | 5.89 | 3.0e-03 | 1.1e-02 |

The last row carries the discarded root of section 4: with it accepted,
`T_adv` at 44.22 R_p is 1393 K and `T_adv/T` is 5.53. The composition and
the heating of that row do not move (`x_HI,adv` 0.937 and `Q_adv/Q` 5.89
either way), which is itself the statement that the discard was in the
energy solve alone.

The 30 R_p run, which is the one the campaign quotes, carries the same
column to the same values: `T_adv/T` = 0.98, 1.03, 1.14, 1.28, 1.46, 1.95,
2.50, 3.05, 3.53 at 2.00, 4.01, 5.99, 7.95, 10.03, 14.96, 20.10, 24.96 and
29.03 R_p, with `Q_adv/Q` = 0.96, 0.98, 1.10, 1.28, 1.51, 2.07, 2.63, 3.15,
3.57 at the same radii. **The split does not begin at the exobase. It
begins where the ionization stops relaxing, at about 4 R_p, and grows
monotonically outward.** Every row of all three runs carries
`adv_T_status = 0` (corrected) outside the eight innermost rows, and
`adv_comp_status = 0` outside the innermost 28, all of which sit below
1.005 R_p; the only other nonzero row is the single discarded root of
section 4. Outside 2 R_p the mass row of the state runs between 3.1e-09 and
2.9e-03 over the three runs, below the fraction `adv_conditional_tol = 1e-2`
at which a row is retained. **None of the validity conditions the
post-process states fires anywhere in the split.**

---

## 3. Which correction moves the temperature

### 3.1 The corrected profile satisfies the equation it states

`post_process_adv.f90` writes the steady internal-energy balance it solves,
upwind-differenced, as

```
rho v de/dr  -  (p/rho) v drho/dr  +  h div(rho v)  =  heat - cool .
```

Re-assembled outside the code from the columns of the files themselves, with
`e = (3/2) p/rho`, `h = e + p/rho`, and `div(rho v)` formed as
`(F_j - F_{j-1}) / [(r_j^3 - r_{j-1}^3)/3]` with `F = rho v r^2`, the
residual of that equation on the `_adv` profile is

| r [R_p] | 4.98 | 9.98 | 20.17 | 28.74 | 35.00 | 44.22 |
|---|---|---|---|---|---|---|
| residual / largest term, `_adv` | 6.3e-05 | 1.0e-06 | 4.3e-08 | 1.5e-08 | 1.4e-14 | 7.5e-11 |
| residual / largest term, solution | 3.2e-02 | 3.1e-02 | 2.9e-02 | 2.9e-02 | 2.8e-02 | 1.1e-01 |

The corrected profile satisfies the equation the post-process states to
6e-5 at 5 R_p and to better than 1e-6 everywhere outside 10 R_p. The
solution satisfies the same equation to about 3 percent at every radius,
which is the size expected of the difference between the conservative
total-energy finite-volume form the wind was solved in, with PLM
reconstruction, and this first-order upwind non-conservative internal-energy
form. **A discretization difference of three percent, not of a factor 3.5:
the factor 3.5 is not in the discretization.**

### 3.2 The enthalpy flux of the mass-flux divergence is not the term

That third term is the one the module header warns about, and it is not what
moves the temperature here: `q_enth / max(q_adv, q_prs)` is 2.6e-05 at
4.98 R_p, 6.5e-05 at 9.98, 5.7e-05 at 20.17 and 5.5e-05 at 28.74 R_p on the
certified state, because `rho v r^2` formed at the cell centers changes by
about 1e-6 of itself from one cell to the next. (It reaches 5.1e-02 in the
last physical cell, where the outer ghost enters the difference, and that is
still two orders below the terms it is compared with.) The post-process
reports its own maximum of that ratio as 1.57e+03 at r = 1.0002 R_p, in the
ghost region, dominant in 7 of 503 rows, all of them there (READ from
`pp.log` of all three runs).

### 3.3 The correction is the ionization advection, and the gas cannot relax

What is left is the composition. The advection correction replaces the
run's photoionization equilibrium, solved cell by cell, by the steady
species balance with the advection term kept, and in this wind that is the
whole answer.
The two timescales the equilibrium assumption needs, against the flow time
`r/v`:

| r [R_p] | t_flow [s] | t_rec [s] | t_ion [s] | t_flow/t_rec | t_flow/t_ion |
|---|---|---|---|---|---|
| 1.99 | 7.2e+05 | 4.1e+05 | 1.1e+07 | 1.7e+00 | 6.7e-02 |
| 4.02 | 2.8e+05 | 1.5e+06 | 1.2e+07 | 1.9e-01 | 2.4e-02 |
| 8.06 | 1.8e+05 | 5.0e+06 | 1.4e+07 | 3.7e-02 | 1.3e-02 |
| 15.00 | 1.7e+05 | 1.4e+07 | 2.0e+07 | 1.2e-02 | 8.6e-03 |
| 29.79 | 2.0e+05 | 4.3e+07 | 2.2e+07 | 4.7e-03 | 9.3e-03 |
| 44.22 | 2.5e+05 | 8.3e+07 | 2.2e+07 | 3.0e-03 | 1.1e-02 |

`t_rec = 1/(n_e alpha_B)` with `alpha_B = 2.59e-13 (T/1e4)^-0.7` (Case B
power-law fit, used here only to set a timescale), and `t_ion = 1/P_HI` with
`P_HI` read back from the run's own equilibrium,
`P_HI = n_HII n_e alpha_B / n_HI`, so the photoionization rate is the run's
and not a second estimate of it. Above 4 R_p a parcel is ionized or
recombines less than three percent of the way to equilibrium while it
crosses the wind. **The local root the run puts the cell on is not the
composition of that gas, and the advected one is.**

This is not a new statement about the code; it is written in
`src/modules/init/parameters.f90` beside the `ionization_transport` key:
"the gas leaves each shell three to seven times faster than it can be
photoionized, so the ionization fraction is not the local root but whatever
the parcel accumulated on the way up" (READ). What is new is that on this
planet the effect is an order of magnitude larger than the hot-Jupiter case
that comment was written from, and that the option which would put the
transported proton inside the solve is refused for an atomic run:
`input_read.f90` requires `Molecular chemistry: True` and
`Molecular carrier transport: True` before `Ionization transport: True` is
accepted (READ). **The LHS 1140 b wind as configured has no self-consistent
route to the composition its own post-process says it has.**

---

## 4. The defect: a root discarded because the iteration stopped improving

### What was wrong

The scalar energy equation of each cell is solved by MINPACK `hybrd1`, and
the code kept the iterate only on `info = 1`:

```fortran
if (info /= 1) then
   sys_x_T(1) = T_in(j)
   n_T_noconv = n_T_noconv + 1
```

`info` is a statement about how the ITERATION ended, not about the iterate.
`info = 4` and `5` are returned when the steps stop improving the residual,
which is what a stalled search and an arrived one look like alike: at the
root the residual sits at the cancellation floor of the terms it is
assembled from and no step can lower it.

### Reproduction (RED)

The 45 R_p post-processing pass, re-run from the archived state with the
delivered binary, reproduces its `_adv` files bit for bit in 1.8 s and
reports `3 cell solves did not converge and kept the equilibrium
temperature`. Instrumented, those three cells are

| loop index j | r [R_p] | info | T_in [K] | iterate [K] | residual | largest term | ratio |
|---|---|---|---|---|---|---|---|
| 312 | 2.3124 | 4 | 3537.8 | 3454.9 | -4.96e-23 | 3.41e-06 | 1.5e-17 |
| 442 | 15.5253 | 4 | 557.3 | 1118.3 | -6.14e-25 | 2.21e-08 | 2.8e-17 |
| 484 | 32.5780 | 4 | 309.2 | 1241.4 | 2.84e-25 | 5.47e-09 | 5.2e-17 |

(`j` is the loop index of the module, which runs from `1-Ng`; the last of
the three is the row that carries `adv_T_status = 2` at r = 32.578 R_p in
`Hydro_ioniz_adv.txt`.)

Every one of the three iterates is the root to full double precision. The
event is at the rounding boundary: a private build of the same source with
the system gfortran 13 instead of the conda-forge gfortran 16.2 does not hit
it in any cell, and its 45 R_p profile is then the one the fix produces.

### The fix

`src/modules/post_process/post_process_adv.f90`, the `info /= 1` branch of
the legacy solve: the iterate is kept when its residual is negligible
against the largest term the equation holds,

```
s   =  max( |mu_up  rho v T|                  advected energy of the cell
            |mu_cell rho v T_up|              what the flow carries in
            |(gamma-1) coeff T|               compression work
            |(gamma-1) mu_cell mu_up dr Q| )  photoheating
|R| <=  1e2 * epsilon * s
```

with `mu_up = mmw(j-1)`, `mu_cell = mmw(j)` and
`coeff = mu_up v (rho_j - rho_{j-1})`: the four terms `T_equation`
assembles, in its own variables. A sum of terms of
size `s` cannot be formed to better than a few machine epsilons of `s`, so
`1e2 * epsilon` is the level at which the equation is an identity in double
precision, and an iterate that is NOT a root stands orders of magnitude
above it (the three above are at 1e-17, fifteen orders below `1e2*epsilon =
2.2e-14`). Kept iterates are counted and reported on their own line; the
"did not converge" count now counts only cells that really did not.

### GREEN, and what it moved

| case | control build | with the fix |
|---|---|---|
| `.L8/r30` PP | no failed solve | byte-identical output |
| `.L8/r45` PP | 3 solves discarded | 0 discarded, 2 roots kept; output differs |
| `.L8/r60` PP | no failed solve | byte-identical output |

The control binary is the measured build with one object swapped: the
delivered object set of the private build, relinked with
`post_process_adv.o` recompiled from the entry text of the file. It
therefore differs from the measured binary in that one file and in nothing
else, which matters because other work is landing in this tree: two full
builds taken 20 minutes apart already march a WASP-121b case apart at the
1e-9 level, and a comparison across them would not have been a comparison of
this change. The control reproduces all three delivered `_adv` files bit for
bit.

What the 45 R_p profile does after the fix:

| r [R_p] | T solution | T_adv control | T_adv fixed | T_adv of the 60 R_p run |
|---|---|---|---|---|
| 28.74 | 339 | 1198.60 | 1198.60 | 1199.05 |
| 32.58 | 309 | 309.21 | 1241.41 | 1237.13 |
| 35.00 | 294 | 447.33 | 1270.16 | 1266.72 |
| 39.70 | 269 | 665.46 | 1328.33 | 1326.94 |
| 44.22 | 252 | 838.05 | 1392.91 | 1387.54 |

Largest `|T_adv(45)/T_adv(60) - 1|` over 2 to 44 R_p: **0.751 before, 0.0019
after**, against **0.0028** for the two solutions over the same interval.
**The `_adv` post-process is as boundary-independent as the solve; it only
looked otherwise because of the discarded root.**

### What it moves in the regression matrix

**Nothing, and the reason is structural before it is measured.** The branch
edited here is the legacy MINPACK route, taken only when
`pp_metal_on .and. use_brent_tsolve` is false. `wasp_full` carries a
`metals.inp` with `pp_metals = 1`, so `pp_metal_on` is true, and
`use_brent_tsolve` is true unless `Brent solver: False` is stated, which that
input does not state: its post-process energy solve takes the bracketing
route and never reaches the branch. The same holds for every metals-on case
of the matrix. The metals-off cases reach the branch and change only if a
cell there returns `info /= 1`; none of the three LHS 1140 b runs outside
`r45` does, and their output is byte-identical across the two builds.

MEASURED on a WASP-121b state carrying that physics: `wasp_full`'s own
`input.inp` and `metals.inp` marched to a fixed 2000 steps, and the
post-processing pass then run twice on that one state, once with each
binary. `output/Hydro_ioniz_adv.txt`, `output/Ion_species_adv.txt` and
`output/OI_levels_adv.txt` are BYTE-IDENTICAL, and neither run printed an
energy-equation line at all: with metals on the branch is never entered,
as the paragraph above says it cannot be.

`src/tests/adv_static_limit` is the only suite whose driver links
`post_process_adv.o` (grep over `src/tests/*/run.sh`). Run against both
private object directories it reports 53 rows and `PASSED` in each, with
identical row verdicts. The regression matrix was not run and no golden was
refreshed.

A second consequence of the same discard: `collisional_validity.py` read the
restarted cold step at 32.58 R_p as a sound-speed drop and reported a
critical point at 32.5 R_p for the 45 R_p run. On the fixed `_adv` profile
there is no critical point in the domain at all (max Mach 0.601), and on the
solution there is one at 40.055 R_p, which is the value the L8 report
measured from the solved state.

---

## 5. What the He I 10830 line inherits

`EXHALE_transit.py` reads `output/Hydro_ioniz_adv.txt` and
`output/Ion_species_adv.txt` with no key or environment variable to choose
otherwise (READ, lines 79-80). The question L8 raises -- how much of the
30 R_p equivalent width of 1.1702 percent A is `_adv` contamination -- is
answered by running the tool on four combinations of the two files, in
scratch copies where the `_adv` names carry the substituted content. Same
binary, same WINERED HIRES-Y kernel (`R = 68000`), same run directory.

| fed to the transit tool | EW [%A] | red depth [%] | FWHM [A] |
|---|---|---|---|
| A `_adv` temperature + `_adv` composition (the L8 number) | 1.1702 | 4.326 | 0.2538 |
| B solution temperature + solution composition | 15.5116 | 59.425 | 0.2462 |
| C solution temperature + `_adv` composition | 1.1685 | 4.469 | 0.2452 |
| D `_adv` temperature + solution composition | 16.1466 | 58.184 | 0.2623 |

**The temperature is worth 0.15 percent of the equivalent width (A against
C, the composition held fixed), 3.3 percent of the depth and 3.4 percent of
the width. The composition is worth a factor 13.8 at the `_adv` temperature
(A against D) and 13.3 at the solution's (C against B).**

The He 2^3S column itself is 0.164 of the equilibrium column in all three
runs (6.76e10 against 4.14e11 cm^-2 at 30 R_p), because the metastable is
fed by recombination and the advected gas has fewer free electrons to
recombine: `n_e(adv)/n_e(solution)` is 0.73 at
2 R_p, 0.32 at 5, 0.20 at 10, 0.15 at 20 and 0.13 at 29 R_p, and the He 2^3S
density itself follows it down, 0.61, 0.11, 0.031, 0.016, 0.013 at the same
radii.

So: **the L8 equivalent width carries no measurable contamination from the
`_adv` temperature, and the item's premise on that point is settled in the
negative.** What it does carry is the whole weight of the choice of
ionization closure, and section 3.3 measures that choice: the advected
composition is the one the timescales of this wind permit and the
equilibrium one is not. Read the other way: at this run's He/H, a line built
from the solved composition would have been an equivalent width of
15.5 percent A against a measurement of 1.108 +/- 0.030 (Cherubim et al.
2026, READ from `LHS1140b/MODELS.md`). The campaign crosses He/H to meet
that measurement, so the closure does not only shift a number, it moves the
He/H the crossing returns.

---

## 6. What the tools should read, and the change made

The two files are two different states. `Hydro_ioniz.txt` and
`Ion_species.txt` are a state whose density, velocity, pressure and
temperature satisfy the momentum equation the wind was solved with.
`*_adv.txt` carries that state's density and velocity beside a temperature
and a composition from a second closure, for which the momentum equation was
never re-solved. Anything formed by mixing them belongs to neither.

- **Momentum-equation quantities -- sound speed, Mach number, critical
  point, the structure scale a Knudsen number is divided by, the mass flux
  -- must come from the solution.** Measured cost of not doing so: no
  critical point in the domain from `_adv`, 40.055 R_p from the solution.
- **The composition is the one question `_adv` answers better**, and the
  mean free path depends on it: `Kn_bulk = 1` sits above the outer boundary
  on the solution's composition and at 29.0 R_p on the advected one, and
  `Kn_bulk = 0.1` at 12.91 against 7.65 R_p. The two bracket the exobase;
  they do not agree, and no file in the run holds the state that would
  settle it.
- **`EXHALE_transit.py` should go on reading `_adv`** and needs no change:
  section 5 measures its temperature sensitivity at 0.15 percent of the
  equivalent width, and its composition -- the quantity that sets the line
  -- is the one the `_adv` file gets right.

`src/utils/collisional_validity.py` is changed accordingly: the SOLUTION is
now the default and `--adv` is the opt-in (the flag was `--eq` with the
opposite default). The docstring gains a section "Which state this reads"
carrying the reason and the two measured numbers, and the report header line
now names the state it read and what that state's composition is. No
formula, constant, or diagnosis in the file is touched.

Nothing else in the repository calls the script (checked across `*.py`,
`*.sh`, `*.f90`, `*.md`), so the flag rename breaks no caller. The `_adv`
numbers quoted in `docs/lhs1140b_width_nonthermal_candidates.md` section on
the collision diagnostic and in section 7 of the L8 report were produced
with the old default and are `--adv` numbers; they are labelled as `_adv` in
both places and stay correct.

---

## 7. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
SC=$EX/LHS1140b/models/.L10

# the post-processing pass alone, from an archived certified state
cp -r $EX/LHS1140b/models/.L8/r45 $SC/r45_x && cd $SC/r45_x
rm -f output/Hydro_ioniz_adv.txt output/Ion_species_adv.txt
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
       -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
OMP_NUM_THREADS=8 $EX/EXHALE.x > pp.log 2>&1     # 1.8 s

# the transit line from a chosen pair of profiles: substitute under the
# _adv names, the only pair the tool reads
cp -f output/Hydro_ioniz.txt output/Hydro_ioniz_adv.txt
. $EX/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=$EX python3 $EX/EXHALE_transit.py > transit.log 2>&1

# the collision diagnostic, both states
python3 $EX/src/utils/collisional_validity.py $SC/r45_x          # the solution
python3 $EX/src/utils/collisional_validity.py --adv $SC/r45_x    # the _adv pair
```

The private build and its control, which differ in one object:

```bash
cd $EX
make OBJDIR=build_L10 EXE=EXHALE_L10.x            # the measured build
cp -a build_L10 build_L10ctl2
gfortran -O3 -fopenmp -Jbuild_L10ctl2 -Ibuild_L10ctl2          -c <entry text of post_process_adv.f90> -o build_L10ctl2/post_process_adv.o
# relink the same object list into EXHALE_L10ctl2.x
```

The WASP-121b impact measurement, on one state and two binaries:

```bash
# march wasp_full's own input to a fixed step count once
( cd <dir with wasp_full input.inp + metals.inp> &&   OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=2000 $EX/EXHALE_L10.x > run.log 2>&1 )
# then the post-processing pass on that state with each binary, as above
```

The energy residual of section 3.1, the Damkohler numbers of section 3.3 and
the equivalent-width integral of section 5 are formed from the file columns
alone; the equivalent-width window is 10832.60 to 10834.20 A in air, the one
`LHS1140b/models/run_case.sh` uses.

---

## 8. Noticed outside the scope, not acted on

- **The same `info /= 1` fallback governs the advection ionization system
  and the metal stage re-solve** (`n_adv_noconv`, `n_metal_noconv` in the
  same file), where the unknowns are a vector and the same argument applies
  to a residual norm against a row scale. Neither fired in any of the three
  LHS 1140 b runs, so there is no failing case to test a change against and
  none was made.
- **The `_adv` energy correction carries an unstated amplitude assumption.**
  Two windows say, without saying so, that the corrected temperature is
  expected to stay within a factor of a few of the run's: the bracketing
  route, which is the default whenever metals are on, scans only
  `[0.05, 4] x T_in` (`T_equation.f90`, `solve_T_brent`), and the MINPACK
  route with metals on additionally rejects any root outside
  `[0.5, 2] x T_in` (`post_process_adv.f90`). MEASURED on this wind,
  `T_adv/T_in` crosses 2 at about 16 R_p and 4 at about 33 R_p, so the same
  planet run with metals on would lose its correction above 33 R_p on the
  default route, and above 16 R_p with `Brent solver: False`, in both cases
  reverting to the run's temperature and reporting only a count. The
  scanning window lives in a file this item does not own.
- **The certified LHS 1140 b solutions are built on a composition their own
  post-process contradicts by a factor 4.5 above 4 R_p** (section 3.3), and
  the key that would close the loop, `Ionization transport`, is refused for
  an atomic run. This is a statement about the campaign, not about the
  post-process, and it belongs to whoever owns the campaign's conditions.

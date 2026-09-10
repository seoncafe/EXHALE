# Audit E: main program, input/output, restart, lower-profile handoff, post-processing, Python tools

Repository root for every path below: `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/`.

---

### The final output file claims the secondary-ionization coupling is on while its own heat/cool columns were produced with it off

Status: CONFIRMED
Severity: P1
Location: `src/EXHALE_main.f90:1782-1787`, consumed by `src/modules/functions/utilities.f90:180-198` (`write_coupling_state_header`) and `src/EXHALE_main.f90:1875-1879` (final `write_output`)

```
1782      ! Ensure the post-processing and final-output rates carry the full physics
1783      ! even if the loop hit the iteration cap before the stage flip fired.
1784      if (use_sec_ion) then
1785         if (.not. sec_ion_active) sec_ion_armed_step = count
1786         sec_ion_active = .true.
1787      endif
```

What the code does: after the marching loop it forces `sec_ion_active = .true.` and stamps
`sec_ion_armed_step = count`. It does NOT re-solve `ioniz_eq`, so the `heat` and `cool` arrays
handed to the final `write_output` are still the ones the last in-loop sweep produced with the
coupling OFF. Meanwhile `write_coupling_state_header` reads the global `sec_ion_active` at write
time (`s = 'F'; if (sec_ion_active) s = 'T'`, utilities.f90:183) and so writes
`# coupling: sec_ion=T sec_ion_step=<last step>` into that same file, and
`write_heat_breakdown_eq`, `write_cool_breakdown_eq` and `post_process_adv` all run afterwards
with the coupling armed.

What it should do: the `heat`/`cool` columns, the channel breakdowns and the coupling header must
all describe one state. Either the sweep is redone once with the coupling armed before the final
write, or the header records the coupling the state was actually relaxed under (which is what the
header exists for: utilities.f90:142-153 states the case, and `load_IC` acts on it at
`src/EXHALE_main.f90:846-852`).

How it was verified: read the call order from the loop exit (line 1780) to `write_output`
(line 1878); nothing between them calls `ioniz_eq`. Measured on the shipped golden of the
`mol_base_handoff` case, whose header reads `# coupling: sec_ion=T sec_ion_step=12000`
(= the step cap, i.e. the coupling never flipped during marching):

```
max  |heat(Hydro_ioniz.txt) - heat_total(Heating_breakdown.txt)| / heat = 8.23e-01  at r = 1.1148
median over the 500 physical cells                                     = 7.23e-01
```

so the two files disagree by a median of 72 percent on the same state. The direction and size are
what the Shull and van Steenberg heating fraction gives in weakly ionized gas.

Reaches: every run that leaves the loop without the staged flip having fired, i.e. every run
stopped by `EXHALE_MAXSTEPS` or `count_max`, or by a NaN. That is 8 of the 11 default regression
cases (all the `mol_*` cases and `lower_profile` are 12000-step snapshots; `run_check.sh`'s own
comment for `mol_sec_ion` says the flip never fires in them). The goldens themselves therefore
carry the inconsistent pair. The restart consequence is the one the header was introduced to
prevent, inverted: `Hydro_ioniz.txt` of such a run, fed back as `Hydro_ioniz_IC.txt`, tells the
next run the state was produced with the coupling armed when it was not, so the restart adopts a
photoionization rate the state was never relaxed under.

---

### A restart overwrites the cell-centre grid with the radii stored in the IC file, leaving `r_edg`, `dr_j` and `j_min` from the run's own grid

Status: CONFIRMED
Severity: P1
Location: `src/modules/files_IO/load_IC.f90:751` (schema-2 path, the default) and `:214` (legacy path)

```
751      r(j) = vals(1)
```
```
214            read(2,*) r(j),                  &
215                      nsp_l(j,isp_HI),       &
```

What the code does: `r` is the global cell-centre array of `global_parameters`, built by
`define_grid` before `load_IC` runs. `scatter_row` (and, on the legacy path, the read statement
itself) assigns the first column of `Ion_species_IC.txt` into it, for every row including the four
ghost rows. `r_edg`, `dr_j`, `j_min` and `j_flux` are NOT reassigned: they keep the values
`define_grid` computed from `Grid type`, `Base grid`, `Grid cells` and `r_max`.

What it should do: the grid is a property of the run, not of the file. A restart should either
leave `r` alone or check the file's radii against it and refuse a mismatch. The current code
silently accepts a file written on a different grid.

How it was verified: read `scatter_row` end to end (no local `r` is declared, so `r` is the
module-associated global) and traced the only guard `load_IC` applies to a foreign file, the row
count `nrec .ne. N + 2*Ng` (`load_IC.f90:157-166`). That guard catches a different `Grid cells:`
and nothing else: two runs with the same `N` but different `Grid type`, `Base grid` or
`Outer radius` pass it, and the restart then integrates on cell centres from one grid and cell
faces, widths and window indices from another. On a same-grid restart the round trip is exact
only because the writer uses list-directed output, which gfortran emits with 17 significant
digits (verified against `backup/regression/golden/wasp_full/Hydro_ioniz.txt`).

Reaches: every run with `Load IC? True`, including the `roundtrip` case. No default matrix case
compares a cross-grid restart, so nothing guards it.

---

### With a lower-atmosphere profile the base LEVEL the run uses is never checked against the profile's `p_match_bar`, and the shipped regression case runs 20 times away from it

Status: CONFIRMED
Severity: P1
Location: `src/modules/files_IO/input_read.f90:1652` (the check that is skipped) and
`:2326` (`apply_lower_atmosphere_profile`)

```
1652   if (base_level_from_base_inp .and. .not. lap_in_use) then
```
```
2326   p_base_bar = lap_p_match_bar
2327   write(*,'(A,ES9.2,A)') '   profile: p_base -> ', p_base_bar, ' bar'
```

What the code does: `apply_lower_atmosphere_profile` adopts the profile's matching pressure as
`p_base_bar` and prints it as the base level. The consistency check that refuses a base level
stated twice (`input_read.f90:1652-1741`, tolerance 1 percent) is explicitly disabled when a
profile is in use, with the comment "A lower-atmosphere PROFILE owns the level itself (it sets
p_base_bar from its own matching pressure), so a profile run is outside this rule." But `n0`, the
quantity that actually places the lower boundary, is then taken from
`Log10 lower boundary number density` and nothing compares the two. `p_base_bar` has only two
consumers (`composition.f90:288`, the chemical-equilibrium H2 fit, and the two reports), so on a
profile run the level the profile states and the level the wind uses are independent numbers.

What it should do: the base is one level. Either `n0` follows from `p_match_bar` the way it
follows from `base.inp`'s `p_base`, or the density key is refused beside a profile, or the two are
compared with the same 1 percent rule the scalar path uses. The profile's `T`, `r`, `q_H2` and
every elemental ratio are read AT `p_match_bar`, so a base placed elsewhere is a base whose state
was taken from a different level.

How it was verified: traced every assignment to `n0` in `input_read` (lines 1163, 1735, 1748) and
confirmed none of them is reached on the profile path without `Base BC: pressure`. Then measured
the shipped `backup/regression/lower_profile` case, whose own run log contains both statements:

```
run.log:10        profile: p_base ->  1.00E-06 bar
EXHALE_setup.out:38  - Base level: n0 =  1.0000E+14 cm^-3 (from the density key) -> p =  2.0035E-05 bar
```

a factor 20.03. At the profile's own 2.0e-5 bar level the table gives `r` = 1.372 R_J and
`q_H2` around 0.7, against the `r` = 1.401 R_J, `T` = 1450 K and `q_H2` = 0.1196 the run adopted
from the 1e-6 bar row. The case README asserts "the matching-level base state taken from the table
at p_match_bar = 1e-6 bar: T0 = 1450.0 K, R0 = 1.401 R_J, `p_base`, ...", so the intent is that
the profile fixes the level; the code does not.

Reaches: `Lower atmosphere profile:` runs, the `lower_profile` regression case and its golden, and
`EXHALE_resolved.out`'s `lower_profile_p_match_bar`, which `element_budget.py` and the closure
driver read as the handoff level.

---

### The `du` convergence functional is the spread of |rho v r^2|, not of rho v r^2

Status: CONFIRMED (functional mismatch), SUSPECTED (false stop)
Severity: P1
Location: `src/EXHALE_main.f90:1309-1316`

```
1309            mom_max = maxval(abs(mom(j_min:N)))
1310            mom_min = minval(abs(mom(j_min:N)))
...
1316            du = abs((mom_max-mom_min)/max(mom_min, 1.0d-30))
```

What the code does: takes the absolute value of the signed mass flux `rho*v*r^2` BEFORE the max and
the min, and normalizes by the minimum rather than by the mean.

What it should do: the quantity every comment, the setup report, `docs/input_schema.md` K13 and the
project notes describe is "the radial spread of `rho*v*r^2`", and the companion diagnostic
`flux_spread` twenty lines below (`EXHALE_main.f90:1368-1369`) computes exactly that, on the same
window, from the signed values normalized by the mean. Two different functionals are given one
name. Taking |F| first makes the metric blind to a sign reversal: a window carrying `+F` over half
its cells and `-F` over the other half returns `du = 0` while the flux is maximally non-constant.

How it was verified: evaluated both forms on the shipped goldens over each case's own
`r >= Escape radius` window:

```
case              r_esc   cells   sign changes   du (as coded)   signed spread / mean
lower_profile      2.00     99         2           8.189e+02          3.613e+00
mol_base_handoff   2.00    109         0           1.290e+00          6.961e-01
wasp_full          1.50     27         0           9.909e-04          9.906e-04
```

The `lower_profile` golden state does change sign twice inside the du window, and there the coded
form reads 227 times the signed spread. That direction is conservative (the stop cannot fire), so
I could not exhibit a state where the coded form is small while the signed spread is large; that
half is SUSPECTED. The functional mismatch itself is CONFIRMED and is present in the shipped
matrix.

Reaches: the default stop condition of every run that does not set `Solver: Newton` or
`Resid tol:`, and the arming test for the JFNK hand-off and the secondary-ionization flip
(`EXHALE_main.f90:1319-1330`, `1631-1636`).

---

### `j_min` is taken from a loop that starts in the ghost cells

Status: CONFIRMED
Severity: P2
Location: `src/modules/init/define_grid.f90:193-196`

```
193      do j = 1-Ng,N+Ng
194         if(r(j).ge.r_esc) goto 111
195      enddo
196
197 111   j_min = j
```

What the code does: searches for the first cell with `r >= r_esc` starting at `1-Ng`, so for an
`Escape radius` at or below the innermost ghost centre (`r(1-Ng)` is slightly below 1 for the
Mixed grid, exactly 1 for Uniform and Stretched) `j_min` is a GHOST index. Every window built on
`j_min` then includes boundary data: `du`, `dtu`, `flux_spread`, the level-stability mean
(`EXHALE_main.f90:1309, 1336, 1368, 1406`) and the wind/layer split of `residual_norms`
(`steady_residual.f90:563-564`).

What it should do: the convergence windows are statements about the solution, so the search should
start at cell 1. `write_output` emits a `# rows` header precisely so consumers drop the ghosts
(section 133.6), and `exhale_io` does; the Fortran windows do not have the same guard.

How it was verified: read `define_grid.f90:180-225` and the four `mom(j_min:N)` uses. The existing
guard only covers `j_min > N`, not `j_min < 1`. Not reachable with the escape radii of the shipped
cases (1.5 and 2.0); reachable with any `Escape radius [R_p]: <= 1.0`.

---

### `u(2,:) = u(2,:) + 1.0e-16` modifies the conserved momentum instead of guarding the division it names

Status: CONFIRMED
Severity: P2
Location: `src/EXHALE_main.f90:1333`

```
1333      	u(2,:) = u(2,:) + 1.0e-16 ! To avoid division by zero
1334
1335	      ! --- Infty-norm
1336	      dtu = max(maxval(abs(1.0-u(1,j_min:N)/u_old(1,j_min:N))), &
1337	      	    maxval(abs(1.0-u(2,j_min:N)/u_old(2,j_min:N))))
```

What the code does: adds 1e-16 (code units) to the momentum density of every cell, ghosts
included, once per step, in place. The division the comment names is by `u_old(2,...)`, which this
statement does not touch: the guard works only through the one-step delay, because `u_old` of the
next step is this step's biased `u`. As written it is an unphysical momentum source with no
corresponding mass or energy change, applied outside the operator split and therefore invisible to
the steady residual (`assemble_residual` never sees it).

What it should do: guard the denominator at the point of use, e.g.
`max(abs(u_old(2,j)), tiny)`, and leave the state alone.

How it was verified: read the statement order in the loop body (`u_old = u` at line 1005,
the retry block, then line 1333). Measured the size of the perturbation against the state on the
shipped goldens, in code units:

```
case               min |u(2)| over physical cells   1e-16 / min|u(2)|
wasp_full                    1.99e-03                   5.0e-14
mol_base_handoff             1.98e-06                   5.1e-11
lower_profile                1.18e-10                   8.5e-07
```

so the effect is small everywhere measured; the objection is that a state variable is modified by a
statement whose stated purpose is a division guard, and that the added momentum is not a term of
any equation the residual or the Newton solver contains.

Reaches: every marching run, every golden. Removing it would change every golden at round-off
level.

---

### `report_unmet_steady_gates` is defined and never called

Status: CONFIRMED
Severity: P2
Location: `src/EXHALE_main.f90:2606-2631`

Dead code: `grep -rn report_unmet_steady_gates --include=*.f90` returns only the definition and
its `end subroutine`. It duplicates part of what `report_marching_stop` (line 2564) prints.

---

### The Mdot cell index can fall outside the array for a small grid

Status: CONFIRMED
Severity: P2
Location: `src/EXHALE_main.f90:1943-1946`

```
1943      j = N - 20
...
1946      Mdot = log10(4.0*pi*rho(j)*v(j)*r(j)*r(j)*n0*v0*mu*R0*R0)
```

`input_read` accepts `Grid cells:` down to 10 (`input_read.f90:951-956`), and the arrays are bound
`1-Ng:N+Ng` with `Ng = 2`. For `N < 19` the index `N-20` is below the lower bound (N = 10 gives
j = -10 against a lower bound of -1). Also unguarded: `v(j) <= 0` makes the `log10` argument
non-positive. Not reachable with the shipped inputs (all use N = 500), so it is a robustness
defect, not a wrong number.

The Mdot expression itself is correct: `rho` is the mass density in units of `n0*mu`
(`utilities.f90:372-439`, `calc_rho` weights each species by `bsp_mass` in m_H), so
`rho*n0*mu*v*v0*4 pi (r R0)^2` is g/s.

---

### `get_word` returns an unallocated deferred-length result for a word that is not there, and the label matcher accepts separators `get_word` cannot split on

Status: CONFIRMED
Severity: P2
Location: `src/modules/files_IO/input_read.f90:2452-2551` (`get_word`), `:2361-2380`
(`lbl_match` / `is_sep`)

What the code does: `get_word` has `character(len=:), allocatable :: get_word` and returns without
assigning it whenever the requested word index exceeds the number of words on the line. Every
optional value in the parser depends on that path, testing `if (len_trim(str) .gt. 0)` afterwards
(`du_th` second value, `Base grid` cell count, `Shapiro filter` period, `Solver` switch,
`Low-Mach damping` Mach threshold, `Flux spread tol` radius, `Viscosity` exponent,
`Reconstruction continuation` fields). Referencing an unallocated allocatable is not conforming
Fortran.

Separately, `is_sep` accepts `=` and TAB as label separators, so `lbl_match` will match
`Planet mass=0.69` or a tab-separated line, but `get_word` splits on the blank character only
(`if (c_char .ne. '')` compares against a blank), so the value word is never found on such a line
and the subsequent `read(str,*)` fails or, for an optional key, silently does nothing.

How it was verified: extracted `get_word` verbatim into a standalone program and built it with
both compilers on this machine:

```
gfortran 16.2 (-fcheck=all):  missing word -> len(result) = 0
ifx 2026.1   (-check all)  :  missing word -> len(result) = 0
```

so both toolchains happen to hand back a zero-length string and no shipped case misbehaves. The
defect is that the parser relies on undefined behavior that neither compiler diagnoses. A
one-character line also returns an empty result (the terminating branch tests `counter == str_len`
after an increment, which cannot be reached from `counter = 1` when `str_len = 1`).

---

### `close(unit = 1)` on a unit `input_read` never opened

Status: CONFIRMED
Severity: P2
Location: `src/modules/files_IO/input_read.f90:1168`

`input_read` opens units 11 (input.inp) and `newunit` handles (base.inp), never unit 1. Closing an
unconnected unit is legal and has no effect, so this is a dead statement left from the positional
reader. Unit 1 is later used by `load_IC`.

---

### The post-process recomputes the photoheating with a one-pass-old ionized fraction, and says the opposite

Status: CONFIRMED
Severity: P2
Location: `src/modules/post_process/post_process_adv.f90:346-355` and `:820-822`

```
350	! to the hydrogen and helium nuclei").  Reused for both PH_heat calls
351	! below (nhi/nheii/nheiii unchanged between them).
...
355	xion = min(max(ne/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)
```

What the code does: `xion` is built at the top of each of the ten post-process passes and reused by
the second `PH_heat_HHe` call at line 820. Between the two calls the advection solve rewrites
`nhi`, `nhii`, `nhei`, `nheii`, `nheiii`, `nheiTR` (lines 592-690) and `ne` is recomputed at line
712, but `xion` is not. So the second call, which produces the `theat` written to
`Hydro_ioniz_adv.txt`, evaluates the Shull and van Steenberg photoelectron partition at the
ionized fraction of the PREVIOUS pass's state while passing the current pass's densities.

What it should do: either recompute `xion` from the updated `ne`, `nh`, `nhe` before the second
call, or state the lag. The quoted comment asserts the densities are unchanged between the two
calls, which is false.

How it was verified: read the pass body from line 340 to line 830 and listed every assignment to
`nhi`/`nheii`/`nheiii`/`ne`/`xion`. The lag vanishes at the fixed point of the ten-pass loop, but
the loop has no convergence test, so how much it has vanished is not measured anywhere.

Reaches: every `_adv` file of a run with helium. Not guarded: the regression harness never compares
`_adv` output (see below).

---

### The regression matrix compares no post-processed output at all

Status: CONFIRMED
Severity: P2
Location: `backup/regression/run_check.sh:157`

```
FILES="output/Hydro_ioniz.txt output/Ion_species.txt"
```

`backup/regression/golden/<case>/` holds exactly those two files plus `run.log`, while every case
directory also produces `Hydro_ioniz_adv.txt`, `Ion_species_adv.txt` and `OI_levels_adv.txt`. So
`post_process_adv.f90` (1105 lines, including the advection ionization systems, the metal
re-solve, the Brent temperature solve and the validity gates) has no golden: a change anywhere in
it cannot fail `make check`. The `_adv` files are what `EXHALE_transit.py` and
`collisional_validity.py` consume by default.

---

### A failed build in `run_check.sh` does not stop the script, and its output is discarded

Status: CONFIRMED
Severity: P2
Location: `backup/regression/run_check.sh:243-244`

```
243      echo "[build] rebuilding EXHALE.x ..."
244      ( cd "$ROOT" && make >/dev/null 2>&1 ) && echo "  build OK"
```

What the code does: `set -e` is in force, but a command that is not the last of an `&&` list is
exempt from errexit, so a failing `make` is swallowed. Verified on this machine:

```
$ bash -c 'set -e; ( false ) && echo "build OK"; echo "CONTINUED after failed build"'
CONTINUED after failed build
```

The `make -q` guard at line 253 then catches the stale binary and exits 2, so no stale result is
certified. What is lost is the diagnosis: the compiler's error text went to `/dev/null`, and the
message printed says "The dependency graph missed something; run 'make distclean'", which
misdiagnoses a plain compile error.

Two smaller hazards in the same file:

* `check_case` runs under `set +e` (line 285) and does `cp "${REGRESSION_EXE:-$ROOT/EXHALE.x}"
  "$HERE/$c/EXHALE.x"` (line 168) with no status test, then `./EXHALE.x`. A failed copy (a busy
  text file, a permission problem) runs whatever binary the case directory already held.
* the summary line at 293 greps `"$HERE"/*/check.log`, i.e. every case directory in the tree, not
  the cases of this run, so a stale `check.log` from a case not being checked can turn
  "all cases byte-identical" into "within tolerance". The FAIL/PASS exit status (lines 289-291) is
  computed correctly over `$CASES`; only the summary sentence is contaminated. This is the same
  cross-run contamination the directory lock was added for.

`make -q` itself is a sound stale-binary test, `fortdep.py` covers the `use` graph correctly for
this source set (checked: 113 sources, no duplicate basenames, no module name defined twice), and
the flag-hash stamp forces a full rebuild on a compiler or flag change. I found no path by which a
`check` verdict is taken from a binary older than its sources.

---

### `EXHALE_transit.py` decides "metals off" from a column count that is never small enough, using a numpy-version-dependent call

Status: CONFIRMED
Severity: P2
Location: `EXHALE_transit.py:872-887`

```
872 # Plain np.loadtxt: this reads ONE row for its width, not for its content,
...
875 _ncol_ion = np.loadtxt(Ioniz_file, max_rows=1).size
876 do_metals = (_ncol_ion >= 26)
```

What the code does: `write_output` emits the 27 metal-ion columns unconditionally (they are zero
when metals are off, `write_output.f90:159-190`), so a schema-2 `Ion_species.txt` always has at
least 34 columns. Measured across every golden in the tree: 34 for the atomic cases, 38 for the
molecular ones. `do_metals` is therefore always true, the "metals-off run: skipping metal
resonance lines" branch is dead, and a metals-off run instead reads three all-zero columns and
prints the misleading "no OI_levels_adv.txt (run predates the O I level output)".

Numerically harmless (zero columns give flat lines), but the guard is fragile in the other
direction: `max_rows=1` counts comment lines on numpy <= 1.22 and rows on >= 1.23, and the file
starts with three or more `#` lines. On the older behavior `_ncol_ion` would be 0 and every metal
line would be silently skipped. numpy 1.26.4 is what is installed here and it warns on every call:

```
UserWarning: Input line 1 contained no data and will not be counted towards `max_rows=1`.
```

The same script already computes the column count robustly 570 lines earlier
(`EXHALE_transit.py:301-307`, first non-comment line, `len(_s.split())`); the two computations
should be one.

---

### The triaxial branch of `EXHALE_transit.py` carries the Balmer treatment the spherical branch fixed, and mixes two annulus radii

Status: CONFIRMED
Severity: P2
Location: `EXHALE_transit.py:1190-1191`, `:1249-1253`, `:1236-1239`

Three defects, all inside `if geometry == 'triaxial':`, which is dead by default
(`geometry = 'spherical'` is hard-coded at line 1177):

1. `nlo_of_r['Halpha 6563'] = interp1d(r, n2_cm*1e6, ...)` and the matching `tri_lines` entry use
   the multiplet oscillator strength `f_Ha` on the TOTAL n=2 population. The spherical path uses
   the sub-level-resolved `f_Ha_2s*n_2s + f_Ha_2p*n_2p` and says why in its own comment
   (`EXHALE_transit.py:651-653`): the 2s and 2p populations depart strongly from the 1:3
   statistical ratio, so the multiplet average is not applicable. The triaxial path keeps the
   superseded form.
2. `prob = ... / (r_hi**2 - 1.0)` averages the transmission over the annulus `[1, r_hi]` with
   `r_hi = recon.r_sub[-1]`, but the very next line weights it with `A_atm`, which is
   `pi*(Rib*Rp)^2` and `Rib = min(r[-1], R_star/Rp)` (line 380). When the star is small enough for
   the cap to bind, the two annuli differ and the depth is wrong by their area ratio.
3. `np.exp(-np.abs(tau))` hides the sign of a negative optical depth instead of refusing it.

---

### The rotational Doppler shift is applied in the wrong direction

Status: CONFIRMED
Severity: P2
Location: `EXHALE_transit.py:774-775`

```
774			acc += np.interp(l_onde*(1.0 + v/c_light), l_onde, prof,
775			                 left=prof[0], right=prof[-1])
```

Gas receding at line-of-sight velocity `v` absorbs at observed wavelength
`lam_obs = lam_rest*(1 + v/c)`, so the shifted profile is `prof(lam/(1 + v/c))`, i.e. the argument
should carry `(1 - v/c)` to first order. The code shifts the other way. The azimuth samples are
`cos` of a uniform grid over the full circle (line 771), which is symmetric under
`phi -> phi + pi`, so the averaged result is unchanged and no number in the repository moves. It is
reported as a sign error in an equation, not as a wrong output.

The same sign convention appears in the Voigt argument itself (`wofz(X - v_x/v_th + ...)`), where
it also cancels because a spherical radial outflow has a line-of-sight velocity distribution
symmetric about zero along the chord. Neither is a defect of the present spherical model; both
would matter for an asymmetric geometry, which is what the triaxial branch above is for.

---

### `exhale_io` reads `input.inp` where the run resolved something else

Status: CONFIRMED
Severity: P2
Location: `examples/exhale_io.py:346-374` (`read_input`), `:428-441` (`mdot_log10`);
`src/utils/collisional_validity.py:511-512`

`write_resolved_config` (`write_setup_report.f90:672-700`) exists so that consumers read the
radius, mass, temperature and He/H the wind actually used rather than the pre-handoff values of
`input.inp`, and its header says so. `EXHALE_transit.py` does read it
(`exhale_transit_lib.py:355-379`). `exhale_io.read_input` does not: `mdot_log10` scales by
`run.inp['Rp_RJ']**2` taken from `input.inp`, and `collisional_validity.collisional_diagnosis`
builds `Rp_cm` and `Mp_g` from the same place, so every length, every mean free path, every
Knudsen number and the reported `log10 Mdot` are wrong by the square of the radius ratio whenever
`base.inp` or a lower-atmosphere profile moves `r_base`. No shipped case exercises it (both
`mol_base_handoff` and `lower_profile` restate the input radius), so the defect is latent.

Two documentation defects in the same file:

* `read_input`'s docstring promises `LX` among the returned keys; the returned dict has no `LX`
  (only `LEUV`).
* `physical_cell_rows`'s docstring says the row-layout header is carried by
  "Only `Hydro_ioniz(_adv).txt` ... today; the other profile files fall through to rule 2". The
  writer emits it on units 2, 3, 4, 71, 73, 74 and 75, and the shipped goldens confirm it:
  `Ion_species.txt` line 3 and `OI_levels_adv.txt` line 6 both carry
  `# rows 504: 2 ghost cells at each end; physical cells are rows 3 to 502`.

---

### `collisional_validity.py`: the electron-ion energy coupling rate carries an ion-fraction factor its own definition does not have, and the sound speed is not the code's

Status: CONFIRMED
Severity: P2
Location: `src/utils/collisional_validity.py:588-592` and `:290`, `:156-163`

```
588    nuE_ei = np.zeros_like(T)
589    for s in ('HII', 'HeII', 'HeIII', 'H2p', 'H3p', 'HeHp'):
590        if s in dens:
591            nuE_ei += (dens[s] / np.maximum(ne, 1.0e-300)) \
592                * energy_equipartition_frequency(T, ne, dens, 'e', s)
```

`energy_equipartition_frequency(T, ne, dens, 'e', s)` already contains `dens[s]` through
`momentum_transfer_frequency` (line 384: `return dens[t] * kb_erg * TK / (mu_g * Dhat)`), so the
sum is `sum_s (n_s/n_e) nu^E_(e,s)`, quadratic in the ion density. The definition stated in the
module header (line 176) is `nu^E_st = 2 mu_st/(m_s+m_t) nu_st`, whose total over the ion species
is `sum_s nu^E_(e,s)`. The extra factor is 1 for a pure hydrogen plasma and below 1 wherever
helium or metals donate electrons, so `tau_E_ei` is over-estimated and the reported
`tau_E/tau_heat` is conservative. It is still not the quantity the header defines.

```
290 gamma_ad = 1.666666666667   # parameters.f90:385, the code's polytropic index
```
and the header (lines 156-163) claims "This is the code's own sound speed, not a re-definition:
`eval_dt.f90:29` computes cs = sqrt(g*p/rho) with g = 1.666666666667 (`parameters.f90:385`)".
Both statements are now false:

* `eval_dt.f90:31-38` branches on `caloric_mixture_active` and uses
  `adiabatic_index_from_state(j, rho, p)` when it is set. `caloric_state_from_composition`
  (`caloric_eos.f90:186-205`) sets that flag for any run with `Molecular chemistry: True` and H2 in
  at least one cell, which is every molecular case. So for a molecular run the script's sonic point
  is computed with 5/3 where the solver used gamma_eff (about 1.4 to 1.5 in an 80 percent H2
  layer), moving the critical radius the whole verdict turns on.
* the cited line numbers no longer point at those statements: `eval_dt.f90:29` is inside the
  variable declarations, `parameters.f90:385` is in the middle of the `ionization_transport`
  comment, and `gamma_ad` is defined at `parameters.f90:683`.

For an atomic run the script's sound speed is exact. I did NOT re-run the diagnostic on a
molecular case to measure how far the sonic point moves.

---

### Column-integrated ledgers in `write_output` sum the ghost cells

Status: CONFIRMED
Severity: P2
Location: `src/modules/files_IO/write_output.f90:653`, `:727`

```
653      do j = 1-Ng,N+Ng
654         dr_cm = dr_j(j)*R0*opa_pf(j)
```

The FUV photon and energy ledger (`absph`, `absen`, `heat_col`) and the infrared exchange ledger
(`ir_emit`, `ir_abs`) integrate over `1-Ng..N+Ng`, so two ghost cells at each end enter a column
integral over the physical atmosphere. The FUV ledger is internally consistent because the closed
form it is checked against is evaluated at `tau_fuv(1-Ng, ib)`, i.e. at the same innermost ghost;
the infrared ledger is not, since its bound `ir_bound` is built at `r(1)`. Same class as the ghost
issue the row-layout header was introduced for.

---

### `EXHALE_transit.py` output files do not state their wavelength frame, and the frame is mixed

Status: CONFIRMED (mixed frame), documentation defect
Severity: P2
Location: `exhale_transit_lib.py:56-63`, `:92-99`; `EXHALE_transit.py:957-961`, `:989-1015`

The line list mixes conventions: He I 10830.34/10830.25/10829.09, H-alpha 6562.8, H-beta 4861.35,
Ca II 3933.663/3968.469 and Na I 5889.951/5895.924 are AIR wavelengths; Ly-alpha 1215.6701,
Mg II 2796.352/2803.531 and the O I 1302.168/1304.858/1306.029 triplet are VACUUM. Each choice is
defensible on its own (the He fit is explicitly labelled `frame='air'` at
`EXHALE_transit.py:1376`, and the O I record writes "vacuum wavelengths"), but the `tpm_<line>.txt`
headers written at `EXHALE_transit.py:1287-1290` say only `lambda[A]`, so a consumer comparing
against digitized observations has no statement of which frame each curve is on. For the He 10830
triplet the air/vacuum difference is 2.98 A, or 82 km/s.

I verified the atomic data itself against NIST values for all of them and found no error: the He
triplet f-values 0.2996/0.1797/0.0599 are in the 5:3:1 statistical ratio and are attached to the
right components; Mg II 0.608/0.303 with A = 2.60e8/2.57e8; Ca II 0.6267/0.3116 with
1.47e8/1.40e8; Na I 0.641/0.320 with 6.16e7/6.14e7; O I 0.0520/0.0518/0.0519 with
3.41e8/2.03e8/6.76e7 attached to the 3P2/3P1/3P0 lower levels.

I also checked, and found correct: the disk-average normalization
`avg = ((A_star - A_atm) + (A_atm - A_planet)*<T>)/A_star * A_star/(A_star - A_planet)`, which
returns 1 in the continuum and therefore an EXCESS over the opaque disk; the annulus mean
`trapz(2*T*r, r)*A_planet/(A_atm - A_planet) = <T>` for `r_grid` running from 1 to `Rib`; and the
Christie, Arras and Li (2013) 2s/2p system in `n2_populations` (the 2x2 inverse, the detailed
balance of the reverse collisional rates, the statistical weights, and the `alpha_2l*ne*nHII`
cascade source). There is NO limb darkening anywhere: the stellar disk is uniform in brightness.
That is a stated-nowhere approximation of the transit depths.

---

## Checks that came out clean

Recording these so the negative results are on the record.

* **`input_read` defaults against `docs/input_schema.md`.** I compared every default in the schema
  tables (K1 to K42, the core block lines 1 to 22, and the `base.inp` table) against
  `parameters.f90` and the parsing code. Every one agrees: `du_th=1e-3`/`du_th_plm=-1`,
  `flux_spread_th=2.0e-5`/`r_flux=1.2`, `he_kzz=0`, `lowmach_damp_eps=-1`/`mach_th=1e-3`,
  `dr_base=2.0e-4`/`N_low_cells=50`, `N=500`, `CFL=0.6`, `shapiro_eps=-1`/`every=4`,
  `count_max=1e6`, `coronal_cutoff_width=0.1`, `stall_tol=1e-6`/`N_stall=2000`,
  `newton_du_switch=1e-2`, `visc_s=0.7`, `caloric_eos` = ladder, `mol_reaction_heat` = on,
  `h2_double_ionization` = chung80, `h2_neutral_dissociation` = on,
  `carrier_transport` = `thereis_oxychem`, `base_p_ubar=1.0` microbar,
  `p_base_bar=1e-6` bar, `lya_star_boost=5.0`, `dv_star_lya=70.0`.
* **Keys parsed but not listed, or listed but not parsed.** The keyword-loop labels and the
  `known_keys` array agree exactly; the only names in one and not the other are the six `base.inp`
  keys, which have their own list. No key is parsed without being warned about, and no key in
  `known_keys` is unreachable.
* **Case sensitivity.** It is inconsistent but matches the schema everywhere I checked:
  `Load IC`, `Do only PP`, `Force start`, `Transonic IC`, `Jlya escape-prob`, `Newton solver`,
  `Brent solver` and `Use only EUV` accept only the capitalized word, while the later keys accept
  both cases. The schema documents the capitalized-only form for each of them. Worth unifying, but
  not a discrepancy against the documentation. One trap in the assignment form used by
  `Molecular reaction heat`, `H2 neutral dissociation`, `Molecular carrier transport`,
  `Ionization transport` and `Coupled carrier solve`
  (`x = (str .eq. 'True' .or. str .eq. 'true')`): any other word, `TRUE` included, silently sets
  the flag FALSE, and for the first two the default is TRUE, so a typo turns a default-on channel
  off with no message.
* **Unit conversions of the inputs.** `R0*RJ`, `Mp*MJ`, `a_orb*AU`, `Mstar*Msun`, `R_star*Rsun`,
  `n0 = 10**n0`, `v0 = sqrt(kb T0/mu)`, `p0 = n0 mu v0^2 = n0 kb T0`, `q0 = p0 v0/R0`,
  `b0 = G Mp mu/(kb T0 R0)`, `n0 = p_base/(kb T0 ntot_bc)` with 1 bar = 1e6 erg/cm^3 and
  1 microbar = 1 erg/cm^3. All consistent. He/H is by NUMBER throughout, `mass_per_H` in m_H
  carries the mass, and `H_nuclei_in_H2_fraction = 2 q_H2 (1+HeH)/(1 + q_H2)` in
  `write_resolved_config` is the correct inversion of `q_H2 = (x2/2)/(1 + HeH - x2/2)`.
* **Column schemas against the columns written.** `Hydro_ioniz(_adv).txt` header names 7 columns
  and 7 are written. `Ion_species(_adv).txt` names `r` + 6 + `n_mion` (+4 molecular, +3 oxygen) and
  writes the same, unconditionally including the metal block: verified against the goldens (34
  columns for the atomic cases, 38 for the molecular ones). `Cooling_breakdown.txt` names
  cols 1-14 plus `n_mion` and writes the same, and `exhale_io.COOL_GAS_CHANNELS` matches it
  entry for entry. `FUV_bands.txt` builds its column header from the band table rather than
  writing it out.
* **Restart round trip of the fields.** `load_IC` restores, by label, every species column the
  writer emits: H/He stages, the He 2^3S metastable (with the He I column carrying the total, so
  the triplet is not double counted, consistent with `calc_rho`'s `bsp_is_excited_level` rule), the
  27 metal ion stages, the four molecular carriers and the three oxygen carriers. `rho` is
  reconstructed from those densities under the run's own mass policy rather than read back, which
  is the right choice. `v`, `p`, `T` are read and adimensionalized. The He/H rescale, the
  abundance rebuild for an absent element, the handoff renormalization and the oxygen seeding all
  conserve their elements by construction. The `# coupling:` header is parsed field by field and
  unknown fields are ignored. The one defect I found on this path is the grid overwrite, above.
* **Lower-profile reader units and matching.** `p` in bar throughout (`p_match_bar`, `p_top_bar`,
  `p_deep_bar` and the `p` column), `r` in R_J (fed into `R0`, which is later multiplied by `RJ`,
  and read back as `r(jc)*R0/RJ` in `eddy_diffusion_on_grid`), `Kzz` in cm^2/s, `q_H2` a volume
  mixing ratio, `X_<El>` an El/H NUCLEI ratio (fed through the same `set_element_abundance` door
  as `metals.inp`). Interpolation is linear in ln p with no extrapolation, exact at a node, and the
  table is refused unless `p` is strictly decreasing, `p_top_bar < p_match_bar`, `p_match_bar` is
  covered and `solution_id` is present. The `K_zz` interpolation is linear in the profile's own
  radius, which the routine states is the same weight as linear in ln p within a bracket. The one
  gap is that `r` is assumed monotonically increasing with level index and never checked.
* **`make check` staleness.** `make -q` after `make` is a sound guard and it is a hard stop;
  `fortdep.py` regenerates the `use` graph whenever any source changes and make restarts on it; the
  flag-hash stamp renames on any change of compiler path, version, `FFLAGS` or `MODFLAG` and forces
  a full rebuild; there are no duplicate source basenames and no module name defined in two files.

---

Files read completely: `src/EXHALE_main.f90` (3063 lines),
`src/modules/files_IO/input_read.f90` (2589), `src/modules/files_IO/load_IC.f90` (758),
`src/modules/files_IO/write_output.f90` (886),
`src/modules/files_IO/write_setup_report.f90` (848),
`src/modules/files_IO/lower_atmosphere_profile.f90` (554),
`src/modules/post_process/post_process_adv.f90` (1105), `EXHALE_transit.py` (1572),
`exhale_transit_lib.py` (589), `examples/exhale_io.py` (450),
`src/utils/collisional_validity.py` (823), `src/utils/fortdep.py` (79), `Makefile` (400),
`backup/regression/run_check.sh` (310).
Read in part, to close an argument: `src/modules/init/define_grid.f90` (the `j_min`/`j_flux`
block), `src/modules/init/parameters.f90` (the default block and the constants),
`src/modules/functions/utilities.f90` (`calc_rho`, `write_coupling_state_header`),
`src/modules/functions/composition.f90` (`he_ground_singlet_density`, `get_species_densities`),
`src/modules/states/caloric_eos.f90` (`caloric_state_from_composition`,
`adiabatic_index_from_state`, the monatomic branch), `src/modules/time_step/eval_dt.f90`,
`src/modules/radiation/util_ion_eq.f90` (the `Cooling_breakdown.txt` writer),
`docs/input_schema.md` (the whole key table).

Not checked: the hydrodynamic interior of `EXHALE_main.f90` lines 900-1400 beyond the control flow
(reconstruction, RK stages, positivity repair, Riemann fluxes) - another auditor has it; the
`src/modules/wind_ae/` tree; `roche_recon.py`, `he_line_metrics.py`, `EXHALE_plots.py`,
`eta_approx.py`, `run_EXHALE.sh`, `src/utils/element_budget.py`, `src/utils/run_lower.py`,
`src/utils/vulcan_to_base.py`; `metals_input_read.f90` and `opacity_input_read.f90` (separate
parsers, out of the listed group); `backup/regression/compare_within_tolerance.py` and
`roundtrip/roundtrip_check.sh`; the OpenMP correctness of the ionization sweep. I ran no EXHALE
binary and rebuilt nothing in the repository: the only compilation was a standalone extract of
`get_word` in the scratch directory. Every number quoted above was measured from files already in
the working copy, or from that probe.

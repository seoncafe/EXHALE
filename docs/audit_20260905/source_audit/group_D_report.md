# Audit group D report: local ionization/chemistry solvers, composition, species table, element diffusion, global parameters

Repository root: `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` (HEAD 35d9dd5).
Probes built in the scratch directory (`jac_probe.f90`, `tol_probe.f90`, `tp_init.f90`),
compiled with the gfortran on PATH (conda-forge 16.2.0) against the existing `build/*.o`.
Nothing under the repository was modified and no production run, `make` or regression case
was executed.

---

## 1. Two different Jupiter radii are applied to the same `input.inp` line

Status: CONFIRMED
Severity: P1 (opt-in path: IC mode: windae, and the standalone wind_ae_ic.x)
Location: src/modules/init/parameters.f90:691; src/modules/files_IO/input_read.f90:1465;
          src/modules/wind_ae/wae_exhale_input.f90:14, 41, 61

    parameters.f90:691        real*8,parameter ::  RJ      = 6.9911d9         ! Jupiter radius (cm)
    input_read.f90:1465          R0     = R0*RJ
    wae_exhale_input.f90:14      real*8, parameter :: RJ   = 7.1492d9
    wae_exhale_input.f90:41      else if (has(line,'Planet radius'))  then; call aft(line,val); read(val,*) RpRJ
    wae_exhale_input.f90:61      wae_par%Rp        = RpRJ*RJ

What the code does: `wae_read_exhale_input` opens the SAME `input.inp` the main reader opens
and parses the SAME key, `Planet radius [R_J]`. `input_read` multiplies it by the global
`RJ = 6.9911e9` cm (Jupiter's volumetric mean radius); the Wind-AE bridge multiplies it by
its own `RJ = 7.1492e9` cm (the IAU nominal equatorial radius). The two differ by 2.261 per
cent, so one run carries two planetary radii: for `Planet radius: 1.38` EXHALE solves at
`R0 = 9.648e9` cm while the Wind-AE initial condition is built at `Rp = 9.866e9` cm.

What it should do: one input key, one radius. The physical question of which R_J to adopt is
separate (transiting-planet radii are quoted against the equatorial 7.1492e9 cm almost
universally, and `docs/p23_published_profiles.md` already records that EXHALE's constant is
the mean radius); the defect is that the code holds both answers at once and hands them to
two halves of the same calculation.

How it was verified: read both readers end to end and confirmed they parse the same key from
the same file. `src/utils/vulcan_driver.py:135-140` shows the same double definition was
already found and fixed on the Python side ("the upstream literal 7.1492E9 is 2.3% larger"),
so the Fortran bridge is the remaining one.

Reaches: `IC mode: windae` (`ic_mode = 4`) and `make wind_ae_ic` / `wind_ae_ic.x`. The
Wind-AE seed is a warm start EXHALE then relaxes onto its own planet, so a converged EXHALE
answer is not wrong; the Wind-AE oracle solution, its reported Mdot and the seed grid in
`inputdata/windae_grid/` are for a planet 2.3 per cent larger. No default regression case
uses `IC mode: windae`, so no golden moves.

Two more constants in the same block differ from the global ones by smaller amounts:
`MJ = 1.8982d30` against `1.898d30` (0.011 per cent) and `MSUN = 1.98842d33` against
`1.989d33` (0.029 per cent). Same class, below the 0.1 per cent threshold.

---

## 2. The Wind-AE bridge mixes the 4 m_H helium convention with the CODATA helium mass

Status: CONFIRMED
Severity: P1 (opt-in path: IC mode: windae / wind_ae_ic.x)
Location: src/modules/wind_ae/wae_exhale_input.f90:71-77

          ! composition: He/H number ratio -> H/He mass fractions
          wae_par%HX(1) = 1.0d0/(1.0d0 + 4.0d0*exh_HeH)
          wae_par%HX(2) = 4.0d0*exh_HeH/(1.0d0 + 4.0d0*exh_HeH)
          ! CODATA H-atom and He-atom masses; wind-ae originally passed
          ! 1.6733d-24 and 6.6464790722d-24 here.
          wae_par%atomic_mass(1) = 1.67353284d-24
          wae_par%atomic_mass(2) = 6.6464790722d-24

What the code does: the two mass fractions are built with helium weighing exactly 4 hydrogen
masses, while the atomic masses handed to the same solver are the true CODATA values,
m_He/m_H = 6.6464790722/1.67353284 = 3.971519. Wind-AE recovers the number ratio from `HX`
and `atomic_mass`, so the He/H it actually solves at is
`(HX(2)/m_He)/(HX(1)/m_H) = 4/3.971519 * exh_HeH = 1.00717 * exh_HeH`.

What it should do: the conversion factor must be the same mass the solver is given, i.e.
`HX(2) = (m_He/m_H) exh_HeH / (1 + (m_He/m_H) exh_HeH)`. As written the routine states a
helium abundance 0.717 per cent above the one the run asked for.

How it was verified: read; arithmetic from the two literals in the same routine.

Reaches: the same two entry points as finding 1.

---

## 3. Helium's mass in the equation of state is 4 m_H, but the mass unit is the H atom

Status: CONFIRMED
Severity: P1 (default path, all H/He runs; 0.17 per cent on the mass per H nucleus)
Location: src/modules/init/species_table.f90:84-87 and 104-107;
          src/modules/init/parameters.f90:672-675

    species_table.f90:84-87
          ! ... The numeric values
          ! deliberately MIRROR THE LITERALS THE CODE USES TODAY (He mass 4.0,
          ! not 4.0026, H2/H2+ = 2, H3+ = 3, HeH+ = 5 m_H), so a future
          ! metadata-driven rewrite can stay byte-identical.
    species_table.f90:104-107
          real*8,  parameter :: bsp_mass(n_bsp) = &
               [ 1.0d0, 1.0d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0,                  &
                 2.0d0, 2.0d0, 3.0d0, 5.0d0,                                &
                 16.999d0, 17.999d0, 28.010d0 ]
    parameters.f90:672-675
          ! Mass of the hydrogen ATOM (m_p + m_e - 13.6 eV/c^2), CODATA 2018.
          ! This is the mass unit of the density normalization ...
          real*8,parameter ::  mu      = 1.67353284d-24   ! Hydrogen atom mass (g)

What the code does: `bsp_mass` is in units of `mu`, the hydrogen ATOM mass. Helium is given
4.0 of them, i.e. 6.6941e-24 g. The measured helium atom mass is 6.6464731e-24 g, so the
code's helium is 0.717 per cent too heavy. Every mass that flows from `calc_rho`,
`comp_mass_per_H`, `comp_rho_bc` and `calc_mmw` carries it: at the solar-like He/H = 0.0793
the gas mass per H nucleus is `1 + 4*0.0793 = 1.31720` against the correct
`1 + 3.971519*0.0793 = 1.31494`, 0.172 per cent high, which propagates into `rho_bc`, `v0`,
`p0`, `q0` and Mdot.

What it should do: `bsp_mass(HeI..HeTR) = m_He/m_H = 3.971519` (and HeH+ = 4.971519) if the
unit is the H atom, or the whole table must move to atomic mass units u and `mu` with it.
What is not defensible is the two conventions coexisting: the metals in the same table
already carry their true atomic weights (`melem_A`: 12.011, 15.999, 55.845 ...), so a metal
nucleus is weighed correctly while a helium nucleus is not.

Documentation defect in the same place: the comment says the code uses "4.0, not 4.0026".
4.002602 is the helium mass in *u*, and the table's unit is `mu` (= 1.0079 u), so the value
the comment names as the correct one is itself wrong: the correct value in this table's own
unit is 3.9715, and the error is 0.72 per cent, not the 0.065 per cent the comment implies.

How it was verified: read `calc_rho`/`calc_mmw` (utilities.f90:372-439, 539-574),
`comp_mass_per_H` and `comp_rho_bc` (composition.f90:226-238, 402-405); the numbers are hand
arithmetic from CODATA 2018 (m(4He) = 4.002602 u, 1 u = 1.66053907e-24 g).

Reaches: every run. Fixing it changes every golden by ~0.17 per cent in the density
normalization, which is above the 0.1 per cent "counts as identical" rule.

---

## 4. The constrained-equilibrium continuation ignores the transported-proton constraint

Status: CONFIRMED
Severity: P1 (opt-in path: Ionization transport: True)
Location: src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90:807-863

    814     h2_is_fixed               = ieq_cell%x_h2_fixed
    815     oxygen_carriers_are_fixed = ieq_cell%x_ox_fixed
    ...
    859     if (h2_is_fixed) species_fixed(is_H2) = .true.
    860     if (oxygen_carriers_are_fixed .and. thereis_oxychem) then
    861         species_fixed(is_OH)  = .true.
    862         species_fixed(is_H2O) = .true.
    863     endif

What the code does: `ion_cell_state` carries three "this partition is owned elsewhere"
flags -- `x_h2_fixed`, `x_ox_fixed`, `x_hp_fixed` -- and the fraction systems honour all
three (System_HeH_mol_metals.f90:281-294). `set_molecular_network_layout` honours the first
two and never reads `x_hp_fixed`. So when `Ionization transport: True` hands a cell its
transported H+ fraction, the constrained continuation still makes n(H+) an unknown
(`species_fixed(is_HII)` is never set), still solves the H+ balance row as a reaction row
(`row_of_species(is_HII) = 1`; `full_network_residual` includes it because
`species_fixed(is_HII)` is false), and returns the LOCAL photoionization/recombination root.
The seed does honour the constraint (`molecular_limit_stages`, lines 1464-1467) but the
solve then walks away from it.

What it should do: `if (ieq_cell%x_hp_fixed) species_fixed(is_HII) = .true.`, exactly as for
H2 and the oxygen carriers -- one fewer unknown and one fewer row, so the system stays square
by the routine's own argument at lines 796-805.

Consequence, traced to the caller: the candidate is judged by `normalized_reaction_residual`
(ionization_equilibrium.f90:2484), which calls `ion_system_HeH_mol_metals`, which DOES
replace row 1 by `x(1) - x_hp_fix`. So the returned candidate is measured against a row it
was never asked to satisfy, and is rejected whenever the transported fraction differs from
the local root -- which is the entire regime the option exists for (parameters.f90:369-381
quotes 0.129/0.327/0.585 local against 0.057/0.078/0.106 transported). The result is not a
wrong composition but a path that cannot succeed: each such cell spends up to
`max_field_solves = 40` MINPACK `hybrd` calls and is then discarded to class 4.

How it was verified: read `set_molecular_network_layout`, `set_rung_partition`,
`network_balance_rows`, `full_network_residual` and the caller's acceptance block
(ionization_equilibrium.f90:1570-1640) end to end; confirmed by grep that `x_hp_fixed` is
read in the module only at `molecular_limit_stages:1464` and in the dump/probe I/O.

Reaches: `Ionization transport: True` (default off, requires `Molecular carrier transport`).
No regression case sets it, so no golden moves.

---

## 5. `element_census` computes the charge density and never checks it

Status: CONFIRMED
Severity: P2 (an invariant that is advertised, computed and then discarded)
Location: src/modules/functions/element_census.f90:10-12, 249-289, 375-382

    10        ! ... Beside them it computes the positive charge density, which under
    11        ! neutrality is the free electron density.
    ...
    250       real*8, dimension(1-Ng:N+Ng)           :: n_chg, rho_comp
    ...
    262       call element_nuclei_and_charge(rho, f_sp, n_nuc, n_chg, rho_comp)

What the code does: `element_census_verify` and `element_census_reservoir` both declare
`n_chg`, both pass it to `element_nuclei_and_charge`, and neither ever reads it again;
`element_census_take` stores `snap%n_chg` and nothing compares it. A grep over `src/` shows
the only place the charge density is ever tested is the unit test
(src/tests/element_census_tests.f90:193-196, checks E1e/E1f). So charge neutrality is never
enforced, or even reported, on a live state -- although the module header presents it as one
of the two things the module computes.

The mass closure has a weaker version of the same problem: `mclos` (line 272) is accumulated
but never increments `n_bad`, and when nothing else fails it is printed only under
`EXHALE_ELEMENT_ASSERT=2` (line 292), so with `=1` a closing run prints nothing at all and
the closure is invisible.

What it should do: either gate `|n_chg - n_e|/n_e` against the same tolerance (the state
carries `n_e` through `calc_ne`, and the unit test already shows the two agree to 1e-12), or
say in the header that the charge density is computed for the tests only.

How it was verified: read both routines line by line;
`grep -rn "n_chg" src/ | grep -v element_census.f90` returns only the test file.

Reaches: `EXHALE_ELEMENT_ASSERT=1|2` runs only (the census is a no-op otherwise).

Two smaller gaps in the same routine, same severity:
* `element_census_reservoir` exempts helium when `he_diffusion` is on (lines 387-390) but
  applies the fixed `melem_ab` reservoir test to the metals unconditionally -- yet
  `he_metal_diffusion` moves the metals relative to hydrogen on purpose. A run with
  `he_metal_diffusion` and the assert on would report every diffused metal as broken.
* The comment at lines 188-190, "which is why bsp_mass(CO) = 1 + melem_A(C) + melem_A(O) - 1",
  is arithmetically confused: bsp_mass(CO) = 28.010 = 12.011 + 15.999 exactly, and CO carries
  no H nucleus at all (bsp_nH(13) = 0).

---

## 6. The threadprivate default-initializer claim the code designs around is false

Status: CONFIRMED (measured)
Severity: P2 (documentation; the design it justifies is safe either way)
Location: src/modules/nonlinear_system_solver/ion_cell_state.f90:64-67, 84-86;
          src/modules/nonlinear_system_solver/System_HeH_mol.f90:170-174

    ion_cell_state.f90:64-67
            ! Zero unless the run supplies a Lyman-Werner band flux. Assigned by
            ! ioniz_eq for every molecular cell (no default initializer: this
            ! type is threadprivate, where an initializer only reaches the
            ! master thread).
    System_HeH_mol.f90:170-174
        ! The array carries no declaration
        ! initializer on purpose: it is threadprivate, where an initializer
        ! reaches the master thread only, and every entry is written per cell by
        ! set_mol_turnover_rates before any residual reads it.

What the code does: three sites state, as the reason for a design choice, that a default
initializer on a `threadprivate` object reaches only the master thread.

What is true: OpenMP requires each copy of a threadprivate variable that has an initializer
to be initialized once before its first reference, and gfortran implements that. Measured
with the compiler this project builds with (tp_init.f90 in the scratch directory,
gfortran 16.2.0 conda-forge, OMP_NUM_THREADS=4):

    threads observed: 4
     thread   0  cell%a =   1.23450000E+04   scal =   7.77000000E+02
     thread   1  cell%a =   1.23450000E+04   scal =   7.77000000E+02
     thread   2  cell%a =   1.23450000E+04   scal =   7.77000000E+02
     thread   3  cell%a =   1.23450000E+04   scal =   7.77000000E+02

The design (assign every field of `ieq_cell` per cell) is correct, and I verified that
`ioniz_eq` does assign every field a residual reads on every branch it can take
(ionization_equilibrium.f90:1001-1004 for the H-only system, 1102-1269 for the He systems,
including the `.not. thereis_HeITR` zeroing at 1136-1145). The defect is that the stated
reason is wrong, and a maintainer who believes it will not add an initializer where one
would be the right defence.

The same file contradicts itself: `adv_rates%xe_metal` (line 187) DOES carry `= 0.0d0` while
its type is threadprivate (line 191).

---

## 7. `System_HeH_metals`'s header says the ionization sweep is serial; the same file says it is not

Status: CONFIRMED
Severity: P2 (documentation)
Location: src/modules/nonlinear_system_solver/System_HeH_metals.f90:16-19 vs 54-60

    16  ! per cell through the module-level arrays below, set by the driver via
    17  ! set_metal_coeffs before each hybrd1 call. This is safe because the
    18  ! ionization-equilibrium cell loop in ioniz_eq is serial (no OpenMP),
    19  ! and lets the system grow with the number of metals ...
    ...
    54  ! These cell-by-cell coefficients are set (set_metal_coeffs) and read inside the
    55  ! ionization-equilibrium cell sweep, which is now OpenMP-parallel over cells.

The header's safety argument is the opposite of the truth; the array-based design is safe
because of the `threadprivate` at lines 59-60, not because the loop is serial. The same stale
claim appears at T_equation.f90:20-21 ("Safe because the post-process temperature loop is
serial") -- that one is still true, `post_process_adv` is serial.

---

## 8. `pi` is truncated to 11 significant digits, against the module's own stated policy

Status: CONFIRMED
Severity: P2 (hygiene; 1.3e-11 relative, no physical consequence)
Location: src/modules/init/parameters.f90:665

    665      real*8,parameter ::  pi      = 3.1415926536d0   ! pi

pi is 3.14159265358979324; the literal is 4.1e-11 low (1.3e-11 relative). The same block
argues, for `gamma_ad` twelve lines later, that a constant should be written so "the value is
5/3 to full double precision rather than to the digits a literal happens to carry"
(lines 676-682). pi is the one constant in the block that does not follow that rule, and it
is the only one with an exact expression available (`4*atan(1.0d0)`). It enters
`Mdot = 4 pi r^2 rho v`, the Coulomb logarithm and the ambipolar Debye length. Every other
constant in the block was checked against CODATA 2018 / IAU 2015 and is correct to the digits
given.

---

## 9. The H I / He I collisional-ionization energies are duplicated as single-precision literals in a second file

Status: CONFIRMED
Severity: P2 (extends the known "literal H I/He I ionization potentials in eval_cool" item to
          a file that item does not name)
Location: src/modules/nonlinear_system_solver/T_equation.f90:118-120;
          src/modules/radiation/util_ion_eq.f90:1515-1516

    T_equation.f90:118-120
       coio =  2.179e-11*ci_rate_HI(TT)*nhi                    & ! HI
               + 3.940e-11*ci_rate_HeI(TT)*nhei                & ! HeI
              + e_th_HeII_erg*ci_rate_HeII(TT)*nheii             ! HeII
    util_ion_eq.f90:1515-1516
        coio(j_lo:j_hi) =  2.179e-11*a_ion_HI(j_lo:j_hi)*nhi(j_lo:j_hi)      & ! HI
              + 3.940e-11*a_ion_HeI(j_lo:j_hi)*nhei(j_lo:j_hi)               & ! HeI

The known item names `eval_cool`; the same two literals also sit in `T_equation`, so the
post-process temperature root and the marching cooling carry two independent copies of the
same physical constant. The 2026-09-05 unification (parameters.f90:726-739) made
`e_th_HeII_erg` a single definition and left the other two untouched -- visible on the very
next line of both files, where the He II term uses the named constant and the H I / He I terms
do not.

Values: I(H I) = 13.598434 eV = 2.178711e-11 erg, so 2.179e-11 is 0.013 per cent high;
I(He I) = 24.587389 eV = 3.939338e-11 erg, so 3.940e-11 is 0.017 per cent high. Both literals
are also DEFAULT REAL (`e-11`, not `d-11`), adding another ~1e-7 relative. Under the standing
sub-0.1-percent rule these move nothing; the duplication is the defect.

---

## 10. `parameters.f90` mis-states the `sigma_tab` column order (two ions missing)

Status: CONFIRMED
Severity: P2 (documentation)
Location: src/modules/init/parameters.f90:1062-1066

          ! Metal photoionization cross sections, one column per photo-ionizable
          ! metal ion in species_table iphot order (1=CI,2=CII,3=OI,4=OII,
          ! 5=NI,6=NII,7=MgI,8=MgII,9=SiI,10=SiII,11=CaI,12=CaII,13=NaI,14=KI,
          ! 15=SI). Replaces the former loose s_ci..s_mgii.

`n_mphot = 17` (species_table.f90:57) and `mion_iphot` assigns 16 = Fe I and 17 = Fe II
(species_table.f90:212-215). The array is allocated `sigma_tab(Nl, n_mphot)` with 17 columns
(set_energy_vectors.f90:248). The comment lists only 15, so a reader indexing `sigma_tab` from
it will not find the iron columns.

---

## 11. `solve_ieq` passes one tolerance to two solvers that give it different meanings

Status: CONFIRMED (measured)
Severity: P2 (cost and reporting, not a wrong answer)
Location: src/modules/nonlinear_system_solver/newton_solver.f90:36-79, 84-137;
          src/modules/radiation/ionization_equilibrium.f90:613

    ionization_equilibrium.f90:613     tol = sqrt(dpmpar(1))
    newton_solver.f90:73                  call hybrd1(fcn, n, x, fvec, tol, info_m, wa, lwa, params)
    newton_solver.f90:131         if (dxmax .lt. 1.0d-11*xscale .or. fnorm .lt. ftol) then

What the code does: MINPACK's `tol` is xtol, "the relative error between x and the solution"
(hybrd1.f90:53-55). `newton_dense` receives the same number as `ftol` and compares it against
`sqrt(sum(fvec*fvec))`, the residual 2-norm in cm^-3 s^-1, which scales as n_H^2 times a rate
coefficient. So the two solvers reached through one call are held to two incommensurable
criteria, and the residual branch of the Newton test is density-dependent.

Measured (tol_probe.f90, H/He system, He/H = 0.0793, T = 6000 K, the same rate set at three
densities; tol = 1.4901e-08):

     n_H= 1.00E+06 newton= T conv= 1  ||fvec||_2 at the returned root = 5.2283E-14
     n_H= 1.00E+10 newton= T conv= 1  ||fvec||_2 at the returned root = 8.9511E-10
     n_H= 1.00E+14 newton= F conv= 1  ||fvec||_2 at the returned root = 2.1772E-07

At the base density the attainable residual floor (2.2e-7) is already above `ftol`, so the
`fnorm .lt. ftol` branch is unreachable there and convergence rests entirely on the step test
`dxmax .lt. 1.0d-11*xscale`; a cell that cannot reach 1e-11 in 60 iterations falls back to a
full `hybrd1` solve on top of the Newton work it just discarded.

Second, smaller point in the same routine: when the Armijo loop exhausts `maxls = 20` without
an accepted step, `lam` has been halved once more than the step that produced `xnew`
(lines 114-121), so the convergence measure `dxmax = maxval(abs(lam*dx))` at line 129
under-reports the accepted step by a factor 2 in that branch.

Reaches: every ionization solve of a metals-on or plain H/He run (`use_newton_ieq` default
true). It changes cost and the meaning of the reported `converged` flag, not the accepted
composition -- `ioniz_eq` re-judges every candidate with `ionization_fractions_physical` and
`normalized_reaction_residual` and does not trust the solver status
(ionization_equilibrium.f90:1441-1446).

---

## 12. `calc_mmw` counts nuclei as particles, so it is wrong in a molecular gas

Status: CONFIRMED (no wrong answer today: its only caller is the atomic post-process)
Severity: P2
Location: src/modules/functions/utilities.f90:539-574

    556     if (thereis_He) then
    557         mass_l  = bsp_mass(1)*nh + bsp_mass(3)*nhe
    558         npart_l = nh + nhe + ne

`nh` and `nhe` are the ELEMENT (nucleus) densities `hydrogen_helium_nuclei_density` returns,
which count two H per H2 and three per H3+. `npart_l` therefore counts an H2 molecule as two
particles and an H3+ as three, and `mass_l` omits every molecular carrier's own mass. In a
fully molecular layer the mean molecular weight comes out low by nearly a factor two. The
routine has no `thereis_mol` guard and no comment saying it is atomic-only; its single caller
is post_process_adv.f90:866, which the known-item list already records as molecule-free, so
nothing is wrong in a current run. This is a trap for the next caller.

---

## 13. The trace-metal element-diffusion path is never exercised by the test program

Status: CONFIRMED
Severity: P2 (test coverage)
Location: src/tests/diffusion_tests.f90:130-134;
          src/modules/functions/binary_element_diffusion.f90:550-648, 1943-2055

    130       he_kzz        = 0.0d0
    131       call eddy_diffusion_on_grid
    132       he_ambipolar  = .false.
    133       he_alphaT     = 0.0d0
    134       he_metal_diffusion = .false.

`he_metal_diffusion` is set false in `setup_column` and never set true anywhere in the file
(grep over the whole test program returns only this line), so the whole trace-metal branch of
`element_diffusion_step` -- `solve_trace_element_in_hydrogen`, the metal `GcoX` assembly with
its own `-dln(psi)/dr` term, the exhausted-cell re-seed at lines 639-645 and the
`trace_ratio_under_zero` measurement -- has no test. T13 ("the metals come back with the
hydrogen") exercises `project_elements`, not the trace transport solve. `he_alphaT` is
likewise never set nonzero, so the thermal-diffusion term of `G` (lines 1480-1485 and 617-622)
is untested in both the helium and the metal equation.

What the test program DOES cover, for the record: T1a/T1b/T3/T4/T5/T6/T7/T9/T10/T11/T12/T13/
T14 -- diffusive equilibrium, the Dirichlet boundary flux, the trace limit, elemental
conservation in a closed column, the two-component mass closure, uniform-X preservation, the
molecular/homopause closure with a `he_kzz` sweep, the three ambipolar limits, grid and
time-step convergence, the divergence-point cell, the three friction limits, the pure-helium
band and the strong-settling bound.

`element_census_tests` covers E1 (stoichiometry of `element_nuclei_and_charge` plus the charge
= `calc_ne` identity), E2 (elemental closure through one carrier write-back, including the CO
collapse), E3 (the H2 photoevent ledger) and E4 (the base H2 fraction, both branches). It does
NOT cover `element_census_verify` or `element_census_reservoir` themselves, the `rho_comp`
mass closure, `he_ground_singlet_density`, `element_ratio_HeH`, or any of the composition
scalars `comp_mass_per_H` / `comp_ntot_bc` / `comp_rho_bc` against `calc_rho` / `calc_ntot`.

A related small inconsistency in the operator itself (binary_element_diffusion.f90:623): the
helium equation anchors its base at the INPUT reservoir,
`X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)` (line 495), while the metal equation anchors at
whatever cell 1 currently holds, `fXbase = nXold(1)/nH_phys(1)`. Self-consistent as long as
cell 1 starts at `melem_ab` and is never touched, but it is a different kind of boundary
condition for the same operator.

---

## Verified and found clean (negative results worth recording)

* Analytic Jacobians vs a central finite difference of the same residual. Measured with
  `jac_probe.f90` (built against `build/*.o`), four configurations of a synthetic cell at
  T = 6000 K, n_H = 1e8, He/H = 0.0793, with one metal element (Mg) deliberately absent so the
  pinned identity rows are exercised, and Na/K/S exercising the two-stage branch. Worst
  row-scaled |J_analytic - J_fd|:

       case  1  metals= F cx= F  n=  3  worst = 1.259E-10  at row 1 col 1
       case  2  metals= F cx= T  n=  3  worst = 1.406E-10  at row 1 col 1
       case  3  metals= T cx= F  n= 23  worst = 4.113E-11  at row 23 col 1
       case  4  metals= T cx= T  n= 23  worst = 4.113E-11  at row 23 col 1
       (h = 1e-6, so ~1e-10 is the difference truncation floor)

  `jac_system_HeH` and `jac_system_HeH_metals` are therefore correct, including the rank-1
  electron coupling, `heh_jac_local`, the metal blocks, the identity rows and the
  `cx_add_to_jac` / `he_h_cx_jac` charge-exchange contributions. `jac_system_H` is a 1x1
  expression and was checked by hand. `System_HeH_TR`, `System_HeH_TR_metals`,
  `System_HeH_mol` and `System_HeH_mol_metals` carry no analytic Jacobian by design and were
  not probed.

* The three pair diffusion coefficients of `binary_element_diffusion` were re-derived from the
  Chapman-Enskog D = 3kT/(16 n mu Omega^(1,1)) with the stated cross sections. Hard sphere
  reproduces Banks & Kockarts; polarization gives Omega = 3C/16 with
  C = 2.21 pi e sqrt(alpha/mu) and hence D = kT/(2.21 pi e n sqrt(alpha mu)), matching line
  1070; Coulomb gives the prefactor 3/(4 sqrt(2 pi mu)) = 0.2992/sqrt(mu), matching line 1094.
  Blanc's law across carriers and stages, the Peclet hybrid, the donor/acceptor choice in
  `element_face_flux`, the M-matrix sign pattern, the tridiagonal Thomas solve, the sign of the
  trace-metal flux and the `idom` pair packing/unpacking ((is-1)*n_hcar + it against
  is = (idom-1)/n_hcar + 1) all check out.

* The ambipolar field. eE = -kT dln(n_e T)/dr and
  dmeff = (m_He - m_c1) - (Zbar_He - Zbar_c1) eE/(m_H g) reproduce the three stated limits by
  hand: 3 (neutral), 2.5 (H+ plasma, eE = m_H g/2) and 5/3 (He++ plasma, eE = 4/3 m_H g). The
  sign of `dlnpsi` is right: the variable holds +dln(psi)/dr and G subtracts it, matching the
  header formula.

* The molecular network rows (`mol_heh_rows`) conserve H nuclei, He nuclei and charge channel
  by channel: the H-nucleus closure 1 - x1 - x4 - x5 - x6 - x7 matches
  n_H = n_HI + n_HII + 2n_H2 + 2n_H2+ + 3n_H3+ + n_HeH+, the He closure includes the HeH+
  nucleus, R11 puts one H2 and one H2+ out of one H3+ and one H (4 nuclei in, 4 out), the
  double-ionization channel carries its factor 2 in row 1 and not in row 4, and the
  Penning/associative split of Q31 and k_ion_H2 is charged once each to rows 1, 5, 7 and 8.
  `hydrogen_helium_row_terms` in `constrained_chemical_equilibrium.f90` mirrors rows 1 and 2
  term for term. The advection-correction systems (`System_implicit_adv_*`) carry the same
  Penning and charge-exchange bookkeeping with the correct normalizations by nucleus count
  (n_he = heh_loc*n_h).

* He 2^3S is counted exactly once as a level of He I in every budget I could reach: `calc_ne`,
  `calc_ntot`, `calc_rho`, `hydrogen_helium_nuclei_density`, `element_ratio_HeH`,
  `element_nucleus_counts`, `element_nuclei_and_charge` and `carrier_counts` all skip
  `bsp_is_excited_level` or never receive the triplet column, and `he_ground_singlet_density`
  is the single place the difference is taken.

* The composition conversions in `composition.f90`: x2 = 2q(1+He/H)/(1+q),
  q_max = 0.5/(0.5 + He/H) and h2_bound_fraction = 0.5 x2/(1+He/H) were re-derived from the
  definitions and are correct.

* Physical constants in `parameters.f90` were checked against CODATA 2018 / IAU 2015:
  `kb_erg`, `kb_eV`, `hp_erg`, `hp_eV`, `c_light`, `Gc`, `erg2eV`, `mu`, `Rsun`, `MJ`,
  `M_earth`, `AU`, `parsec` all agree to the digits given. `gamma_ad` is the exact rational. A
  grep for duplicated literals of k_B, m_H, h, c, G and the eV across `src/` found every copy
  numerically identical to the global one; the only value-level disagreements are the
  R_J / M_J / M_sun set of finding 1. (Cool_coeff.f90:3445, 3471 use the truncated 8.61733e-5
  for k_B in eV/K, but only inside the legacy Abel et al. 1997 fit, whose coefficients were
  published against the authors' own conversion, so keeping the literal is the right call
  there.)

* `mion_ethr` was checked ion by ion against the adopted NIST ionization potentials (C I 11.26,
  C II 24.38, O I 13.62, O II 35.12, N I 14.53, N II 29.60, Mg I 7.646, Mg II 15.035,
  Si I 8.152, Si II 16.35, Ca I 6.113, Ca II 11.87, Na I 5.139, K I 4.341, S I 10.36,
  Fe I 7.902, Fe II 16.199): all correct. `melem_Z`, `mion_stage`, `mion_z2`, `mion_iphot`,
  `melem_i0` and `melem_top` are internally consistent, and the extraction at
  ionization_equilibrium.f90:1867-1879 correctly never writes an i0+2 column for a two-stage
  element (which would land on the next element's neutral stage).

* The `hybrd1` retry ladder in `ioniz_eq` does NOT silently accept a failed status: the MINPACK
  exit code is used only to RANK candidates (ok_rank = 2 versus 1), and acceptance requires
  `ionization_fractions_physical` and `normalized_reaction_residual <= 1e-6` independently; a
  non-finite residual is forced to class 4 and counted. Attempt ordering, the "strict
  improvement only" tie-break, the restoration of `x_root_best` and the `cx_metal_base`
  set/reset around every merged solve are all correct.

* `T_equation` / `solve_T_brent`: the upward log-spaced scan from max(0.05 x_guess, 1/T0) takes
  the FIRST sign change, which does select the lowest root and structurally avoids the spurious
  hot root; TT = max(TT, 1.0d0) floors only the arguments of the cooling fits, not the energy
  balance, and the molecular branch is the legacy expression divided by (gamma - 1) with
  u = T/(gamma-1) replaced by the caloric energy, so x_h2 = 0 reproduces the legacy row exactly
  (the branch test is `x_h2 .gt. 0`, so an atomic cell is bit-identical).

* `fdjac1`'s `unit_step_floor`: the interface carries it, both callers pass it explicitly
  (`hybrd1` false, the constrained solve true), the dense and banded branches share one
  `difference_step`, and hybrd1's ml = mu = n-1 makes the dense branch the only reachable one.
  No interface misuse. Every residual declares `params` no larger than the `params(60)` its
  callers actually pass (`ionization_equilibrium.f90:501`, `post_process_adv.f90:205`), and
  none of them reads it.

* OpenMP: `met_*`, `cx_kc`, `mk*`, `ok*`, `oj*`, `mol_inv_turnover`, `ieq_cell`, `adv_cell`,
  `teq_cell`, `sys_x`, `sys_sol`, `wa`, `info` and the constrained-equilibrium layout state are
  all threadprivate; `cx_metal_base` is `copyin`; the `nt_calls` / `nt_fallback` counters use
  `!$omp atomic`; the once-only dump guard is deliberately shared and taken inside
  `!$omp critical`. `equilibrium_from_molecular_limit` resets every piece of the state of each rung at
  entry (lines 546-553) so a cell's answer cannot depend on which cell the thread solved before
  it. I found no shared-write race in this group.

* Name shadowing: the global `mu` (H atom mass) is shadowed by locals named `mu` only in
  `wae_eqns.f90`, `wae_soe.f90` and `wae_continuation.f90`, none of which `use
  global_parameters` -- checked by grep. `binary_element_diffusion` states and honours the
  `tscale`-not-`t0` rule. `newton_solver` imports `global_parameters` with an `only:` clause so
  its local `info` does not touch the global one. No new shadowing found.

---

Files read completely: src/modules/nonlinear_system_solver/System_H.f90, System_HeH.f90,
System_HeH_TR.f90, System_HeH_metals.f90, System_HeH_TR_metals.f90, System_HeH_mol.f90,
System_HeH_mol_metals.f90, System_implicit_adv_H.f90, System_implicit_adv_HeH.f90,
System_implicit_adv_HeH_TR.f90, ion_cell_state.f90, ion_residual_core.f90, newton_solver.f90,
T_equation.f90, constrained_chemical_equilibrium.f90, constrained_equilibrium_probe.f90,
hybrd1.f90, fdjac1.f90; src/modules/functions/composition.f90, element_census.f90,
utilities.f90, binary_element_diffusion.f90; src/modules/init/parameters.f90,
species_table.f90; src/tests/element_census_tests.f90. src/tests/diffusion_tests.f90 was read
by its test-by-test structure and by grep for the switches each test sets, not line by line.

Read as far as the callers required, not audited as their own group:
src/modules/radiation/ionization_equilibrium.f90 (the cell sweep, the acceptance ladder,
`normalized_reaction_residual`), src/modules/radiation/charge_exchange.f90 (the He/H pair,
`cx_add_to_fvec/jac/turnover`, `cx_dens_lin`, `cx_fvidx`),
src/modules/wind_ae/wae_exhale_input.f90, src/modules/files_IO/lower_atmosphere_profile.f90
(`eddy_diffusion_on_grid` only), src/modules/files_IO/input_read.f90 (line 1465 only).

Not checked: the MINPACK internals (hybrd, dogleg, qrfac, qform, r1mpyq, r1updt, enorm, dpmpar)
beyond the caller interfaces; the rate-coefficient VALUES in mol_rates, oxygen_rates,
Cool_coeff and the cx_rate table against their published sources (group C/E territory, and the
brief's known list already covers several of them); h2_photo_channels; the accuracy of
q_h2_equilibrium; whether the O I column written to f_sp is consistently the FREE atomic oxygen
everywhere downstream of ioniz_eq; diffusive_photochemistry (another group); the
diagnostic-only bodies of cce_probe_from_dump, write_continuation_dump and
write_element_flux_profile beyond a read for correctness of the quantities they report. I ran
no EXHALE case and no regression matrix, so no finding here carries a measured effect on a
production solution.

# Audit report, group C: lower atmosphere and molecular physics

Repository `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` (HEAD 35d9dd5).
READ-ONLY audit; nothing under the repository was modified. Probes were compiled
and run only in the scratch directory.

## Summary of what checked out and needs no finding

Verified against the published sources and by probe, all correct:

* **Koskinen et al. (2022) Table 1, R1-R23** transcribed into `mol_rates.f90`,
  compared line by line against `pdftotext -layout` of
  `references/Koskinen_2022_ApJ_929_52.pdf` p. 19. Coefficients, exponents,
  the `(300/T)` form, the two three-body entries returned as two-body
  equivalents with an explicit `n` factor, and the temperature each fit is
  given in all match. The two deliberate departures (R12 by detailed balance,
  R15 on Cohen and Westberg instead of Ham) are documented at the code site.
  R14's exponential uses `Te` where the paper prints `T`; EXHALE carries one
  temperature so the two readings coincide numerically, and the header says so.
* **The H + H <-> H2 equilibrium constant.** `K_eq = Lam(2m) q_int exp(D0/kT) /
  [4 Lam(m)]^2` is the correct concentration equilibrium constant: the nuclear
  spin weights 3:1 inside `q_rovib_H2` cancel the 2x2 per atom on the
  denominator and the residual electronic factor 1/4 is the standard one.
  Probe (`auditC/probe_r12.f90`, linked against the production `mol_rates.f90`)
  reproduces every verification number the R12 comment quotes:
  `k_diss(detailed balance)/k_diss(Baulch 1992)` = 0.112, 0.356, 0.613, 0.823,
  0.974, 1.178, 1.052 at 1000, 1500, 2000, 2500, 3000, 5000, 8000 K against
  the comment's 0.11, 0.36, 0.61, 0.82, 0.97, 1.18, 1.05; `q_rovib(2000 K)` =
  50.2245 against 50.22; D0 = 4.478075 eV. I also re-derived the Visscher
  cross-check by hand: `K_eq(Visscher dilute)/K_eq(code)` = 0.9927 at 2000 K
  and 1.0816 at 3000 K against the comment's 0.993 and 1.082.
* **`lower_column.f90` against its own validation gate.** Probe
  (`auditC/probe_lc.f90`, production module) with Koskinen Model A
  (M = 8.6813e28 g, R_1bar = 2.5559e9 cm, T_eq = 1140 K, He/H = 0.0793,
  1 bar to 1 microbar) gives r_base/R_1bar = 1.34311 and q_H2 = 0.83836,
  q_H = 0.02657, against the paper's 1.34 and ~0.84. The inversion
  `x2 = 2 q (1+fhe)/(1+q)` and the composition it builds are algebraically
  identical to Koskinen eq. (12), which I derived independently.
* **H2 infrared line data** in `molecular_infrared_data.f90` (parsed all 1833
  lines and 302 levels with `auditC/parse_h2.py`). Spot checks: 0-0 S(1) has
  dE/k = 844.608 K, A = 4.761e-10 s^-1, T_u = 1015.1 K, g_u = 21; 0-0 S(0)
  A = 2.943e-11; 1-0 S(1) dE/k = 6780.8 K, A = 3.47e-7, T_u = 6951.3 K.
  Level energies (0,1) = 170.492, (0,2) = 509.9, (0,3) = 1015.1, (0,4) =
  1681.6 K and weights 1, 9, 5, 21, 9 are the standard values.
* **The infrared net-exchange algebra.** I derived the two-level net absorption
  from the Einstein relations and it reproduces
  `absr = f_u A dE k_B nbar (exp(dE/T)-1)` exactly, including the stimulated
  term, so the T = T_rad, W = 1 fixed point holds for the line sum. Simpson in
  ln(lambda) with 25 nodes (24 intervals, correct even count) and the
  `B_nu * nu` Jacobian is right.
* **Self-shielding table normalization.** Extracted the three cubes from
  `h2_self_shielding_table.f90` and reproduced the header's own checks:
  sigma_pump at the bottom of the column axis is 2.589e-17 (700 K) to
  2.877e-17 cm^2 (3200 K) against the quoted 2.59e-17 to 2.88e-17; p_single
  0.1427 to 0.1592 against 0.1425 to 0.1592; p_single(1300 K) = 0.1476 against
  0.1475. Trilinear interpolation index order matches the `(n_col, n_temp,
  n_dens)` reshape; the bracket routine clamps at both ends.
* **Band photon energies** in `water_photolysis.f90`: 2hc/(l1+l2) reproduces
  1.71912e-11, 1.63403e-11, 1.48187e-11, 1.05803e-11 erg for B1-B4 and
  1.96483e-11 for LW from the stated band edges. The cell-mean photolysis rate
  `s N tr exp(-tau_out) (1-exp(-dtau))/dtau` shares the absorbed photons
  between H2O and OH exactly in proportion to their contributions to dtau; I
  verified the identity.
* **Carrier operator conservation.** The write-back reproduces `nH_free`,
  `nO_free`, `nC_free` cell by cell provided the limiter has run, which it has
  (`limit_to_element_budget` is called at the end of `solve_carriers`, over
  `j = 1..N+Ng`, the same range the write-back covers); `bsp_mass` of every
  carrier is exactly the sum of its nuclei masses (16.999 = 1 + 15.999,
  17.999 = 2 + 15.999, 28.010 = 12.011 + 15.999), so mass is conserved with
  the elements. The face-flux derivatives `dJl`/`dJr` match the residual's
  flux term analytically, the settling coefficient `G = (m_i - mbar) g/(kT)`
  has the right sign for `Dphi > 0`, and the advective term is the correct
  non-conservative form for the fraction per unit mass paired with the
  hydrodynamic mass row.
* **OpenMP in the carrier solve.** Both parallel loops write only index `j` of
  shared arrays; `ieq_cell`, the `mk*`/`ok*`/`oj*` rate coefficients and
  `mol_inv_turnover` are all in `!$omp threadprivate` lists with the `save`
  attribute; `rnorm` is a `reduction(max:)` and the argmax is found by a
  serial scan afterwards. No race found.
* **No shadowing of `global_parameters` names.** I extracted 220 declared names
  from `parameters.f90` and matched them case-insensitively against every local
  and dummy declared in the three lower-atmosphere modules that `use
  global_parameters`. The only hit is `n` in `locate_log`
  (`molecular_infrared_cooling.f90:273`), and that module imports only
  `kb_erg, hp_erg, c_light`, so `N` is not in scope. `mol_rates` locals are
  clean against its six imports. `lower_column`'s result variable `mu` does not
  collide because that module imports nothing.
* Argument order of `mol_heh_rows` matches the carrier's call site, and row 4
  is production minus loss in cm^-3 s^-1, which is what `carrier_source` and
  the residual assume.

---

## Findings

### The Lyman-alpha photolysis band carries the stellar flux with no atomic-hydrogen attenuation
Status: CONFIRMED
Severity: P1 (wrong answer in an opt-in path: `Oxygen chemistry: True` with a
non-zero `Stellar Lya flux`)
Location:
`src/modules/radiation/util_ion_eq.f90:179-196` and `:355-374`;
`src/modules/lower_atmosphere/water_photolysis.f90:251-255` and `:403-411`
```
	case (ib_B2)
		F = F_Lya_star
	...
	F = dayside_dilution()*F
```
```
				if (ib .ne. ib_LW) tr_out = 1.0d0
				dtau = max(tau_b(j,ib) - tau_out, 0.0d0)
				j_h2o(j,ib) = water_photolysis_rate(fuv_band_flux(ib),    &
				                            ib, tau_out, dtau, tr_out)
```
```
      double precision function fuv_band_optical_depth(ib, N_H2O, N_OH)   &
                                result(tau)
      tau = sigma_H2O_band(ib)*max(N_H2O, 0.0d0)                          &
          + sigma_OH_band(ib) *max(N_OH,  0.0d0)
```
What the code does: band B2 is the H I Lyman-alpha line (1215.67 A, the cross
section is evaluated at line centre and `e_photon_flat_band(ib_B2)` is the
line's photon energy). Its incident flux is `F_Lya_star`, the stellar flux at
the planet's orbit, times the run-wide dayside factor and nothing else. The
only attenuation applied to it on the way down is
`tau_b = sigma_H2O N(H2O) + sigma_OH N(OH)`, and the H2 line transmission
`tr_lines` is forced to 1 for every band but LW. So the H2O and OH photolysis
rates in the Lyman-alpha band are computed as if the atomic hydrogen between
the star and the molecular base were transparent at line centre.

What it should do: Lyman-alpha is a resonance line of the dominant absorber.
The star-ward H I column above a molecular base is 1e19 to 1e21 cm^-2, giving a
line-centre optical depth of order 1e6 to 1e8. The code already computes that
attenuation for its own purposes: `lya_rt.f90:213-226` builds the penetrating
stellar beam as `Jstar = xi F_Lya_star/(4 pi Dnu_star) T_star` with
`T_star = erfc(x1/(sqrt(2) X_s))`, `x1 = sqrt(a tau/sqrt(pi))` from the H I
Lyman-alpha optical depth of the same cell, and stores the result in
`Jlya_arr`. The photolysis band ignores it. The module's own header states the
opposite of the physics here: "Atomic H does not: the bands lie longward of the
912 A Lyman edge" (`water_photolysis.f90:251-252`) is true of the H I
photoionization continuum and false of the band that *is* the H I resonance
line. The comment then treats only the higher Lyman-series lines inside the LW
band as the H-line question, and Lyman-alpha itself is never raised.

How it was verified: read `fuv_band_flux`, `fuv_band_optical_depth`,
`fuv_lw_photon_field` and `water_photolysis_rate` end to end; grepped every use
of `F_Lya_star` in the tree (only `lya_rt.f90:224` and `util_ion_eq.f90:187`);
confirmed `tau_b` has no H I term and `tr_out` is 1 for B2. Size, from the
configuration in `examples/18_oxygen_chemistry/input.inp`
(LW 343.0, B1 137.9, B2 4809.2, B3 694.8, B4 658782.2 erg cm^-2 s^-1): the
unattenuated H2O photolysis rate of each band is 1.66e-4, 5.35e-5, 4.5e-3,
2.07e-4 and 7.7e-2 s^-1, so B2 is the second largest channel there, about 6 per
cent of the unattenuated total. For a star whose FUV is Lyman-alpha dominated
(any M dwarf, and most active G/K stars) B2 becomes the largest channel and the
error is unbounded, because the neglected factor is `T_star`, which the code
elsewhere lets fall to zero.

Reaches: every run with `Oxygen chemistry: True` and a non-zero
`Stellar Lya flux`, through `j_h2o`/`j_oh` -> `set_oxygen_coeffs` -> `oj3, oj4,
oj5, oj7` -> the H2O/OH/O balance, and through
`heat_per_water_dissociation` into the energy equation. No regression case
moves: the default matrix contains no oxygen-chemistry case at all (see the
test-coverage finding below); the affected configuration is
`examples/18_oxygen_chemistry/`.

---

### `equilibrium_constant_conc` overflows to infinity and the oxygen equilibrium seed becomes NaN below about 10 K
Status: CONFIRMED (by probe)
Severity: P2 (numerical hygiene; the trigger temperature is far below any
converged state I could find, so I could not close the reachability argument)
Location: `src/modules/lower_atmosphere/oxygen_rates.f90:770-791` and
`:929-959`
```
      Teq   = t_shomate(ith_H2, T)
      n_std = p_std_cgs/(kb_erg*Teq)
      Kc    = exp(-dG/(R_gas_J*Teq))*n_std**dn
```
```
      if (kc1 .le. 0.0d0 .or. kc2 .le. 0.0d0) return
      ...
      lw(3) = lw(2) + log(kc1) + log(ratio)                ! H2O
      lwmax = max(lw(1), max(lw(2), lw(3)))
      ...
      	w(i) = exp(max(lw(i) - lwmax, -700.0d0))
```
What the code does: for the O1 pair `OH + H2 <-> H2O + H`, dG(298) is about
-63 kJ/mol and stays near that value down to the 10 K floor `t_shomate` clamps
to. `exp(-dG/(R Teq))` therefore exceeds `exp(709)` and returns `+Infinity` for
any gas temperature at or below about 10.4 K. The guard on the next routine
tests `kc1 .le. 0`, which `+Infinity` passes; `log(Inf)` makes `lw(3) = +Inf`,
`lwmax = +Inf`, and `max(Inf - Inf, -700) = max(NaN, -700)`, which gfortran
evaluates to NaN. Both returned fractions are then NaN.

What it should do: the routine's own comment says it is "Evaluated through
logarithms because K_c(O1) K_c(O2) (n_H2/n_H)^2 overflows a double at the cool,
molecular end", so the overflow of the *product* was anticipated but not the
overflow of `K_c` itself. The equilibrium constant should be carried as
`ln K_c` (or the exponent clamped) so that the cold limit returns f_H2O -> 1
rather than NaN.

How it was verified: probe `auditC/probe_ox.f90`, compiled with `gfortran -O2`
against the production `oxygen_rates.f90` and run at n_H2 = 5e12, n_HI = 1e12:

```
     T        Kc(O1)         Kc(O2)        f_OH        f_H2O    f_O
     5.00       Infinity    4.82510E-42         NaN         NaN         NaN
     8.00       Infinity    1.70623E-41         NaN         NaN         NaN
    10.00       Infinity    4.67627E-41         NaN         NaN         NaN
    10.50    3.21269+306    4.04291E-39     0.00000     1.00000     0.00000
```
and `max(NaN, -700.0d0)` returns NaN on this compiler (tested).

Reaches: `oxygen_chemical_equilibrium_fractions` is the molecular-basin seed of
the cell solve and of a restart. Callers:
`ionization_equilibrium.f90:3157`, `constrained_chemical_equilibrium.f90:1379`,
`set_IC.f90:272`, `load_IC.f90:522`. A NaN seed enters hybrd1. The rate path
`rate_from_detailed_balance` is benign (`k/Inf = 0`). I did not find a
temperature floor that forbids T <= 10 K, and I did not run the code, so I
cannot say whether any configuration reaches it.

---

### The Lyman-Werner dissociation rate is evaluated at the cell's inner face while the H2O and OH rates of the same beam use the exact cell mean
Status: CONFIRMED
Severity: P2
Location: `src/modules/radiation/util_ion_eq.f90:337-346` against `:352-374`;
`src/modules/functions/utilities.f90:491-502`
```
			k_lw(j) = lyman_werner_dissociation_rate(                     &
			              fuv_band_flux(ib_LW),                           &
			              NH2col(j), T_K(j), nH_nuc(j), tau_b(j,ib_LW))
```
What the code does: `calc_column_dens_one` accumulates from the top down and
`Ncol(j)` therefore contains the *whole* of cell j
(`Ncol(j) = Ncol(j+1) + nsp(j)*dr_j(j)*R0*opa_pf(j)`). `k_lw(j)` is
`sigma_diss(NH2col(j))`, i.e. the rate a molecule would see at the cell's
inner (planet-ward) face, applied to the whole cell. In the same routine, and
sharing the same beam, the H2O and OH rates were deliberately rewritten as the
exact cell mean, `tau_out = tau_b(j+1)` plus `(1-exp(-dtau))/dtau`, and
`water_photolysis.f90:447-458` records why: the one-point rule "under-counts
the absorption everywhere else, one-signed" and was measured 30 per cent low in
the Ly-alpha band. `sigma_diss` falls by four decades across the shielding
transition, and it falls fastest exactly at the H2 front, where a cell can span
a large fraction of a decade.

What it should do: the same cell-mean treatment, or at least the half-cell
column, so that the two absorbers of one beam are discretized consistently. As
written, `k_lw` is one-signed low wherever `sigma_diss` varies inside a cell.

How it was verified: read `calc_column_dens_one` and both branches of
`fuv_lw_photon_field`; confirmed from the tabulated `sigma_diss` that the
column dependence is steep (f_shield falls from 1.0 to 4.2e-2 over
log N = 12 to 16.6 at 1300 K, i.e. about one decade of rate per 0.5 decade of
column near the knee). I did not measure the resulting error on a run.

Reaches: every molecular run with `Stellar LW flux` (regression case
`mol_lyman_werner`, and `mol_carrier`/`mol_ir_bands` if they carry a band
flux). Fixing it would move those goldens.

---

### The molecular reaction-heat ledger deposits more than the H(n=2) threshold as thermal energy in two dissociative-recombination channels
Status: CONFIRMED (the code path; the size of the effect is not measured)
Severity: P2
Location: `src/modules/lower_atmosphere/molecular_reaction_heat.f90:224-250`
```
         ! R5  H2+ + e  -> H + H
              k5 *nh2p*ne(j)          *( h_H2p )                           &
         ...
         ! R16 HeH+ + e -> He + H
            + k16*nhehp*ne(j)         *( h_HeHp )                          &
```
What the code does: the species-enthalpy table is defined with the zero at
ground-state H, ground-state He and free electrons at rest
(`species_enthalpies`, lines 268-284), so every reaction heat is the
exoergicity into *ground-state* fragments. For R5, `h_H2p = IP(H2) - D0(H2) =
15.4259 - 4.4781 = 10.948 eV` is deposited as thermal energy per dissociative
recombination of H2+; for R16, `h_HeHp = IP(H) - D0(HeH+) = 13.5984 - 1.8452 =
11.753 eV`. Both exceed the 10.199 eV excitation energy of H(n = 2), so
ground-state H + H (and He + H) is not the only accessible exit channel of
either recombination; a fragment left in n = 2 radiates 10.2 eV out of the gas
as Lyman-alpha, leaving 0.75 eV (R5) or 1.55 eV (R16) as kinetic energy.

What it should do: the module already makes exactly this distinction elsewhere
("the energy goes into the fragments, not a photon" is asserted for R6, whose
9.25 eV *is* below the n = 2 threshold, and R7 at 4.77 eV likewise), and the
Lyman-Werner block of `ionization_equilibrium.f90:2055-2081` carries a
branching for the analogous fluorescent return. R5 and R16 need either a
product-state branching or an explicit statement that the deposited heat is an
upper bound. The header's exclusion list (lines 55-76) does not mention product
excitation.

How it was verified: recomputed every enthalpy in `species_enthalpies` and
every reaction heat in `species_enthalpy_report` by hand, and checked them
against the endothermicities the rate fits imply (R10's `exp(-21900/T)` = 1.89
eV against the ledger's 1.8275 eV = IP(H2) - IP(H); R11's `exp(-20000/T)` =
1.72 eV against 1.698 eV; R8's +1.698 eV against the header's "+1.70 eV"). All
consistent, so the ledger is internally exact; the issue is the assumed product
state, not the arithmetic. I did not measure the two terms on a run; both
channels are minor in the molecular base because H2+ and HeH+ are destroyed by
the two-body reactions R8 and R18 far faster than they recombine, so I expect
the effect to be small there. The same ground-state assumption applies to the
vibrational excitation of the H2+ formed by R23 (9.16 eV) and of the H2 formed
by R6 (9.25 eV), which is thermalized at the layer's density and so costs
nothing there.

Reaches: `Molecular reaction heat` (default ON), i.e. every molecular
regression case.

---

### `oxygen_rates.f90` states it is not in the build, and it is
Status: CONFIRMED
Severity: P2 (documentation)
Location: `src/modules/lower_atmosphere/oxygen_rates.f90:5-9`
```
      ! THIS MODULE IS NOT WIRED INTO THE CODE.  It is the milestone-M1
      ! deliverable: coefficients, their sources and their validity ranges,
      ! plus the thermodynamic data that generates the reverse rates.  It is
      ! deliberately absent from SRC in the Makefile; M2 connects it.  The
      ! standalone driver that checks it is src/tests/a2_m1/.
```
What the code does: `Makefile:133` lists
`src/modules/lower_atmosphere/oxygen_rates.f90` in `SRC`, and the module is
used by `water_photolysis.f90`, `System_HeH_mol.f90`,
`diffusive_photochemistry.f90`, `ionization_equilibrium.f90`,
`constrained_chemical_equilibrium.f90`, `set_IC.f90`, `load_IC.f90` and
`write_output.f90`. The header is the state of the tree at milestone M1 and has
not been updated since M2 connected it.

How it was verified: `grep -n oxygen_rates Makefile` and a tree-wide grep of
`use oxygen_rates`.

Reaches: nothing at run time; it misleads a reader into thinking a live
physics module is a deliverable stub.

---

### The "FUV shell average" key that two comments describe does not exist
Status: CONFIRMED
Severity: P2 (documentation)
Location: `src/modules/init/parameters.f90:1208-1214` and
`src/modules/lower_atmosphere/lyman_werner.f90:153-157`
```
      ! lower, because the slant columns away from the substellar point are
      ! longer: 0.20-0.33 for the Lyman-Werner band over a hot-Uranus
      ! molecular layer (docs/e2_lw_geometry.md sec. 5.3).  That refinement
      ! is the separate, default-off "FUV shell average" key; this function
      ! is the convention it replaces when asked for.
```
What the code does: there is no such key. `input_read.f90` contains no
occurrence of "shell" in any form, there is no variable that carries a shell
average, and `dayside_dilution()` is the only dilution in the code. The
lyman_werner header makes the same claim: "that refinement is the separate,
default-off shell-average key".

How it was verified: case-insensitive grep for `shell` over `src/` and `docs/`;
the only hits are the two comments themselves and the two prose references to
`docs/e2_lw_geometry.md`.

Reaches: a reader who believes the 0.20-0.33 correction is available as an
option. The shielded-band shell average is genuinely omitted, and the run has
no way to ask for it.

---

### `Tcode` is threaded through four carrier routines and never used
Status: CONFIRMED
Severity: P2
Location: `src/modules/lower_atmosphere/diffusive_photochemistry.f90:570, 662,
847, 2253` (declarations at 572, 664, 851, 2289)
```
      subroutine carrier_state(rho, Tcode, f_sp, fc, ntot, nrho, wfac,   &
                               TK, mbar, nH_free, nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      ...
      do j = 1-Ng, N+Ng
         ntot(j) = bg_cell(j)%ntot
         TK(j)   = bg_cell(j)%T_K
      enddo
```
What the code does: `Tcode` appears only in the four signatures and the four
declarations; the temperature the carrier chemistry, the diffusivities and the
sound speed use is `bg_cell(j)%T_K`, i.e. the temperature the last ionization
sweep left. `carrier_steady_residual` is documented as answering "Is the state
(rho, v, T, f_sp) a steady state OF THE CARRIER EQUATION" while `rho` is used
and `T` is not.

What it should do: either read `Tcode*T0` (and drop the `bg_cell` read), or
drop the argument. As it stands, a caller that passes a temperature the sweep
has not seen gets a residual assembled at the previous temperature with no
diagnostic.

How it was verified: `grep -n Tcode` over the file returns exactly the twelve
lines above (four calls, four declarations, four signatures) and no use. In the
steady solver the omission is currently harmless: `eval_residual`
(`steady_newton.f90:565`) calls `ioniz_eq(T, ...)` before
`carrier_steady_residual` at line 614, so `bg_cell%T_K` is the trial T.

Reaches: no wrong number today; it is a trap for the next caller.

---

### `adv_corr` can be read uninitialized if the first call is made with `weno_mode = 2`
Status: CONFIRMED (the code path); SUSPECTED unreachable in the current caller order
Severity: P2
Location: `src/modules/lower_atmosphere/diffusive_photochemistry.f90:1063-1080`
```
      if (.not. allocated(adv_corr)) allocate(adv_corr(1:N,n_carrier_max))
      ...
      if (weno_mode .eq. 2) return
      adv_corr = 0.0d0
```
What the code does: the array is allocated before the mode test and is filled
only after it, so the first ever call at `weno_mode = 2` leaves it allocated
and undefined. `carrier_residual` then adds `ntv(j)*adv_corr(j,ic)` to every
row (line 1657) guarded only by `allocated(adv_corr)`.

What it should do: `adv_corr = 0.0d0` on allocation, or set a "filled" flag.

How it was verified: read the routine and traced every assignment of
`weno_mode` (`grep -rn weno_mode src/`). `weno_mode = 2` is set in exactly one
place, `steady_newton.f90:2115`, nine lines after `weno_mode = 1` and an
`eval_residual` call that reaches `carrier_advection_correction` through
`carrier_steady_residual`, so today the mode-1 fill always precedes the mode-2
reuse. The marching loop and `relax_photochemical_composition` run at mode 0.

Reaches: no configuration I could construct; it is one reordering away from
reading uninitialized memory into the residual of every cell.

---

### Dead and duplicated constants in `water_photolysis.f90`
Status: CONFIRMED
Severity: P2
Location: `src/modules/lower_atmosphere/water_photolysis.f90:287-288, 318-319,
364, 388`
```
      real(dp), parameter :: wl_lya_A   = 1215.67d0
      real(dp), parameter :: e_lya_erg  = 1.63403379d-11
      ...
      logical,  save :: wp_ready = .false.
```
What the code does: `wl_lya_A` and `e_lya_erg` are declared, are not public,
and are never referenced. `e_lya_erg` is a second copy of the number already
carried as `e_photon_flat_band(ib_B2)`, which is the project's own named
hazard (two definitions of one constant that cannot be caught by the compiler
when they drift). `wp_ready` is set by `water_photolysis_init` and never read,
so the guard its name implies does not exist: if `water_photolysis_init` is not
called, `eth_H2O_band` and `eth_OH_erg` stay at 0 and
`heat_per_water_dissociation` returns the whole photon energy as fragment
kinetic energy, i.e. the bond energy is deposited twice. (`input_read.f90:1807`
does call it under `thereis_oxychem`, so this is latent.) The `use oxygen_rates,
only:` list on line 287 imports `fuv_band_lo_A` and `fuv_band_hi_A`, neither of
which the module uses.

How it was verified: grep of each symbol over the file.

---

### `molecular_infrared_data` exports a public module variable named `i`
Status: CONFIRMED
Severity: P2
Location: `src/modules/lower_atmosphere/molecular_infrared_data.f90:24-25, 53`
```
      module molecular_infrared_data
      implicit none
      public
      save
      ...
      integer :: i
```
What the code does: `i` is the implied-do index of the `data` statements and is
a public, saved module variable. `molecular_infrared_cooling.f90:116` does
`use molecular_infrared_data` with no `only:` list, so `i` is in scope in every
routine of that module. `implicit none` cannot catch a routine that forgets to
declare its own `i`, because `i` is declared, and the shared copy would then be
written from inside the OpenMP cell sweep. Today no routine relies on it
(`planck_band_integral` declares its own `i`, which shadows the module one),
so nothing is wrong at run time.

What it should do: the index belongs inside the `data` statements' own scope,
or the module variable should be `private`.

How it was verified: read both modules; checked every routine of
`molecular_infrared_cooling` for an undeclared loop index.

---

### `carrier_mass_amu` hardcodes species-table indices under a comment that says it does not
Status: CONFIRMED
Severity: P2
Location: `src/modules/lower_atmosphere/diffusive_photochemistry.f90:1238-1254`
```
      ! Mass of a transported carrier [m_H], read from the species table so
      ! that this module and calc_rho cannot disagree about it.
      double precision function carrier_mass_amu(ic) result(m)
      select case (ic)
      case (ic_H2)
         m = bsp_mass(7)
```
What the code does: the masses are read as `bsp_mass(7)`, `bsp_mass(11)`,
`bsp_mass(12)`, `bsp_mass(2)` and `bsp_mass(13)` with literal positional
indices, while `species_table` exports named indices for the same species
(`isp_H2` and so on) and the `bsp_*` block is ordered by a separate convention.
The values are correct today (2.0, 16.999, 17.999, 1.0, 28.010), but inserting
a species into the `bsp` list silently changes every carrier's diffusion
coefficient and settling term with no compile error. The comment claims the
opposite property.

How it was verified: read `species_table.f90:97-108` and matched every index.

Reaches: the molecular diffusion coefficients and the settling drift of all
five carriers, i.e. `mol_carrier` and any run with
`Molecular carrier transport: True`.

---

### The deferred second-order advection correction does not cancel the donor-cell truncation exactly on a quadratic
Status: CONFIRMED
Severity: P2 (documentation)
Location: `src/modules/lower_atmosphere/diffusive_photochemistry.f90:230-235`
```
      ! for v >= 0 (and the mirror image for v < 0), with sigma the limited
      ! slope.  On a linear profile the two sigma cancel and on a quadratic
      ! they cancel the donor-cell truncation term exactly, so the scheme is
      ! second order;
```
What the code does: `carrier_slope` returns the van Leer harmonic mean
`2ab/(a+b)`, not the arithmetic (central) mean. On a uniform grid with
`f = r^2` and `r_j = j h`, the one-sided differences are `a = (2j-1)h` and
`b = (2j+1)h`, the harmonic mean is `2jh - h/(2j)`, and
`donor + (sigma_j - sigma_{j-1})/2 = 2jh + h/(4j(j-1))`. The truncation term is
cancelled to `O(h/j^2)`, not exactly; with the arithmetic mean it would be
exact. The scheme is second order either way, so only the claim is wrong.

How it was verified: worked the algebra above; read `carrier_slope`
(lines 1017-1031) and `carrier_advection_correction`.

---

### The FUV band set leaves three 1-Angstrom gaps and drops the 1202-1230 A continuum
Status: CONFIRMED
Severity: P2
Location: `src/modules/lower_atmosphere/oxygen_rates.f90:382-386`
```
      real(dp), parameter :: fuv_band_lo_A(n_fuv_band) =                  &
           (/  912.0d0, 1110.0d0, 1202.0d0, 1231.0d0, 1451.0d0 /)
      real(dp), parameter :: fuv_band_hi_A(n_fuv_band) =                  &
           (/ 1110.0d0, 1201.0d0, 1230.0d0, 1450.0d0, 2304.0d0 /)
```
What the code does: LW and B1 are contiguous at 1110 A, but 1201-1202,
1230-1231 and 1450-1451 A belong to no band. The same edges are the
user-facing definition of the four flux keys
(`parameters.f90:388-395`), so a user integrating a stellar spectrum band by
band drops 3 A of the 912-2304 A range. Separately, B2's flux key is
`F_Lya_star`, the line flux; the 1202-1230 A *continuum* under the line is
therefore not represented at all. The module header states the second point
("the band's flux really is the line") but the header's own claim, "ONE
INTERVAL, ONE FLUX, ONE BEAM ... otherwise the energy of the interval is
counted twice", is about double counting and does not cover the gaps.

How it was verified: read the arrays and `fuv_band_flux`; the band edges are
used nowhere at run time (only in the two standalone test printouts), so
nothing is miscomputed inside the code; the loss is in the definition the user
must integrate to.

---

### The tabulated fluorescent-trapping ratio does not tend to 1 in the optically thin limit
Status: SUSPECTED
Severity: P2
Location: `src/modules/lower_atmosphere/h2_self_shielding_table.f90:93-104`
and the generator `src/utils/h2_shielding_table_line_by_line.py:63-85`
```
def trapping_ratio(run_dir, T, ln):
    """p_eff/p_single of the CLOUDY run at this node, on the column axis."""
    ...
            pe = np.minimum(ps * trapping_ratio(run_dir, T, ln), 1.0)
            lsig[:, it, jn] = np.log10(sp * pe)
```
What the code does: `sigma_diss = sigma_pump * p_eff` with
`p_eff = p_single * (p_eff/p_single)_CLOUDY`. Extracting the three cubes, the
ratio `p_eff/p_single` at the *bottom* of the column axis
(N_H2 = 1e12 cm^-2, T = 1300 K) is 1.324, 1.208, 1.173 at n_H = 1e12, 1e13,
1e14 cm^-3, rising only to 2.05, 1.50, 1.39 at N_H2 = 1.55e21. At
N_H2 = 1e12 the Lyman-Werner lines are optically thin in the star-ward
direction (the table's own `f_shield` is 1.0000 there, and a strong line with
f = 0.01 and b = 3.3 km s^-1 gives tau_0 ~ 3e-3), so no fluorescent photon can
be re-absorbed on the star-ward side and the ratio should be 1. The table
therefore raises `sigma_diss` by 17 to 32 per cent at the thin end, and the
lookup is indexed by the star-ward column alone, so a cell with little H2 above
*and* below it still receives that enhancement.

What it should do: either the trapping axis is not the star-ward column (it is
the depth into a CLOUDY slab whose far side is thick, which would make the
ratio a property of that slab and not of the cell), or the ratio must tend to
1. The module's header already flags the geometry ("ILLUMINATED ON ONE FACE AND
CLOSED ON THE OTHER ... the layer is spherical and open outward") and quotes a
slab/sphere over-prediction of "up to 20 per cent where H2 is already
negligible", which is the same order as what I measure; I cannot tell from the
tree whether the two statements are the same finding, because the `.rat` CLOUDY
outputs and the deck that produced them are not in the repository.

How it was verified: parsed the three cubes out of the generated Fortran
(`auditC/`, exact reshape order checked against the declaration) and evaluated
the ratio across the column and density axes; reproduced the header's own
"13 per cent between n_H = 1e12 and 1e14 at fixed T = 1300 K" (I measure
1.3237/1.1732 = 1.128). The header's wording "p_eff/p_single -- the
SUPPRESSION of the branching by re-absorption" is also the wrong sign for a
ratio that is everywhere above 1.

Reaches: `k_lw` in every molecular run with a Lyman-Werner band flux, and
through it the 0.4 eV fragment heating and the fluorescence heating.

---

### What the A2 and E1 tests print but do not assert
Status: CONFIRMED
Severity: P2
Location: `src/tests/a2_m1/a2_m1_rate_check.f90`,
`src/tests/a2_m2/a2_m2_kinetics_check.f90`,
`src/tests/a2_m2/a2_m2_g3_run_check.f90`,
`src/tests/e1_h2/e1_h2_channel_check.f90`

None of the four programs returns a non-zero exit status under any condition,
and none of them compares a computed number against a tolerance except in the
one case noted below. Specifically:

* `a2_m1_rate_check`: sections 1 (every adopted rate against the other
  network's transcription), 2 (the whole Shomate table), 3 (the round-trip
  residual of `rate_from_detailed_balance`, which the header calls out as
  something that "must be at machine precision"), 4 (the three published
  detailed-balance tests, including the dH/dS comparison against the Baulch
  data sheet's -62.9 kJ/mol and -10.9 J/K/mol) and 5 (the measured budget) are
  all `write` statements. The *only* verdict in the file is
  `conservation_report`, which prints `CONSERVATION: PASS` or `FAIL` for a
  stoichiometry table hand-written inside the test (`cnt`, `irc`, `ipr` at
  lines 400-430): it tests the test's own table, not `oxygen_rates`. Dead code:
  `react3`/`prod3` are declared and assigned at lines 240-242 and never read.
* `a2_m2_kinetics_check`: section 2 prints the two G3 row residuals scaled by
  the largest term, which is exactly the quantity a gate would threshold, and
  thresholds nothing. Sections 0, 1 and 3 (the whole FUV band and threshold
  table, the heat per dissociation, the `<E>/<hv>` ratios) are printouts.
* `a2_m2_g3_run_check`: computes `nbad`, the number of cells whose solved
  partition deviates from chemical equilibrium by more than 1e-6, and the worst
  deviation, and prints both. It exits 0 whatever they are. It also hardcodes
  the `Ion_species.txt` layout (41 fields, column 2 = H I, column 35 = H2) and
  assumes that file and `Oxygen_chemistry.txt` have identical row counts and
  ordering; both are true today (both writers loop `1-Ng, N+Ng` and n_mion =
  27) and both are silent coupling to an output schema that carries a
  `# columns` header the test does not read. `nH2(k)` is indexed with the
  `Oxygen_chemistry.txt` row counter with no bound check.
* `e1_h2_channel_check`: test 1's channel-sum closure, test 2's `nbad` count of
  channels violating the H-nucleus and charge invariants (the A3 ledger),
  test 3's continuity across all ten joins, test 4's reproduction of the
  pre-E1 two-way split and test 5's rate-weighted ledger are all printed and
  none is compared against a bound.

Not covered at all by any of the four: the FUV attenuation and the shared-beam
split (the `a2_m2` header says explicitly "A run measures the discretization of
that identity; it is not reproduced here"); the carrier transport operator
(element closure, the limiter, the write-back, the base boundary); the
Lyman-Werner rate and its self-shielding table; the H3+ and molecular infrared
cooling; the H2 thermochemistry table.

Additionally: the regression matrix has **no oxygen-chemistry case**.
`DEFAULT_CASES` in `backup/regression/run_check.sh:187` is
`wasp_full wasp_he23off wasp_full_newton mol_base_handoff mol_metals
mol_lyman_werner mol_diffusion mol_ir_bands mol_sec_ion mol_carrier
lower_profile`, and a grep for `Oxygen chemistry` over every
`backup/regression/*/input.inp` returns nothing. So `oxygen_rates.f90`,
`water_photolysis.f90`, the OH/H2O/CO carriers of
`diffusive_photochemistry.f90` and the `oj*`/`ok*` rows are entirely outside
the byte-identical regression, and the only configuration that exercises them
is `examples/18_oxygen_chemistry/`. Every finding above that touches the oxygen
path can be fixed without any golden moving, which is convenient and is also
the reason a regression could not have caught them.

---

## Files read completely

* `src/modules/lower_atmosphere/mol_rates.f90` (647 lines)
* `src/modules/lower_atmosphere/oxygen_rates.f90` (974)
* `src/modules/lower_atmosphere/water_photolysis.f90` (509)
* `src/modules/lower_atmosphere/lyman_werner.f90` (549)
* `src/modules/lower_atmosphere/h2_self_shielding_table.f90` (1138: header and
  all executable code read line by line; the three generated numeric cubes were
  read programmatically and their contents checked against the header's claims)
* `src/modules/lower_atmosphere/h2_vibrational_relaxation.f90` (233)
* `src/modules/lower_atmosphere/h3p_cooling.f90` (252)
* `src/modules/lower_atmosphere/molecular_infrared_cooling.f90` (434)
* `src/modules/lower_atmosphere/molecular_infrared_data.f90` (2225: header and
  declarations read line by line; the `data` blocks parsed and spot-checked
  programmatically)
* `src/modules/lower_atmosphere/molecular_reaction_heat.f90` (329)
* `src/modules/lower_atmosphere/lower_column.f90` (226)
* `src/modules/lower_atmosphere/diffusive_photochemistry.f90` (2397)
* `src/tests/a2_m1/a2_m1_rate_check.f90` (486)
* `src/tests/a2_m2/a2_m2_kinetics_check.f90` (190)
* `src/tests/a2_m2/a2_m2_g3_run_check.f90` (115)
* `src/tests/e1_h2/e1_h2_channel_check.f90` (223)

Read in part, as callers or interfaces: `util_ion_eq.f90` (the FUV/LW field,
`fuv_band_flux`, the molecular cooling call site), `ionization_equilibrium.f90`
(the LW and molecular heating terms, the FUV field call, `nh`),
`System_HeH_mol.f90` (`mol_heh_rows`, `set_mol_coeffs`, `set_oxygen_coeffs`,
the threadprivate list), `steady_newton.f90` (`eval_residual`, the `weno_mode`
sequence, the carrier scales), `write_output.f90` (the Ion_species,
Lyman_Werner, Oxygen_chemistry and FUV_bands writers), `utilities.f90`
(`calc_column_dens_one`, `hydrogen_helium_nuclei_density`),
`binary_element_diffusion.f90` (`hard_sphere_pair_diffusion`),
`species_table.f90` (the `bsp_*` block), `parameters.f90` (the constants,
`dayside_dilution`, the FUV keys), `Cool_coeff.f90` (`base_sky_fraction`),
`grav_field.f90`, `input_read.f90` (the init calls and the FUV keys),
`Makefile`, `backup/regression/run_check.sh`,
`src/utils/h2_shielding_table_line_by_line.py`.

Sources consulted: `references/Koskinen_2022_ApJ_929_52.pdf` (Table 1 p. 19 and
eqs. 11-13 p. 4, via `pdftotext -layout`).

Probes written and run (in the scratch directory only, linking the production
modules): `probe_lc.f90` (lower_column against Koskinen Model A),
`probe_r12.f90` (mol_rates R12/K_eq/q_rovib against the comment's own
verification table), `probe_ox.f90` (oxygen_rates equilibrium constants and the
chemical-equilibrium fractions at low T), `parse_h2.py` (H2 line and level
data), and a parser for the three self-shielding cubes.

## Not checked

* The offline generators and their inputs: `cooling_data/molecular_infrared_bands.py`,
  `src/utils/h2_shielding_table_line_by_line.py` beyond its trapping step,
  `src/utils/h2_shielding_table_from_cloudy.py`, the CLOUDY decks under
  `src/utils/h2_shielding_cloudy_decks/`, and the `.rat` / `abs_T*.npz` files
  they read. The numbers in `h2_self_shielding_table.f90` and
  `molecular_infrared_data.f90` were checked for internal consistency and
  against the headers' own claims, not re-derived from HITEMP, the Photochem
  k-coefficients or the Abgrall line list.
* The H2O and CO band-mean cross sections were not compared against HITEMP or
  against the Photochem correlated-k files (not read).
* Individual coefficients of `oxygen_rates` against the Baulch (2005) and
  Atkinson (2004) PDFs: those two evaluations are not in `references/`, so the
  O1/O2/O9/O10/O12 rate expressions were checked only for internal consistency
  (room-temperature magnitudes, and agreement with the two reference-network
  transcriptions quoted in the module) and not against the published data
  sheets. O1 at 300 K gives 6.3e-15 and O10 at 300 K gives 1.4e-12 cm^3 s^-1,
  both the accepted values.
* The Shomate coefficients themselves against the NIST WebBook pages (no
  network access assumed); only the derived quantities the tests print were
  reasoned about.
* `System_HeH_mol`'s rows beyond row 1, 4, 9 and 10 and the argument-order
  check at the carrier's call site.
* I did not build EXHALE, did not run `EXHALE.x`, and did not run the
  regression matrix, so no finding here carries a measured effect on a
  converged model.

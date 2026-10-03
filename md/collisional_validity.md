# Collisional validity of a hydrodynamic escape solution

Tool: `src/utils/collisional_validity.py` (Python post-processing; it reads
an existing run directory and changes nothing).

A hydrodynamic wind solution is a continuum solution, and it is a solution
of the physical problem only where the gas is collisional on the scale over
which the flow varies -- above all through its critical point, which is
where the topology, and with it the mass flux, is fixed.  Where that fails,
the result is reported **unvalidated**, not overestimated: which way a
collisionless critical region moves the answer is not decided by the
continuum solution that assumed it away.

---

## 1. Collision model

The momentum-transfer collisions are the same four limits, in the same
Chapman-Enskog first approximation, that the binary element-diffusion
operator uses, so a mean free path quoted here and a diffusion coefficient
quoted by the wind cannot drift apart.

| pair | limit | coefficient | source |
|---|---|---|---|
| neutral-neutral | rigid sphere; the prefactor is a collision diameter d = 2.99 A in the Chapman-Enskog first approximation | `D = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2)/n` | Banks & Kockarts (1973); `binary_element_diffusion.f90` `hard_sphere_pair_diffusion` |
| ion-neutral, non-resonant | induced dipole (Langevin) **and** rigid core, combined as `1/D = 1/D_pol + 1/D_hs`. This is an interpolation between the polarization and rigid-core limits (it recovers each where the other is negligible), not a collision integral of the combined potential; in the crossover it can overstate the friction by up to a factor of 2 against the stronger estimate alone | `D_pol = kT/(2.21 pi e n (alpha_n mu)^(1/2))` | `binary_element_diffusion.f90` `polarization_pair_diffusion`, `ion_neutral_pair_diffusion` |
| ion-ion, ion-electron | screened Coulomb, Spitzer `ln Lambda` | `D = 3(kT)^(5/2)/[4(2 pi mu)^(1/2) n (Z_s Z_t e^2)^2 ln Lambda]` | Spitzer (1962) sec. 5.2, Paquette et al. (1986) sec. II; `binary_element_diffusion.f90` `coulomb_pair_diffusion`, `coulomb_logarithm` |

Python cannot call the Fortran, so the four routines are transcribed, with
the Fortran function or parameter named at each definition and every
constant taken from that file rather than re-chosen: `e_esu =
4.803204713e-10` (CODATA 2018), `1.52e18` (Banks & Kockarts),
`alpha(H) = 4.5 a_0^3` (exact nonrelativistic), `alpha(He) = 1.383 a_0^3`,
`alpha(H2) = 5.315 a_0^3` (Schwerdtfeger & Nagle 2019), and the helium mass
`m_He/m_H = 3.9715259` of the species table.  The two implementations are
meant to be read side by side.  (The Banks & Kockarts prefactor put through
`D = 3/(8 n d^2) (kT/(2 pi mu))^(1/2)` gives d = 2.99 A; a 2.7 A diameter
would give 1.86e18.)

**From D to a collision frequency.**  The friction force density between two
species is `F_s = (n_s n_t kT)/(n D_st) (v_t - v_s) = n_s mu_st nu_st
(v_t - v_s)`, which is the definition of the binary coefficient, so

```
nu_st = n_t k T / (mu_st Dhat_st),      Dhat_st = n D_st
```

and the total carrier density cancels: only the partner density enters.
`nu_st` is the Chapman-Enskog momentum-transfer frequency, which for rigid
spheres of diameter `d` is `(4/3) n_t pi d^2 gbar_st`, with
`gbar_st = sqrt(8kT/(pi mu_st))` the mean relative speed.  In a gas of like
particles the mean free path below is therefore 4/3 **shorter** than
Maxwell's `1/(sqrt(2) n pi d^2)` and 1.9x shorter than the elementary
`1/(n pi d^2)`.  The threshold used for the verdict (0.1) carries that
factor with room to spare, and the convention is stated so a number from
here is not compared against a differently defined one.

**Species.**  H I, H II, He I, He II, He III, e, plus H2, H2+, H3+, HeH+
when a run tracks them. The diagnostic adds HeH+ and electrons to the
Fortran H/He friction lists; HeH+ is not a friction carrier in the production
element diffusion operator. The masses use the same species table:
an H2 is one partner of mass
2 m_H).  He I 2^3S is **not** a
separate species -- it is an excited level inside the He I column
(`species_table.f90` `bsp_is_excited_level`), and counting it again would
count those atoms twice.  **Metals are excluded as collision partners**,
and so are the oxygen-chemistry molecules OH, H2O and CO.  They are kept in
the electron budget, where charge neutrality has to close exactly, and the
tool computes and prints their largest mass and particle fractions on the
profile it reads (atomic weights from `species_table.f90` `melem_A_u`).
Measured on the solution profiles: 1.1e-3 by mass and 2.1e-4 by particle
count on `LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH2.09`;
2.6e-3 and 3.0e-4 on `LHS1140b/archive_20260830/exhale/flux_closure/hi/k06`;
1.1e-2 and 8.3e-4 on the HD 209458 b state
`backup/phase_d_baseline/new_kzz1e9_d3b` (the `_adv` pairs give the same
maxima to the digits quoted).

**Omitted channels, all toward a longer mean free path, i.e. a larger
Knudsen number, i.e. a more conservative verdict within the retained
collision model.** This monotonic statement holds with retained coefficients
fixed. It does not prove an upper bound on the physical mean free path:
the polarization-plus-core interpolation and the estimates for resonant
pairs still require independent validation. Omitting a collision
channel can only lower a collision frequency, because the frictions of
separate partners add; every mean free path here is therefore an
overestimate of the one with the channel included.  (i) The metals and
oxygen-chemistry molecules as partners (above).  (ii) *Resonant* charge
exchange (H+ + H, He+ + He) is not a channel: like-element ion-neutral
pairs are computed with the non-resonant (polarization + core) coefficient,
whose friction is the smaller one, so the H I and H II paths in the
partially ionized layer are upper limits.  The size of the resonant
momentum-transfer cross section was not checked against a published table
in this repository.  (iii) Electron-neutral momentum transfer is not in the
model, so the electron Knudsen number is meaningful only where the gas is
ionized -- reported with that caveat, and never used for the bulk, which is
built from the heavy particles.

## 2. Definitions

```
nu_s     = sum_t nu_st                          total momentum-transfer rate
lambda_s = vbar_s / nu_s,  vbar_s = sqrt(8kT/(pi m_s))     (harmonic sum:
                                             1/lambda_s = sum_t 1/lambda_st)
1/L^2    = (d ln rho/dr)^2 + (d ln T/dr)^2 + (|dv/dr|/c_s)^2 + (1/r)^2
L        = max( 1/sqrt(that), dr_cell )         local structure scale
Kn_s     = lambda_s / L
lambda_bulk = mass-weighted mean of lambda_s over the heavy species
```

**Which fields enter, and why the sum of squares.**  The continuum closure
expands about a local Maxwellian, and a local Maxwellian is fixed by
`(n, T, u)`.  Those three gradients -- logarithmic for `n` and `T`, measured
against the thermal speed for `u`, see below -- plus the spherical
divergence the geometry adds, are the complete set, and the change of the
state over a distance `d` is the vector `(d/H_rho, d/H_T, ...)` whose length
is the RMS above.  Two properties follow, both wanted: `L` is never longer
than the shortest individual scale, since every term enters `1/L^2` with a
positive sign; and a field that goes logarithmically flat contributes zero
instead of removing itself from a minimum.  A minimum rule behaves the
opposite way -- when the field it is currently taken from flattens, `L`
jumps to the next one, discontinuously and by whatever factor separates
them.

**Pressure is not in the list.**  `p = n k T` is a derived field; adding it
would count the density and temperature gradients a second time wherever
they do not cancel, and it contributes nothing exactly where they do -- at
an isobaric front, which is precisely where it fails as a structure scale.
That failure was in this tool until 2026-08-28, when `L` was
`min(H_p, L_v, r)`: across the heating peak of the LHS 1140 b runs
`d ln p/dr` passes through a broad near-zero, because `rho` falls and `T`
rises with nearly the same log slope.  In
`LHS1140b/archive_20260830/exhale/heh0p55`, `H_p`
ran from 9.13e6 cm at 1.0244 `R_p` to **3.05e8 cm** at 1.0486 `R_p`, while
at 1.0453 `R_p` -- where it is already 2.55e8 -- `H_T` is 1.80e7 and
`H_rho` 1.68e7 cm, fifteen times shorter; `L` followed the pressure,
and `Kn` showed a spurious dip -- `Kn(H I)` 3.67e-5 at 1.0244, 6.82e-6 at
1.0453, 3.56e-5 at 1.0648, a 5.4x hole in a quantity that is monotonic on
either side of it.  With the definition above `Kn(H I)` rises
monotonically from 1.02 to 3 `R_p` (4 sign changes of `d Kn/dr` before, 0
after).

**Velocity enters as `|dv/dr|/c_s`, not as `|(1/v) dv/dr|`.**  The
first-order term the closure drops is the viscous stress, whose size
relative to the pressure is `~ lambda |dv/dr| / vbar`: a change of the bulk
velocity distorts the distribution function in proportion to the *thermal*
speed, not to the local bulk speed.  Normalizing by `v` diverges at every
stagnation point -- the cell-centered velocity of the base cells crosses
`v = 0` repeatedly, through the collocated two-cell odd-even mode those cells
carry (a mode of the cell-centered field, not a wave of the solution: the
Riemann face mass flux there is the wind's own, MEASURED 2026-09-24,
`md/lhs1140b_element_row_terms_20260924.md`) -- which is what an ad-hoc
Mach-number floor used to patch; that
floor is gone.  `c_s` stands in for `vbar` so that `L` remains one length
common to all species, the two differing by an O(1) factor.

**The floor at the local cell width** is a statement about what a discrete
solution can carry: no structure exists below one cell, so a gradient
claiming one is measuring the mesh.  It binds in 3 of 500 cells in one of
the four runs below, all at `r < 1.002 R_p`, where `Kn ~ 1e-5` -- four
orders of magnitude below the verdict threshold -- and in no critical
region.  The residual cell-to-cell scatter of `Kn` at `r < 1.01 R_p` is the
same odd-even mode of the cell-centered fields of the base cells, present in
`rho`, `T` and `p` alike; it is not a defect of the definition and it is far below any
radius the verdict uses (every critical region here starts at the heating
peak, 1.05-1.11 `R_p`).

| quantity | definition |
|---|---|
| **exobase** `r_exo` | `Kn_bulk = 1`, first crossing from below. If `Kn_bulk < 1` throughout, reported as *above the outer boundary*, never extrapolated. |
| **critical point** `r_s` | the sonic point, `v = c_s`, `c_s = sqrt(gamma p/rho)`, `gamma = 5/3`; `p`, `rho` are the cgs columns of the output file, so the mean molecular weight is whatever the solution carries. In atomic gas this is the **code's own** sound speed (`eval_dt.f90`, `cs = sqrt(gamma_ad*p/rho)`, `gamma_ad = 5/3` in `parameters.f90`). Where H2 is present and the caloric EOS is active the code uses `gamma_eff(T, composition) < 5/3` (`caloric_eos.f90`) and this tool does not, so its `c_s` there is too large by `sqrt(5/3 / gamma_eff)`; the effect on the verdicts was not measured. |
| **critical region** | peak volumetric heating to the sonic point -- the region that heats and accelerates the wind and sets its topology. |
| **collisional** | `max Kn < 0.1` across the critical region, over the bulk and over every species carrying at least 1e-3 of the mass. |
| **coupling times** | `tau_s = 1/nu_s` against the flow time `r/|v|`; electron-ion energy coupling `nu^E_st = 2 mu_st/(m_s+m_t) nu_st` against the heating time `(3/2) n_tot k T / Q`. |

Coupling time and Knudsen number are **not** the same test and can point
different ways: `tau_s/tau_flow = Kn_r x Mach x (v/vbar)`-like, so a
strongly subsonic flow gives many collisions per flow time even at
`Kn ~ 1`.  The Knudsen number is what decides continuum validity, because
the terms the continuum closure drops (viscous stress, heat flux, the
non-Maxwellian tail) scale with the gradient over a mean free path, not with
the collision count per flow time.  Both are reported.

## 3. Data read

`output/Hydro_ioniz.txt` and `output/Ion_species.txt` -- the **solution**,
because the sound speed, the Mach number, the critical point and the
structure scale are properties of the state that satisfies the momentum
equation, and the `_adv` pair carries that state's density and velocity
beside a temperature from a closure the momentum equation was never
re-solved for.
`--adv` reads the `_adv` pair instead, which is the question to ask of the
COMPOSITION, since the mean free path is set by how much of the gas is
neutral.  Until 2026-09-13 the default was the `_adv` pair and the flag was
`--eq` with the opposite sense; every table below states which pair it was
measured on.  The
loaders are `examples/exhale_io.py`; no second column parser was written.
The `Ng = 2` ghost cells at each end (`parameters.f90` `Ng`) are dropped.

**Planet radius and mass.**  The radius the profiles are scaled by, and the
mass in the escape parameter, are those the run was solved with.
`exhale_io.load_run` reads them from the run's `EXHALE_resolved.out`
(`planet_radius_RJ`, `planet_mass_MJ`, written by `write_resolved_config` in
`write_setup_report.f90`) and falls back to `input.inp` only when that file
is absent, with a printed note; a resolved file without those two keys, or
whose radius disagrees with the `R0[cm]` the profile header states, is
refused.  The lower-atmosphere profile handoff moves the base radius to its
matching level (`input_read.f90`, `lap_value_at_match('r', ...)`): on
`LHS1140b/models/molecular_photochem_gj1132_kzzprofile/HeH2.09` the solved
radius is 0.1625085 R_J against `Planet radius` 0.157692 R_J in `input.inp`
(3.05% larger), which moves `log10 Mdot` from 8.033 to 8.059 (the run's own
report: 8.06).  The report prints the radius and mass used and their source.
`log10 Mdot` is `exhale_io.mdot_log10`, i.e. `4 pi rho v r^2` twenty cells
from the top with the `2D approximate method` factor, so it can differ by
a few per cent from a number quoted from a flux-window median.

## 4. Results

Run as
`python3 src/utils/collisional_validity.py <case_dir> [<case_dir> ...]`.
The tables of this section are historical: they were measured with an
earlier version of this tool (input.inp radius, helium mass 4 m_H) on the
LHS 1140 b outputs now under `LHS1140b/archive_20260830/exhale/` and on
`backup/phase_d_baseline/new_kzz1e9_d3b`.  The current tool on those
outputs gives the same verdicts but not the same digits (exobase of
`heh0p55` 23.84 against 23.73 `R_p`, of `flux_closure/hi/k06` 19.40
against 18.84), so each number describes the state as it was measured.

| case | planet | peak heat | T max | critical point | exobase | `Kn = 0.1` at | max Kn in crit. region | verdict |
|---|---|---|---|---|---|---|---|---|
| `LHS1140b/archive_20260830/exhale/heh0p55` (diffusion off) | LHS 1140 b | 1.053 | 1.490 (4532 K) | **none in domain**, max Mach 0.540 at 29.05 | 23.73 | 6.72 | 2.3 (H I) | **criterion not met** |
| `LHS1140b/archive_20260830/exhale/heh2p13_diff_kzz1e9` | LHS 1140 b | 1.050 | 1.398 (5330 K) | **none in domain**, max Mach 0.565 | 25.44 | 7.10 | 2.0 (H I) | **criterion not met** |
| `LHS1140b/archive_20260830/exhale/flux_closure/hi/k06` (converged closure) | LHS 1140 b | 1.105 | 1.340 (3944 K) | **none in domain**, max Mach 0.550 | 18.84 | 5.80 | 3.1 (H I) | **criterion not met** |
| `backup/phase_d_baseline/new_kzz1e9_d3b` | HD 209458 b | 1.111 | 1.500 (8353 K) | **4.085** (2574 K) | above 4.15 (`Kn_top = 2.0e-3`) | -- | **0.019** (H I) | **criterion met** |

Radii in `R_p`.  The LHS 1140 b domain ends at 30 `R_p`, HD 209458 b's at
4.15 `R_p`.

`Kn_bulk` along the radius:

| `r/R_p` | 2 | 4 | 6 | 8 | 9.5 | 15 | 20 |
|---|---|---|---|---|---|---|---|
| `heh0p55` | 0.0040 | 0.029 | 0.078 | 0.144 | 0.203 | 0.461 | 0.747 |
| `heh2p13_diff_kzz1e9` | 0.0037 | 0.026 | 0.069 | 0.129 | 0.182 | 0.413 | 0.661 |
| `flux_closure/hi/k06` | 0.0056 | 0.040 | 0.108 | 0.203 | 0.288 | 0.675 | 1.11 |
| HD 209458 b `new_kzz1e9_d3b` | 6.4e-4 | 2.0e-3 (at 4.085) | -- | -- | -- | -- | -- |

Species split in the critical region (`heh0p55`): H I 2.3, He I 1.7,
H II 3.8e-5, He II 1.9e-5, He III 5.0e-6, e 6.9e-4.  The ions and electrons
are held by Coulomb collisions two to five orders of magnitude more tightly
than the neutrals; **the neutral hydrogen sets the Knudsen number**, and it
is the species the He 10830 diagnostic is about.

Coupling times (`heh0p55`): `tau_HI` = 0.5 s at 1.5 `R_p`, 942 s at 8
`R_p`, 1.6e4 s at 20 `R_p`, against a flow time `r/|v|` of 1.3e6, 1.9e5,
1.8e5 s -- so the gas is still momentum-coupled several hundred times per
flow time at 8 `R_p` even where `Kn = 0.14`.  Electron-ion energy coupling
is fast everywhere it matters: `max tau^E_ei/tau_heat` over the critical
region is 0.030, 0.023, 0.018 (LHS 1140 b cases) and 8.3e-5 (HD 209458 b),
so the single-temperature energy equation is not what fails.

### The headline result

**LHS 1140 b's wind never reaches its critical point inside the
computational domain.**  In all three representative solutions the flow is
subsonic out to 30 `R_p` (max Mach 0.54-0.57), so the mass flux is set by
the outer boundary condition rather than at a critical point -- and by 30 `R_p` the
gas is already collisionless (`Kn_bulk` = 1.4-2.4 there, exobase at
18.8-25.4 `R_p`).

[Note: the three runs above predate the 2026-09-13 re-solve of
the LHS 1140 b catalog with `Well balanced: True`, whose outer wind is
denser. MEASURED with this tool on current states (the certified
`molecular_scalar_gj1132_kzz1e9/HeH2.13` and `HeH9.7`,
`atomic_scalar_gj1132_kzz1e9/HeH2.13`, `atomic_scalar_gj1132_wellmixed/HeH0.40`,
and the uncertified `molecular_scalar_gj1132_kzz1e9/HeH0.55` state of the
2026-09-21 alternation audit): still no critical point in the domain (max
Mach 0.45-0.52), but `Kn_bulk` at the top cell (29.03 `R_p`) is only
0.055-0.142 and the exobase lies ABOVE the domain; the largest species
Knudsen number in the critical region is 0.41-0.59, and `Kn_bulk` passes
0.1 at 8.7 `R_p` on the audit state. So the outer wind of the current
states is transitional, not collisionless, and the verdict (criterion not
met) stands on that ground. The numbers of this section describe the runs named
in its table.]

Against the p-winds retrieval, whose isothermal Parker solution places the
sonic point at 8-9.5 `R_p`: that radius is **inside** the exobase computed
here (18.8-25.4 `R_p`), but `Kn_bulk` there is already 0.13-0.29, i.e. well
past the collisional threshold.  So the answer is neither "the
sonic point is outside the exobase" nor "the critical region is
collisional": the critical region of LHS 1140 b is **transitional**.  The
scale separation the continuum equations need has been lost by the radius
where either code sets the mass flux, while the gas is not yet free
molecular.

HD 209458 b is the control: sonic point at 4.085 `R_p` with `Kn` reaching
only 0.019 anywhere between the heating peak and the critical point, and no
exobase inside the domain -- the collisional criterion is met through its
critical point.  The contrast is a factor of ~100 in `Kn`, and it is what
one expects from the two winds' scales, not a marginal call.

### Sensitivity to which profile pair is read

The `_adv` and the solution profiles differ in more than the ionization
split -- `post_process_adv` also corrects the temperature, hence the
pressure and the sound speed -- so the diagnostic is not indifferent to the
choice. Read on the solution, which is the default since 2026-09-13:

| case | critical point | exobase | `Kn = 0.1` at | `Kn_bulk` at the top | max Kn |
|---|---|---|---|---|---|
| `heh0p55` | none (max Mach 0.822) | above 30 | never (`Kn_bulk < 0.1`) | 0.097 | 0.98 |
| `heh2p13_diff_kzz1e9` | none (max Mach 0.824) | above 30 | 28.77 | 0.101 | 0.86 |
| `flux_closure/hi/k06` | none (max Mach 0.823) | above 30 | 7.85 | 0.191 | 1.25 |
| HD 209458 b | 3.769 (max Mach 1.342) | above 4.15 | -- | 5.2e-4 | 0.013 |

The verdicts do not move: LHS 1140 b has no critical point in the domain on
either pair and reaches `Kn ~ 1` in the region where the outer boundary sets
its flux, HD 209458 b is collisional through its critical point on both. The
exobase radius does move -- the equilibrium profiles keep the far field more
ionized, and Coulomb collisions there are far stronger than the neutral ones
-- so an exobase radius quoted from this tool must state which pair it came
from. The tables above are the `_adv` pair (`--adv`) unless marked
otherwise.

### Formal static-Maxwellian (Jeans) comparison

The Jeans (1925) flux of a static Maxwellian, evaluated at the diagnosed
exobase or, flagged, at the outer boundary:

```
Phi_J,s    = n_s vbar_s/4 (1 + lambda_J,s) exp(-lambda_J,s)
vbar_s     = sqrt(8 k T/(pi m_s))          (mean speed)
lambda_J,s = G M_p m_s/(k T r_exo)
Mdot_Jeans = 4 pi r_exo^2 sum_s m_s Phi_J,s
```

With the most probable speed `v_mp = sqrt(2kT/m_s)` the prefactor reads
`n_s v_mp/(2 sqrt(pi))`.  It is the outward flux of an isotropic Maxwellian
through a surface, integrated over `v_r > 0` and speeds above the escape
speed (checked against direct numerical integration of that flux).  The
historical rates below were computed with `vbar_s/(2 sqrt(pi))`, the
most-probable-speed factor applied to the mean speed, which overstates every
formal rate by `2/sqrt(pi) = 1.128`.

Historical values (the `_adv` pair, earlier tool, states named in the
results table).  `Mdot_hydro` and the other columns describe those states as
measured; the formal Jeans rate and the ratio carried the factor 1.128 and
are superseded:

| case | `r_exo` | `T(r_exo)` | `lambda_J(H)` | `lambda_J(He)` | `Mdot_Jeans` | `Mdot_hydro` | ratio |
|---|---|---|---|---|---|---|---|
| `heh0p55` | 23.73 | 1233 K | 0.84 | 3.36 | ~~1.8e7 g/s~~ superseded | 10^7.764 | ~~3.2~~ superseded |
| `heh2p13_diff_kzz1e9` | 25.44 | 1158 K | 0.83 | 3.33 | ~~2.0e7 g/s~~ superseded | 10^7.806 | ~~3.3~~ superseded |
| `flux_closure/hi/k06` | 18.84 | 853 K | 1.53 | 6.11 | ~~1.0e7 g/s~~ superseded | 10^7.474 | ~~3.0~~ superseded |
| HD 209458 b | (no exobase; evaluated at 4.15) | 2508 K | 10.85 | 43.4 | ~~1.1e6 g/s~~ superseded | 10^10.056 | ~~1.0e4~~ superseded |

The same four states measured with the current tool (`--adv`, corrected
prefactor, the solved radius from `EXHALE_resolved.out`, helium mass
`m_He/m_H = 3.9715`), on the outputs at the paths of the results table:

| case | `r_exo` | `T(r_exo)` | `lambda_J(H)` | `lambda_J(He)` | `Mdot_Jeans` | `log10 Mdot_hydro` | ratio |
|---|---|---|---|---|---|---|---|
| `heh0p55` | 23.84 | 1239 K | 0.81 | 3.23 | 1.69e7 g/s | 7.776 | 3.5 |
| `heh2p13_diff_kzz1e9` | 26.12 | 1163 K | 0.79 | 3.14 | 1.92e7 g/s | 7.824 | 3.5 |
| `flux_closure/hi/k06` (solved radius 1.0299 x input) | 19.40 | 854 K | 1.41 | 5.59 | 1.04e7 g/s | 7.519 | 3.2 |
| HD 209458 b `new_kzz1e9_d3b` | (no exobase; evaluated at 4.15) | 2508 K | 10.61 | 42.1 | 1.25e6 g/s | 10.075 | 9.5e3 |

The expression assumes a stationary Maxwellian at the evaluation radius
(zero bulk drift), a point-mass potential and independent escape of each
species; the tool applies it to all heavy species, ions included, without
the ambipolar electric field.  Its ratio to the continuum mass flux is a
comparison of two formulas on the same profile.  It is not a bound on a
kinetic escape rate, not an error estimate for the continuum solution, and
not by itself evidence for a particular escape regime.  On LHS 1140 b,
`lambda_J(H) = 0.8-1.4` at the exobase: the escaping part is not a small
tail of the Maxwellian there, and the usual reading of the Jeans rate does
not apply.  On HD 209458 b, where `lambda_J(H) = 10.6` and the collisional
criterion is met through the critical point, the continuum rate is about
1e4 times the formal Jeans rate at the outer boundary.

### What the 2026-08-28 change of `L` moved

Same runs, same outputs, same collision model; only the structure scale of
section 2 changed, from `min(H_p, L_v, r)` to the state-gradient RMS.
Nothing that carries a conclusion moved.

| quantity (`_adv` pair) | old `L = min(H_p, L_v, r)` | new `L` |
|---|---|---|
| HD 209458 b, max Kn in the critical region | 0.0214 | 0.0186 |
| HD 209458 b, sonic point | 4.085 `R_p` | 4.085 `R_p` (`L` does not enter) |
| HD 209458 b, `Kn_bulk` at the outer boundary | 2.28e-3 | 2.01e-3 |
| LHS 1140 b, max Kn in the critical region (H I) | 1.9-3.0 | 2.0-3.1 |
| LHS 1140 b, max `Kn_bulk` | 1.3-2.3 | 1.4-2.4 |
| LHS 1140 b, exobase | 19.8-26.7 `R_p` | 18.8-25.4 `R_p` |
| LHS 1140 b, `Kn_bulk` at 8-9.5 `R_p` | 0.14-0.30 | 0.13-0.29 |
| LHS 1140 b, `Kn_bulk` at 1.1 `R_p` | 4e-5 - 1.2e-4 | 2.5e-4 - 4.2e-4 |
| LHS 1140 b, formal `Mdot_Jeans` | superseded (prefactor 1.128 too large) | superseded (prefactor 1.128 too large) |
| max Mach, `log10 Mdot`, coupling times | -- | unchanged (independent of `L`) |

Two directions, both expected.  Inside 1.1 `R_p` the new `Kn` is up to 3x
larger, because that is where the old `H_p` sat on its isobaric plateau.
Outside the heating peak the new `Kn` is 10-15 % *smaller* on HD 209458 b
and 5-10 % larger on LHS 1140 b: where `rho` and `T` fall together, as in
the HD 209458 b outer wind, `|d ln p/dr| = |d ln rho/dr + d ln T/dr|`
exceeds the RMS of the two, so the old scale was shorter there by up to
sqrt(2) -- correctly conservative by accident, on a field that is not an
independent state variable.

## 5. What this does and does not settle

Settled by measurement: the collision model, the species-resolved Knudsen
numbers and coupling times, the exobase and the critical point of the four
solutions above, and the verdicts in the table.

Not settled, and not claimed: what the LHS 1140 b mass-loss rate actually
is.  That needs a kinetic or transitional-flow calculation (DSMC, or a
13-moment closure carried through the transitional region), which this
repository does not contain.  Also untouched: whether the He 10830 line
metrics -- which are formed at 1.1-3 `R_p`, where `Kn <= 0.019` in every case
here (0.0003 at 1.1, 0.004-0.006 at 2, 0.012-0.019 at 3) -- are affected at all.  They are formed well inside the collisional
region; what is unvalidated is the *wind solution's mass flux*, not the
line-forming layer.

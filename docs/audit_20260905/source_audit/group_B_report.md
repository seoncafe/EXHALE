# Audit group B: radiation, photoionization, heating, cooling, atomic rates
EXHALE v1.00, HEAD 35d9dd5. Read-only audit, 2026-09-05.

All numbers below were produced by probes in
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/ff16b3fd-4a9b-42d7-8a3e-10c562168851/scratchpad/auditB/`
(`grid_probe.py`, `attn.py`, `metal_probe.py`), which re-implement the exact
grid construction of `set_energy_vectors`, the exact cross-section functions of
`cross_sec.f90` and the exact `J_inc` normalization, and compare the code's
quadrature against a 4e5-point reference integral. Nothing in the repository was
run or modified.

---

### Threshold bin at 13.6 eV charges the H I cross section over energies where it is zero: P_HI is 9% too high

Status: CONFIRMED
Severity: P0 (default path: every regression case with `Include He23S? True` or an active low-IP metal, i.e. all eight matrix cases)
Location: `src/modules/init/set_energy_vectors.f90:56-67`, `:78-92`, `:118-121`

```
56			NlTR      = 0
57			e_sub_low = e_th_HI
58			if (thereis_HeITR) then
59				NlTR      = 20
60				e_sub_low = e_th_HeTR
...
80				do j = 1,NlTR
81					e_v(j) = e_sub_low*	&
82						  (e_th_HI/e_sub_low)**((j-1.0)/(NlTR*1.0))
...
119		de_v(1)      = 0.5*(e_v(2) - e_v(1))
120		de_v(2:Nl-1) = 0.5*(e_v(3:Nl) - e_v(1:Nl-2))
121		de_v(Nl)     = 0.5*(e_v(Nl) - e_v(Nl-1))
```

**What the code does.** The photon grid puts `e_v` at the ionization thresholds
themselves (13.6, 24.6, 54.4 eV are grid points by construction) and gives each
point the central-difference width `de_v(k) = (e_v(k+1) - e_v(k-1))/2`. Every
integral is then the midpoint sum `sum(f*de_v)`. At a threshold point the lower
half of the bin lies BELOW the threshold, where the cross section is identically
zero, but the sum charges the whole width at the threshold value of sigma.

Below 13.6 eV the grid has only `NlTR = 20` logarithmic points spanning
4.8-13.6 eV (or 5.139-13.6 eV when Na I sets the floor, which is what the
regression `metals.inp` does), so the last spacing below the H I edge is 0.69-0.85 eV
against 0.16 eV above it. The bin centred on 13.6 eV is then 0.404-0.426 eV wide
where only 0.081 eV of it is physical.

**What it should do.** A threshold must sit on a bin EDGE, or the straddling bin
must carry the cross section averaged over its physical part. `sum(sigma(E) F(E)/E dE)`
is not approximated by `sigma(E_th) F(E_th)/E_th * (full bin)`.

**How it was verified.** `grid_probe.py` / `metal_probe.py` rebuild the grid and
integrate against a fine reference (default tutorial spectrum, PLind = -1,
[13.6, 123.98, 1240] eV):

| quantity | grid floor 4.80 eV | grid floor 5.139 eV | grid floor 13.6 eV |
|---|---|---|---|
| P_HI code/exact | **1.0951** | **1.0891** | 1.0003 |
| P_HeI code/exact | 1.0157 | 1.0157 | 1.0157 |
| P_HeII code/exact | 1.0299 | 1.0299 | 1.0299 |
| Heat_HI code/exact | 0.9998 | - | 0.9998 |

The heating is untouched because the photoelectron share `1 - E_th/E` vanishes at
the threshold; the error is confined to the RATE, i.e. to the ionization balance,
n_e, and everything that scales with them. `attn.py` shows how it survives
attenuation:

```
N_HI= 0.0e+00  code/exact 1.0951
N_HI= 1.0e+17  code/exact 1.0720
N_HI= 3.0e+17  code/exact 1.0378
N_HI= 1.0e+18  code/exact 1.0018
```

so the full 9-10% bias applies through the whole ionized wind and dies only once
the 13.6 eV bin is individually optically thick.

The +1.6% (He I) and +3.0% (He II) entries are the same defect on the smooth part
of the grid, where the spacings on either side of the threshold are comparable;
they are present with or without the sub-Lyman grid.

Reaches: `PH_heat_HHe`/`PH_heat_H` -> `P_HI`, `P_HeI`, `P_HeII` -> every
ionization solve, `dp_bc`, the cooling through n_e. Every current regression case
(all use `Spectrum type: Power-law`; `wasp_he23off` still has `NlTR = 20` because
its `metals.inp` carries Na, Ca, Mg and Fe, all with sub-13.6 eV neutral edges).

---

### The sub-Lyman-continuum field is the EUV power law extrapolated two decades below its normalization band; it is 500x too faint for He 2^3S photoionization in the same run whose stellar Teff the code already uses

Status: CONFIRMED (code path and numbers); the choice of which spectrum is right is a modelling decision, but the run carries two incompatible descriptions of one band
Severity: P0 (default path with `Spectrum type: Power-law` + `Include He23S? True`; changes the He 10830 observable)
Location: `src/modules/radiation/J_inc.f90:40-54`, `src/modules/init/set_energy_vectors.f90:78-83`, `:258-261`

```
J_inc.f90:45      JEUVnorm = J_EUV*P1/(e_mid**P1 - e_low**P1)
J_inc.f90:50   if(E.lt.e_mid) then
J_inc.f90:51      J_inc = JEUVnorm*E**PLind
```

**What the code does.** `JEUVnorm` is fixed by requiring the integral of `J_inc`
over `[e_low, e_mid]` to be `J_EUV = 10^LEUV/(4 pi a^2)`. `set_energy_vectors`
then evaluates `J_inc` on grid points down to `e_sub_low` (4.8 eV for the He
triplet, or the lowest active low-IP metal edge), i.e. it extrapolates the
normalized EUV power law two decades below the band it was normalized on, and
that extrapolated field is the entire radiation the He 2^3S metastable and the
low-IP metals see.

**What it should do.** At 4.8-13.6 eV (2580-912 A) the stellar flux is
photospheric, not coronal. The same `input.inp` states `Stellar Teff` and
`Stellar radius`, and `excited_hydrogen.f90:gamma_n2_balmer` already integrates
exactly `pi B_nu(T_eff) (R_star/a)^2` over 3.4-13.6 eV for the H(n=2) Balmer
continuum. One run therefore describes one band two ways.

**How it was verified.** With the `wasp_full` regression input (T_eff = 6459 K,
R_star = 1.458 Rsun, a = 0.02544 AU, LEUV = 30.42, LX = 29.46, PLind = -1):

| E [eV] | pi B_nu (R/a)^2 [erg cm^-2 s^-1 eV^-1] | code power law | ratio |
|---|---|---|---|
| 4.8 | 2.24e8 | 1.36e5 | **1.6e3** |
| 6.0 | 5.06e7 | 1.09e5 | 464 |
| 8.0 | 3.30e6 | 8.17e4 | 40 |
| 10.0 | 1.77e5 | 6.54e4 | 2.7 |
| 13.0 | 1.78e3 | 5.03e4 | 0.035 |

Integrating the code's own `sigma_HeI23S` (two-wing VFKY96 form, 4.68 Mb at
threshold) over 4.8-13.6 eV with the dayside dilution of that case (Rate/2,
xi = 0.5):

```
P(He 2^3S)  with the code's extrapolated power law : 8.98e-02 s^-1
P(He 2^3S)  with pi B_nu(6459 K) (R*/a)^2          : 4.45e+01 s^-1     (x496)
```

At 0.09 s^-1 the metastable photoionization is comparable to its other sinks
(A31 = 1.27e-4 s^-1, n_HI Q31 ~ 1e-3..1, n_e(q31a+q31b) ~ 1e-2..10); at 44.5 s^-1
it would dominate everywhere. The same band carries the neutral-edge photoionization
of Na I (5.14 eV), Ca I (6.11), Mg I (7.65) and Fe I (7.90).

Second, smaller consequence of the same extrapolation: the flux the run actually
integrates over its grid is **1.4654 x** the nominal `J_XUV = (10^LX+10^LEUV)/(4 pi a^2)`
that `write_setup_report` prints (measured, tutorial input, 4.8 eV floor); with
the floor at 13.6 eV it is 0.9963 x. `J_XUV` is also what `wae_exhale_bridge`
hands Wind-AE as `Ftot_t`.

Reaches: `P_HeITR` and hence n(2^3S) and the He 10830 transit depth; the neutral
fractions of every low-IP metal; `examples/`, `benchmarks/` and all regression
cases (all power-law). A LOADED SED is NOT affected: `sed_read` lowers `e_low`
to the same floor and reads the file's own flux there.

Known-item note: `docs/lower_atmosphere_coupling.md:626-631` records the
extrapolation in the context of the Lyman-Werner band (which is handled by a
separate user key). Its consequence for He 2^3S and the low-IP metals, and the
size of it, are not recorded anywhere I found.

---

### C I electron-impact ionization: the Voronov (1997) row is transcribed with P and X shifted, and the rate is 1.5-2.6x too high

Status: CONFIRMED
Severity: P1 (metals are opt-in through `metals.inp`, but every metals-on case carries carbon)
Location: `src/modules/radiation/Cool_coeff.f90:3004-3007`

```
3004   real*8, parameter :: dE = 11.26, A = 6.85e-8, P = 0.193,        &
3005                        X = 0.25,   K = 0.25
3006   U(j_lo:j_hi) = dE/(kb_eV*T(j_lo:j_hi))
3007   a_ion_coeff_CI(j_lo:j_hi) = A*(1.0+P*sqrt(U(j_lo:j_hi)))/(X+U(j_lo:j_hi))*U(j_lo:j_hi)**K*exp(-U(j_lo:j_hi))
```

**What the code does.** Evaluates `k = A (1 + 0.193 sqrt(U)) U^0.25 exp(-U)/(0.25 + U)`.

**What it should do.** The Voronov (1997, ADNDT 65, 1) row for neutral carbon is
`(dE, P, A, X, K) = (11.3, 0, 6.85e-8, 0.193, 0.25)`, i.e.
`k = 6.85e-8 U^0.25 exp(-U)/(0.193 + U)` with NO `(1 + P sqrt(U))` factor.
The code has taken X's value (0.193) as P and K's value (0.25) as X.
Independent confirmation in the workspace: `p-winds/p_winds/carbon.py:160-162`

```python
    ionization_rate_ci = 6.85E-8 * (0.193 + energy_ratio_ci) ** (-1) * \
        energy_ratio_ci ** 0.25 * np.exp(-energy_ratio_ci)
```

Every other Voronov row in this file is transcribed correctly (checked C II, N I,
N II, O I, O II, Mg I, Mg II, Si I, Si II, Ca I, Ca II, Na I, K I, S I, Fe I,
Fe II, and H I / He I / He II in `voronov_ci`); C I is the only one shifted.

**How it was verified.** Direct evaluation of both forms:

```
T=  2000  ratio code/Voronov 2.558
T=  5000  1.982
T=  8000  1.774
T= 10000  1.690
T= 15000  1.560
T= 20000  1.481
T= 50000  1.286
```

For scale, `alpha_rec_metal('CI', 1e4 K)` = 8.5e-13 cm^3 s^-1 (Badnell RR+DR,
hand-evaluated from the coefficients in this file), so the collisional term moves
from 2.4% to 4.1% of the recombination rate at 1e4 K; it becomes the dominant
ionization term above ~2e4 K. I did NOT measure the resulting change in the
carbon ionization structure or in the C I / C II line cooling, because that needs
a run.

Reaches: `ion_coeff_by_ion_range(im_CI, ...)` -> `aion_m(:,im_CI)` ->
`met_b0(iel_C)` in every metals-on ionization solve, and
`ionization_balance_at_fixed_ne` / `normalized_reaction_residual`. Regression
cases `wasp_full`, `wasp_he23off`, `mol_metals`, `lower_profile` would move.

---

### Opacity model 'P': the pressure-broadening factor multiplies the column densities only, so the photons removed from the beam are not the photons that ionize

Status: CONFIRMED
Severity: P1 (opt-in: `opacity.inp` with `OPACITY_MODEL = P`)
Location: `src/modules/radiation/opacity_models.f90:251-262`;
`src/modules/functions/utilities.f90:461-466`, `:473`, `:498-500`, `:523-527`;
`src/modules/radiation/util_ion_eq.f90:770-772`, `:1002-1010`

```
utilities.f90:473	    dr = dr_j(j)*R0*opa_pf(j)
util_ion_eq.f90:770		int_f  = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
util_ion_eq.f90:771		int_1  = int_f*s_hi/e_v
util_ion_eq.f90:1002		acc_q = s_hi *nhi  (j)
```

**What the code does.** `opa_pf(j) = opacity_pT_factor(p)` is applied to the
radial step inside `calc_column_dens`, `calc_column_dens_one` and
`calc_column_dens_metals`, i.e. it scales `tauE`. The photoionization integrand
(`int_1`, `int_15`, `int_2`, `int_TR`, `int_h2`, `int_m`), the photoheating
integrand (`acc_HI`, `acc_HeI`, ...) and the absorbed-energy integrand (`acc_q`)
all use the unscaled `s_*` / `sigma_tab`.

**What it should do.** A cross-section multiplier enters the optical depth and
the absorption rate identically; a photon taken out of the beam has to ionize
something and deposit its energy. As written, model 'P' removes
`f(p)` times more photons from the beam than it converts into ionizations and
heat, so with `f(p) = 10` at the base nine tenths of the absorbed XUV energy
disappears from the ledger, and the heating efficiency `q = heat/q_abs` is
computed against the unscaled `q_abs` as well.

**How it was verified.** Read end to end: `opacity_pT_factor` has exactly two
call sites (`ionization_equilibrium.f90:782`, `post_process_adv.f90:373`), both
of which only fill `opa_pf`; `opa_pf` is read only by the three column routines
and by two diagnostic column dumps in `write_output.f90`. There is no site where
it multiplies a cross section.

Reaches: any run with `opacity.inp` selecting `OPACITY_MODEL = P`. No regression
case ships an `opacity.inp`, so no golden moves.

---

### The photon field of a cell is the field at its inner face: the column density includes the whole emitting cell

Status: CONFIRMED
Severity: P2 (bounded first-order-in-dtau bias; inherited verbatim from ATES-Code-main)
Location: `src/modules/functions/utilities.f90:461-483` (and `:498-500`, `:523-527`)

```
461	N1(N+Ng)  = dr_j(N+Ng)*R0*nhi(N+Ng)*opa_pf(N+Ng)
473	    dr = dr_j(j)*R0*opa_pf(j)
476	    N1(j)  = N1(j+1)  + nhi(j)*dr         ! HI
```

**What the code does.** `N1(j)` contains the whole of cell j, so
`exp(-tauE)` in `PH_heat_*` is the transmission down to the INNER face of cell j.
Every cell is therefore illuminated by a field attenuated by its own full optical
depth rather than by half of it (or, exactly, by the cell average
`(e^{-tau_out} - e^{-tau_in})/dtau`).

**What it should do.** The column to the cell CENTRE, or the exact cell-averaged
attenuation factor. The present form is a one-sided underestimate of the rate and
the heating in every cell.

**How it was verified.** Measured on the `wasp_full` golden
(`Ion_species.txt` + `Hydro_ioniz.txt`, R_p = 1.543e10 cm): the 13.6 eV
line-centre optical depth of one cell reaches 54 at the base and is 3.1e-2 at the
cell where the cumulative depth passes 1 (r = 1.3875 R_p). So the bias is ~1.5%
at the ionization front and negligible in the thin wind; deep inside the base it
is exponentially large in relative terms but the absolute rate is zero there. The
practical consequence is a front position shifted by up to half a cell.

Documentation note: this is identical to
`ATES/ATES-Code-main/src/modules/functions/utilities.f90:83,98`, so it is
inherited and not introduced by EXHALE.

---

### eval_cool applies the He I ground-state collisional coefficients to the total neutral helium and then adds the 2^3S terms on top

Status: CONFIRMED
Severity: P2 (measured effect <= 2.2e-4 of the He I channel)
Location: `src/modules/radiation/util_ion_eq.f90:1515-1519`, `:1568-1592`

```
1515	coio(j_lo:j_hi) =  2.179e-11*a_ion_HI(...)*nhi(...)      & ! HI
1516		  + 3.940e-11*a_ion_HeI(...)*nhei(...)               & ! HeI
...
1518	if (present(nheiTR)) coio(...) = coio(...) + (e_th_HeTR/erg2eV)*aion_HeITR(...)*nheiTR(...)
...
1568	coex(j_lo:j_hi) = coeff_coex_rate_HI(...)*nhi(...) + coeff_coex_rate_HeI(...)*nhei(...) + ...
1591		coex(...) = coex(...) + ( coeff_coex_HeI23S_10830(...) + ... )*nheiTR(...)
```

**What the code does.** `nhei` is the state vector's He I column, which
`composition.f90` documents as carrying the 2^3S population inside it
(`he_ground_singlet_density` is "the only place the difference is taken"). So the
metastable is charged 24.6 eV of collisional-ionization cooling and the
ground-state Cen-1992 collisional-excitation coefficient, and then its own 4.8 eV
ionization and 10830 A / singlet-conversion terms are added.

**What it should do.** Use `he_ground_singlet_density(nhei, nheiTR)`, exactly as
`PH_heat_HHe:706` and `photoheating_of_composition:1163,1167` do for the same
species. The comment at `:1583-1585` states the intent ("that channel is already
carried by the Cen-1992 He I coex term above"), which is only true if `nhei` there
is the ground singlet.

**How it was verified.** Read the state-vector convention in
`composition.f90:47-54` (the `he_ground_singlet_density` interface block); measured the size on the `wasp_full` golden:
max n(2^3S)/n(He I) = 2.157e-4 at r = 1.610 R_p, 3.9e-12 at the base. So the
double count is below the 0.1% movement threshold everywhere and no golden moves.

---

### He II recombination cooling is charged at alpha_B while the triplet path recombines He+ at a different total coefficient, and the comment that bounds the mismatch does not describe that path

Status: CONFIRMED
Severity: P2
Location: `src/modules/radiation/ionization_equilibrium.f90:935-938`;
`src/modules/radiation/Cool_coeff.f90:3409-3413` (`lambda_rec_HeII`)

```
ioniz_eq:937	! The He II recombination *cooling* (rec_cool_HeII = kT alpha_B) is left
ioniz_eq:938	! unchanged; the mismatch is <= y alpha_1 kT, negligible.
```

**What the code does.** `lambda_rec_HeII(T) = kb_erg*T*alpha_rec_HeII_B(T)`, i.e.
kT times the case-B coefficient, always. With the triplet tracked,
`HeITR_coeffs` overwrites `rcheiiB` with `rec_HeII_11S` and the balance also
carries `rcheiTR = rec_HeII_23S`; with `use_he_rec_coupling` on (default),
`he_rec_coupling` further replaces `rcheiiB` by `y_net*alpha_1 + 0.25*alpha_B`.

**What it should do.** The electron gas should lose kT per recombination the
solver actually performs. The quoted bound `<= y alpha_1 kT` describes only the
atomic branch, where `rcheiiB_new = alpha_B + y_net*alpha_1`.

**How it was verified.** Hand-evaluated at T = 1e4 K from the coefficients in
this file: alpha_B(He II) = 2.81e-13, Mao alpha_1 = 1.57e-13,
`rec_HeII_11S` = 1.54e-13, `rec_HeII_23S` = 2.10e-13. In the triplet branch the
solver removes He+ at 3.64e-13 (no coupling) or ~4.3e-13 (with the coupling and
y_net ~ 1) cm^3 s^-1 while the cooling is charged at 2.81e-13 -- 23% to 35% low
on that channel. He II recombination cooling is a small share of `reco` at
He/H = 0.08, so the effect on total cooling is well under a percent; I did not
measure it in a run.

Reaches: every `Include He23S? True` run (the whole regression matrix except
`wasp_he23off`).

---

### Chemical-equilibrium retry seed applies the (1 - x_H2) factor twice to the proton fraction

Status: CONFIRMED
Severity: P2 (a starting point only; it cannot change a converged root)
Location: `src/modules/radiation/ionization_equilibrium.f90:3143`, `:3156`

```
3143	x(1) = x(1)*(1.0d0 - x_h2)
...
3156		nhis = max(1.0d0 - x_h2, 0.0d0)*(1.0d0 - x(1))*ieq_cell%nh
```

**What the code does.** After line 3143, `x(1)` is already the H+ fraction PER H
NUCLEUS. The neutral atomic H density is therefore `n_H (1 - x_H2 - x(1))`, but
line 3156 forms `n_H (1 - x_H2)(1 - x(1))`, which applies `(1 - x_H2)` a second
time inside the bracket.

**What it should do.** `nhis = max(1 - x_h2 - x(1), 0) * ieq_cell%nh`.

**How it was verified.** Read; `nhis` is the atomic-H density handed to
`oxygen_chemical_equilibrium_fractions`, which sets the OH/H2O seed of attempt 2
of the molecular cell solve. Since the hybrd1 root is what is accepted, this
changes only which basin the retry starts in (and only with the oxygen chemistry
on, which no regression case uses).

---

### The pure-hydrogen secondary-ionization path gives helium's share of the ionization budget to hydrogen

Status: CONFIRMED
Severity: P2 (documented as a design choice, but the size is not stated; `thereis_He = .false.` is rare, and the same limit is reached in the He path wherever He I is fully ionized)
Location: `src/modules/radiation/electron_energy_degradation.f90:663-670`;
call site `src/modules/radiation/util_ion_eq.f90:486-489`

**What the code does.** `w_tot = n_Hw*heh_neutral_svs85*f_HI + n_HeI*f_HeI`, and
`c_ion_HI = f_ion * heh_neutral_svs85 * f_HI / w_tot` with
`f_ion = f_HI + f_HeI`. With `n_HeI = 0` this gives `c_ion_HI = f_ion/n_HI`, i.e.
the whole Shull & van Steenberg ionization budget, He channel included, drives
hydrogen secondary ionizations.

**What it should be measured against.** At x = 0 the amplitudes are
`f_HI = 0.3908` and `f_HeI = 0.0554`, so the H I secondary-ionization rate is
`(0.3908 + 0.0554)/0.3908 = 1.142` times the SvS85 H I channel. The energy is
conserved by construction (the module proves the identity), so the 14% comes out
of what would otherwise have been excitation/heat; `f_heat` is clamped against
the same enlarged `f_ion`.

**How it was verified.** Algebra on the routine plus the limits the module's own
header states. Not measured in a run.

---

### Dead and stale items

Status: CONFIRMED
Severity: P2

1. `src/modules/radiation/util_ion_eq.f90:88-90` -- `f_sing_HeI` and
   `Ee_sing_HeI` are declared as parameters and referenced only from comments
   (`:2684`, `:2736`); `he_rec_coupling` rebuilds both quantities inline
   (`fs_HI`, `Ees_HI`). Dead parameters.
2. `src/modules/functions/h2_photo_channels.f90:115` -- `h2_channel_threshold` is
   public and has no consumer anywhere in `src/`.
3. `src/modules/functions/h2_photo_channels.f90:299-303` -- the header of
   `h2_channel_stoichiometry` says "This is the single table every source term is
   derived from; nothing downstream should write these numbers out by hand." The
   source terms are in fact written out by hand in
   `System_HeH_mol.f90:245-246,271,279-280` (`P_H2_di + 2 P_H2_dd`, etc.); the
   only caller of `h2_channel_stoichiometry` is `src/tests/e1_h2/`.
4. `src/modules/radiation/Cool_coeff.f90:2184` -- the
   `coronal_excitation_cutoff` header says the softest O I exponential is
   `exp(-930.111/T)`; `cool_OI_chianti` (`:1252-1259`) carries `exp(-847.937/T)`.
   The C I number quoted in the same sentence (2351.38) does match.
5. `src/modules/radiation/util_ion_eq.f90:1712-1716` -- "the iscool ions, in
   canonical order, are exactly CI,CII,OI,OII,NI,NII,MgI,MgII".
   `species_table.f90:235-240` also flags Ca II, Na I, Fe I and Fe II as
   `mion_iscool`, and the loop directly below the comment sums all twelve.
6. `src/modules/radiation/Cool_coeff.f90:667-681` -- `gbar_ff` uses
   `Ry/k = 157807 K`. van Hoof et al. (2014) define gamma^2 with the Rydberg
   ENERGY, 157887.5 K; 157807 K is the Hui & Gnedin (1997) H I ionization
   temperature this module uses elsewhere. The resulting shift in
   `log10(gamma^2)` is 2.2e-4 dex against a 0.1 dex table step -- numerically
   irrelevant, but the comment names the wrong constant.
7. `src/modules/radiation/util_ion_eq.f90:1667-1711` -- the legacy
   `use_2lev_cool` branch sums C I, C II, O I, O II, Mg I, Mg II, Ca II (17),
   Na I (19) and Fe II (26) but omits Fe I (25), whose `c_metal(:,im_FeI)` was
   computed a few lines above and is a non-zero coolant. N I / N II are also
   omitted, but they are zero in that branch by construction. Opt-in comparison
   branch only.
8. `src/modules/files_IO/metals_input_read.f90:93-176` -- `n_set` is incremented
   only for abundance lines; the six option keys `cycle` before it, so the
   "N metal abundance(s) set" line under-reports what the file changed.

---

### Two different spectra describe 3.4-13.6 eV inside one run, and the Balmer continuum carries no attenuation

Status: CONFIRMED
Severity: P2 (an approximation the source, Christie et al. 2013, also makes; recorded because it interacts with the sub-Lyman finding above)
Location: `src/modules/radiation/excited_hydrogen.f90:132-133`, `:284-345`

**What the code does.** `gamma2_bal` and `hpe2_bal` are SCALARS: the n=2
photoionization rate and its photoelectric heating per atom are computed once per
update from `pi B_nu(T_star_eff)(R_star/a)^2` over 3.4-13.6 eV, with no radial
attenuation. Over the overlapping part of that band (`e_sub_low`..13.6 eV) the
XUV grid simultaneously carries the extrapolated EUV power law, which is what
He 2^3S and the low-IP metals absorb (and which IS attenuated).

**What it should do.** One band, one field. Christie et al. (2013) footnote 5
justifies dropping the attenuation for the n=2 columns specifically ("minimal
attenuation of Gamma_2s and Gamma_2p due to the small columns of n = 2
hydrogen"), which does not extend to the low-IP metals, whose columns are not
small.

**How it was verified.** Read both paths; the flux ratio at 4.8-13 eV is the
table in the second finding above. `use_excited_H` is set whenever the input
states `Stellar Teff` and `Stellar radius`, which `wasp_full` and every planet
example do, so both descriptions are live simultaneously in the default runs.

---

## Checks that came out clean

Verified against the published source and found correct, so that the report says
what was checked and not only what failed:

* **Verner et al. (1996) form.** `sigma_VFKY96` implements Eq. (1) with
  `x = E/E_0 - y_0`, `z = sqrt(x^2 + y_1^2)`, `Q = 5.5 - 0.5P` (no orbital-l
  term); the high-energy limit comes out as E^-3.5, which is the correct
  asymptote. `sigma_Verner96` is the VY95 form with `Q = 5.5 + l - 0.5P`, and the
  two agree when `y_0 = y_1 = 0, l = 0`, as the header claims. Argument order
  matches the declaration for all 17 metal ions.
* **Badnell (2006) RR fits** for H II, He II and He III checked line by line
  against `references/Badnell_2006_ApJS_167_334.pdf` (rows `1 0 8.318(11) 0.7472
  2.965(0) 7.001(5)`, `2 0 1.818(10) 0.7492 1.017(1) 2.786(6)`, `2 1 5.235(11)
  0.6988 7.301(0) 4.475(6) 0.0829 1.682(5)`): exact.
* **Voronov (1997)** rows for H I, He I, He II, C II, N I, N II, O I, O II,
  Mg I, Mg II, Si I, Si II, Ca I, Ca II, Na I, K I, S I, Fe I, Fe II: all correct
  (C I is the single exception reported above).
* **Christie, Arras & Li (2013) Table 2** checked against
  `references/Christie_2013ApJ_772_144.pdf`: R2 alpha_B, R3 C_1s->2s, R4
  C_1s->2p, R5 C_2s->2p (including the `sqrt(T)` and `T0 = 1.02 K`), R8/R9 the
  2s/2p split, R10/R11 the A values -- all match `hydrogen_n2_rates.f90`.
  The 2x2 solve in `n2_populations` is the correct Cramer solution of the
  documented system, and the statistical weights are module parameters that the
  dummy names can no longer shadow.
* **Line-centre opacity** `line_center_opacity_lte`:
  `lambda^3/(8 pi^{3/2}) (g_u/g_l) A n_l (1 - e^{-Ek/T})/v_th` is the standard
  Doppler-core result (derived independently); `hc/k = 1.43877736 cm K` correct.
* **`C_lya = 1.49736e-2 f_lya`** reproduces `sqrt(pi) e^2/(m_e c) = 1.49733e-2`.
* **Neufeld/Harrington wing escape.** `beta = pi^{-1/4} sqrt(a/tau)` is exactly
  half the single-flight escape probability
  `2 pi^{-1/4} sqrt(a/tau)` obtained by integrating the damping wing beyond
  `x_1 = sqrt(a tau/sqrt(pi))` (derived independently), i.e. escape through one
  face; `x_1` in the stellar-penetration term uses the same definition, so the
  two are consistent. `line_escape_probability_one_face` is continuous at its
  branch point `tau_c = sqrt(pi) e^{a^2/4} = 6.967` (both branches give 0.03067).
* **Three-level fine-structure equilibrium.** `fine_structure_populations_3level`
  is the correct Cramer solution; the level ordering for the INVERTED O I term
  (level 1 = 3P2, A31 = A(44 um), A32 = A(145 um)) is consistent between
  `oxygen_ground_term_collisions`, `cool_OI_ne_func`,
  `oxygen_ground_term_populations` and the `beta_fs`/`nbar_fs` slot order in
  `cool_OI_ne_range` and `oxygen_ground_term_levels`; the `g_u` used in
  `electron_impact_deexcitation` is the upper level of each pair in all four ions.
  `fine_structure_line_opacity` uses the right lower level, partition function
  and `g_u/g_l` for all eight lines.
* **SED unit conversion.** `F_XUV = F_lambda * hp_eV*c*1e8/E^2` is
  `F_lambda |d lambda/dE|` with lambda in Angstrom and E in eV; `de_v` in
  `read_sed` is positive because the file is ordered in decreasing energy; the
  `LEUV`/`LX` split loop exits correctly on that ordering; `sum(de_v)` telescopes
  to `e_v(Nl) - e_v(1)` exactly.
* **Secondary ionization.** The share identity
  `n_HI c_HI + n_H2 c_H2 + n_HeI c_HeI = f_ion` holds identically at every
  tabulated energy; each absorber block divides its own excess energy
  `e_v - e_th(absorber)` by the TARGET's threshold (13.6 / 24.6 / 15.4), which is
  the correct pairing; `f_heat` is clamped to `1 - f_ion`; the H2 dissociative
  secondary yield `dissoc_ion_per_H2p = 1/22` adds to (not shares) `Psec_H2`, as
  the Dalgarno harmonic-mean statement requires; the SvS85 amplitudes 0.3908 /
  0.4092 / 1.7592 and 0.0554 / 0.4614 / 1.6660 are correct; the
  `a(x) = 0.5 (x 1e4)^0.15` branch of `h2_vibrational_share` joins its neighbours
  exactly at x = 1e-4 (0.5) and x = 1e-7 (0.1774), which is the typographical
  correction the header claims.
* **Bremsstrahlung** uses the ion NET charge (`Z_ion^2 gbar(Z_ion)`), He III with
  the factor 4 and the Z=2 Gaunt table, molecular ions deliberately excluded and
  documented; `1.426e-27 sqrt(T)` is the standard prefactor.
* **Charge exchange.** All 64 descriptor arrays have consistent lengths (23 A + 2
  B + 6 C + 32 D + 1 E) and the donor/acceptor stage pairs decode to the labelled
  reactions; `cx_fvidx` maps to the right residual rows for every active reaction
  (no active row addresses the pinned X++ identity row of a two-stage element);
  `he_h_cx_jac` is the exact derivative of `he_h_cx_fvec`; `he_h_cx_fvec_adv` correctly divides each rate by the nucleus count of its own row. The O I / H+ pair reassignment satisfies
  detailed balance to 9.4% at 1e4 K
  (k_f/k_r = 0.788 against (8/9) exp(-227/T) = 0.869), as its comment states.
* **H2 channel split.** `d = 2 rho (1 + r_pm)/(r_pm (1 + 2 rho))` is the correct
  inversion of `rho = (d r_pm/2)/(1 + r_pm(1-d))`; `sum(sig) = sigma_tot` by
  construction; the `off` default reduces to the old `f_di` branching exactly.
* **Metal recombination dispatch.** All seventeen `rec_*_range` wrappers request
  the correct daughter label (`rec_CII -> 'CI'`, etc.).
* **OpenMP.** `PH_heat_HHe`'s `PRIVATE` list covers every scratch vector and
  scalar written in the loop body; the shared output arrays are written only at
  the loop index. `eval_cool`'s block decomposition is disjoint and its
  `reduction(+:ec_t)` targets a module `save` array used for timing only.
  `cx_kc` is threadprivate and lazily allocated on first use inside each thread;
  `cx_metal_base` is `copyin`; `ieq_cell` is threadprivate
  (`ion_cell_state.f90:135`); `sys_x`/`sys_sol`/`wa` are lazily allocated per
  thread. `count == 0` sweeps run serially because of the neighbour warm start.
  I found no uninitialized threadprivate read.
* **Numerical hygiene.** `absorbed_photon_fraction`'s series is the correct
  truncation of `1 - e^{-tau}`; `cool_table_value` and `cool_table_value_2d` use
  `.not.(pos > 1)` so a NaN temperature lands on the table edge;
  `residual_decade` brackets a non-finite residual into the top bin;
  `coronal_excitation_cutoff` is C1 at the floor; `he_rec_coupling`'s `1/nhig`
  divisions cancel against the `u2`/`um` ratios in the n_HI -> 0 limit.
* **`dayside_dilution()`** returns 0.25 / 0.5 / 1 for Rate/4 / Rate/2 / alpha and
  is applied once to the XUV grid, the five FUV bands, the stellar Ly-alpha beam
  and the Balmer continuum; I found no path that applies it twice or omits it.

---

Files read completely: `src/modules/radiation/util_ion_eq.f90`,
`src/modules/radiation/ionization_equilibrium.f90`,
`src/modules/radiation/Cool_coeff.f90` (all executable statements and comments;
the numeric bodies of `cool_logL_FeII_ne`, `cool_logL_FeI` and `gff_avg` were
inspected for shape and endpoints, not digit by digit),
`src/modules/radiation/electron_energy_degradation.f90`,
`src/modules/radiation/excited_hydrogen.f90`,
`src/modules/radiation/hydrogen_n2_rates.f90`,
`src/modules/radiation/lya_rt.f90`,
`src/modules/radiation/charge_exchange.f90`,
`src/modules/radiation/sed_read.f90`,
`src/modules/radiation/opacity_models.f90`,
`src/modules/radiation/J_inc.f90`,
`src/modules/functions/cross_sec.f90`,
`src/modules/functions/h2_photo_channels.f90`,
`src/modules/init/set_energy_vectors.f90`,
`src/modules/files_IO/opacity_input_read.f90`,
`src/modules/files_IO/metals_input_read.f90`.
Read in part, as callers: `src/modules/functions/utilities.f90`
(`calc_column_dens*`), `src/modules/functions/composition.f90`
(`he_ground_singlet_density`), `src/modules/init/parameters.f90` (constants,
defaults, `dayside_dilution`), `src/modules/files_IO/input_read.f90` (spectrum
and energy-band parsing, `use_excited_H`), `src/modules/init/species_table.f90`
(ion order, `mion_iscool`), `src/modules/nonlinear_system_solver/T_equation.f90`
(metal-cooling dispatch), `src/modules/nonlinear_system_solver/System_HeH_mol.f90`
(H2 photo source rows), `src/modules/files_IO/write_setup_report.f90`.

Not checked:
* The numeric content of the large tables was not re-derived: `cool_logL_FeII_ne`
  (41x29), `cool_logL_FeI` (41), `gff_avg` (161), the Samson & Haddad and Backx
  `sigma_H2` tables (72 + 7 rows), the Chung `f_di` table (69 rows), the Barragan
  O2+ + H table (26 rows), the Dalgarno Tables 4/5/7 coefficient blocks, and the
  CHIANTI fine-structure Upsilon / k_H polynomial coefficients. I checked their
  form, ranges, endpoints and units, not every digit.
* The Verner et al. (1996) metal coefficient sets were checked for internal
  consistency (argument order, threshold, form) but not against the published
  Table 1, which is not in `references/`.
* `rr_mao` (Mao & Kaastra 2016) coefficients: the functional form is as the code
  documents it, but the paper is not in `references/` and the seven coefficients
  per ion were not verified.
* The Shull & Van Steenberg (1982) Ca I recombination power law
  `6.78e-13 (T/1e4)^-0.8` was not verified against a source.
* The Oklopcic-style `rec_HeII_23S` / `rec_HeII_11S` power laws were not checked
  against their source; only their consistency with the recombination-cooling
  coefficient was (see the He II finding).
* No run was made, so no finding carries a measured change in a converged wind:
  the sizes quoted are of the coefficients, integrands and rates themselves.
* The molecular network residual rows, the H3+ / molecular infrared modules, the
  Lyman-Werner and water-photolysis modules, and the steady solver were outside
  this group and were read only where a caller relationship required it.

# EXHALE `input.inp` schema (authoritative)

This document is the authoritative schema of EXHALE's `input.inp` file. It began
as §5.6 Inc 0 of `docs/refactor_plan_system_composition_parser.md` ("document the
schema"), the prerequisite for the increments that converted the Fortran
positional reads to anchored label matching and unified the Python loaders.
Those increments have since been carried out; §3.0 describes the parser as it
now behaves, and §3.1--3.2 are retained as a record of the positional design it
replaced.

**How to find a key in the source.** Every key is located by its label, not
by a line number: `lbl_match(line, '<key>')` in
`src/modules/files_IO/input_read.f90` for the reader, and the literal the Tk
writer emits in `src/utils/EXHALE_interface_functions.py` for the GUI. This
document carries no line-number citations into either file: a label survives
the edits that move a line, and a citation that no longer points at the key it
names is worse than none.

The single source of truth for the file format is the Fortran parser
`src/modules/files_IO/input_read.f90` (`subroutine input_read`, plus its
`get_word` word-splitter). Everything below is derived from that code and
cross-checked against the Python readers, the GUI writer, and the shipped
example inputs. Where the manual and the code disagree, the code wins, and the
mismatch is called out in the "Discrepancies and fragilities" section.

## 1. What the file is and who reads it

`input.inp` is a plain-text file, read line by line. It has two structurally
different regions:

1. A **core block** (the planet/star parameters and the numerics selectors).
   Each of these keys is mandatory; a missing one aborts with `error stop 1`.
   They are located by anchored label match (`req` -> `find_lbl` ->
   `lbl_match`), so their order in the file does not matter, but the *value* is
   still taken by word position on the matched line via `get_word(line, n)`. A
   few are **conditional**: the spectrum-property line depends on
   `Spectrum type`, the energy-band line is skipped for a monochromatic
   spectrum, and the X-ray luminosity line is read only when X-rays are
   included. (Historically this block was read strictly by line order; §3.1
   records that design.)

2. An **optional keyword-extension block**. Every line of the file is scanned in
   a `do` loop and matched with the same anchored `lbl_match`. Any of these keys
   may be omitted (the code holds a default for each), may appear in any order,
   and blank / `#`-comment lines are skipped. An unrecognized non-blank line
   draws a warning. These keys carry the newer physics and solver options.

Four independent consumers read this file:

- **Fortran** (`input_read.f90`): the authoritative reader; drives the
  simulation.
- **`examples/exhale_io.py`** (`read_input`): a label-based loader used by the
  analysis notebooks. It splits each line on the first `:` and matches keys by
  string prefix, so it is robust to line order, but it surfaces only a subset of
  the parameters.
- **`EXHALE_transit.py`** (top level of `EXHALE/`) via
  `exhale_transit_lib.read_input_params`: the transit post-processor. Every
  field it needs (Rp, Mp, T0, a_orb, Mstar, LEUV, `2D approximate method`) is
  matched by label through `find_input_label`; the positional counter it used
  before Inc 2 is gone.
- **`src/utils/EXHALE_interface_functions.py`**: the Tk GUI writer
  (`start_func`). It writes the core block only (through `Force start`); it does
  not emit any keyword-block line.

Two companion configuration files are also read at startup and are covered in
the appendices: `metals.inp` (its mere presence turns metals on) and
`opacity.inp`. A third optional file, `base.inp`, is read by
`read_base_inp` and can override a few core values after `input.inp` is parsed;
it is summarized at the end of the main table notes.

Word positions below refer to whitespace-separated tokens counted from 1, as
`get_word` counts them. In the "Read by" column, "GUI" means the Tk writer in
`EXHALE_interface_functions.py` emits the line, "IO" means
`exhale_io.py:read_input` reads it, and "TR" means `EXHALE_transit.py` reads it
with its `num` counter. The Fortran parse site is found by the key's own label,
as the header states.

## 2. Main table

### 2a. Core block (fixed order; consumed first)

Rows are in the exact order `input_read.f90` reads them.

| # | Line / label | Match mode | Type | Units | Default (if absent) | Read by | Notes |
|---|---|---|---|---|---|---|---|
| 1 | `Planet name:` | positional, word 3 | string | - | mandatory | GUI; IO raw `Planet name` | Sets `p_name`. Value is a single token (word 3); a name with spaces would be truncated. |
| 2 | `Log10 lower boundary number density [cm^-3]:` | positional, word 7 | real | log10(cm^-3) | **optional since 2026-09-03** | GUI; IO `n0_log10` | `n0` is stored as the exponent, then raised to `10**n0`. One of the three inputs that can state the base level, and startup refuses two that disagree, see `p_base` in section 2c and "the base level has one source" in `input_read.f90`. Not read by TR. |
| 3 | `Planet radius [R_J]:` | positional, word 4 | real | R_Jupiter = 7.1492e9 cm | mandatory | GUI; IO `Rp_RJ`; TR num==2 word 4 | Converted to cm (`R0*RJ`). The unit is the IAU 2015 nominal equatorial radius, the one transit depths are quoted against (`parameters.f90` `RJ`). |
| 4 | `Planet mass [M_J]:` | positional, word 4 | real | M_Jupiter = 1.8982e30 g | mandatory | GUI; IO `Mp_MJ`; TR num==3 word 4 | Converted to g (`Mp*MJ`). The unit is the IAU 2015 nominal mass (`parameters.f90` `MJ`). |
| 5 | `Equilibrium temperature [K]:` | positional, word 4 | real | K | mandatory | GUI; IO `T0`; TR num==4 word 4 | Base temperature `T0`. Can be overridden later by `base.inp`. |
| 6 | `Orbital distance [AU]:` | positional, word 4 | real | AU | mandatory | GUI; IO `a_AU`; TR num==5 word 4 | Converted to cm (`a_orb*AU`). |
| 7 | `Escape radius [R_p]:` | positional, word 4 | real | R_planet | mandatory | GUI; IO `r_esc` | Inner escape radius. Not read by TR. |
| 8 | `He/H number ratio:` | positional, word 4 | real | - | mandatory | GUI; IO `HeH` | `HeH > 0` sets `thereis_He`. Can be overridden by `base.inp` and zeroed by a sub-He monochromatic spectrum. Not read by TR. |
| 9 | `2D approximate method:` | positional, word 4 (+ word 6 if `alpha`) | keyword string | - | mandatory | GUI; IO raw `2D approximate method`; TR content match | Word 4 is one of `Mdot/4`, `Rate/2`, `Rate/4`, `alpha`. `Rate/2` is rewritten to `Rate/2 + Mdot/2`, `Rate/4` to `Rate/4 + Mdot`. If word 4 is `alpha`, word 6 is read into `a_tau`. Sets the day-night flux factor. |
| 10 | `Parent star mass [M_sun]:` | positional, word 5 | real | M_sun = 1.98842e33 g | mandatory | GUI; IO `Mstar_Msun`; TR num==9 word 5 | Converted to g (`Mstar*Msun`). The unit is the IAU 2015 nominal mass (`parameters.f90` `Msun`). |
| 11 | `Spectrum type:` | positional, word 3 | string | - | mandatory | GUI; IO raw `Spectrum type` | Word 3 must be `Load`, `Power-law`, `Planck`, or `Monochromatic`; anything else is a fatal `error stop 1`. The GUI writes `Load from file..`, whose word 3 is `Load`. Selects the conditional line 11a/11b/11c/11d below. **The type states the spectrum of the WHOLE photon grid**, the XUV and the band below 13.6 eV alike, where the He 2^3S metastable (4.80 eV) and the low-IP metals absorb: no band is built from a type the input did not select (`docs/development_plan_20260905_rev3.md` section 10.5 decision 13). `write_setup_report` states the type, the source of the sub-13.6 eV band, and the integrated grid flux against the nominal `(10^LX + 10^LEUV)/(4 pi a^2)` for every run. |
| 11a | `Spectrum file:` (only if `Load`) | positional, word 3 | string (path) | - | conditional | GUI; `sed_read.f90` (`read_sed`) | Sets `sed_file`, `do_read_sed = .true.`. **The table must reach the floor of the photon grid.** That floor is the lowest ionization threshold of any active absorber: 4.80 eV = 2583.0 A with `Include He23S? True`, or the lowest active neutral-metal threshold when a low-IP metal is in `metals.inp` (K I at 4.341 eV = 2856.1 A is the lowest the species table has; the shipped metal sets bottom out at Na I, 5.139 eV = 2412.6 A), and 13.6 eV when neither is active. A table that stops above it leaves that absorber with no field over part of its own band, so `read_sed` stops the run with a fatal `error stop 1`: the message names the file, the band not covered in eV and in A, each absorber left without a field with its threshold, and the two remedies (drop the absorber, i.e. `Include He23S? False` or remove the element from `metals.inp`; or supply an SED reaching the stated wavelength). There is no key to continue (`docs/development_plan_20260905_rev3.md` section 10.5 decision 17). The inventory of SED files and which of them cover which floor is `inputdata/sed/README.md`. |
| 11b | `Power-law index:` (only if `Power-law`) | positional, word 3 | real | - | conditional | GUI | Sets `PLind`, `is_PL_sed = .true.`. The law is normalized on `[e_low, e_mid]` and evaluated wherever the grid reaches, so below `e_low` it is an extrapolation of an EUV fit -- that is what the type means. Measured on `backup/regression/wasp_full` (index -1): with the grid floor at 4.80 eV (He 2^3S on) the grid integrates to 1.420 times the nominal XUV flux, and with the floor at 13.6 eV to 0.995, so the extrapolated sub-Lyman band carries 0.425 of the nominal flux. |
| 11d | `Stellar Teff [K]:` + `Stellar radius [R_sun]:` (only if `Planck`) | keyword block (K3, K4) | reals | K, R_sun | conditional | - | The photospheric blackbody `pi B_nu(T_eff) (R_star/a)^2` per unit photon energy on every point of the grid (`planck_stellar_flux_eV`, `J_inc.f90`), diluted like every other stellar beam. No property line of its own: it reads the two optional stellar keys, and a `Planck` run missing either of them is a fatal `error stop 1`. Its absolute scale comes from `T_eff` and `R_star`, not from `LX`/`LEUV`, so the integrated grid flux differs from the nominal XUV flux (113x on the `wasp_full` parameters) and the setup report prints both. A real photosphere is not a blackbody -- line blanketing and the Balmer jump depress the near-ultraviolet -- so this field is an upper bound below 13.6 eV and falls far below the true EUV of an active star above it; a run that needs both bands right needs `Load`. |
| 11c | `Photon energy [eV]:` (only if `Monochromatic`) | positional, word 4 | real | eV | conditional | GUI | Sets `e_low`, `is_monochr = .true.`. If below the He I threshold, He is removed and `HeH` is zeroed. |
| 12 | `Use only EUV?` | positional, word 4 | bool-ish | - | mandatory | GUI | Word 4 `== 'False'` sets `thereis_Xray = .true.` (X-rays INCLUDED). `True` leaves X-rays off. See the double-negative note in §6.7. |
| 13 | `[E_low,E_mid(,E_high)] = [ ... ]` (skipped if monochromatic) | positional, words 4, 6, 8 | reals | eV | conditional | GUI | If X-rays off: read `e_low` (word 4) and `e_mid` (word 6); `e_top` defaults to `1.24e3`. If X-rays on: also read `e_top` (word 8). The GUI writes `-`-separated values; examples use `,`; both parse because only word positions 4/6/8 matter. |
| 14 | `Log10 of X-ray luminosity [erg/s]:` (only if X-rays on) | positional, word 6 | real | log10(erg/s) | conditional | GUI | Sets `LX`; `LX = 0` when X-rays off. |
| 15 | `Log10 of EUV luminosity [erg/s]:` | positional, word 6 | real | log10(erg/s) | mandatory | GUI; IO `LEUV`; TR content match | Sets `LEUV`. |
| 16 | `Grid type:` | positional, word 3 | string | - | mandatory | GUI | `Uniform`, `Mixed`, or `Stretched`. |
| 17 | `Numerical flux:` | positional, word 3 | string | - | mandatory | GUI | `HLLC`, `ROE`, or `LLF`. |
| 18 | `Reconstruction scheme:` | positional, word 3 | string | - | mandatory | GUI | `PLM`, `WENO3`, or `PLM+WENO3` (two-stage, `recon_two_stage`). The GUI list offers only `PLM`/`WENO3`; `PLM+WENO3` is hand-added. This choice governs how the later `du_th` values are used. |
| 19 | `Include He23S?` | positional, word 3 | bool-ish | - | **optional**, default on | GUI | The one core-block line that may be omitted: `thereis_HeITR` defaults to `.true.` in `parameters.f90`, so a file without the line gets the triplet. Word 3 `True`/`true` sets it, `False`/`false` clears it (the deliberate opt-out, which the HeITR-off regression case uses); any other word leaves the default. Forced off later if He is absent. |
| 20 | `Load IC?` | positional, word 3 | bool-ish | - | mandatory | GUI | `True` sets `do_load_IC`. A restart restores every species density the file carries, with one exception: an element whose reservoir a lower-atmosphere handoff states (`<El>_H_base`, or an elemental ratio of the `Lower atmosphere profile:` file) has its loaded column renormalized onto that `El/H` by a single factor common to all its ionization stages, reported in the run log. Elements the handoff does not state are loaded unchanged. |
| 21 | `Do only PP:` | positional, word 4 | bool-ish | - | mandatory | GUI | `True` sets `do_only_pp` and clears `force_start`, and arms `sec_ion_active` before the loop so the single equilibrium solve uses the same physics as the post-process (`Update_EXHALE_stage1.md` sec. 94). |
| 22 | `Force start:` | positional, word 3 | bool-ish | - | mandatory | GUI | `True` sets `force_start` and clears `do_only_pp`. Historically the last positional line; with anchored matching its position no longer matters. |

### 2b. Keyword-extension block (optional; scanned after the core block)

Every line is matched by `lbl_match(line, 'KEY')` in the keyword `do` loop, in
the order shown. Any line may be omitted. None of these are read
by the Python loaders except where noted.

| # | Key substring | Match / value word | Type | Units | Default (if absent) | Sets | Notes |
|---|---|---|---|---|---|---|---|
| K1 | `Domain mode` | word 3 == `Spherical` | flag | - | `.false.` (Roche) | `spherical_domain` | Spherical potential out to `Outer radius`. |
| K2 | `Outer radius` | word 4 | real | R_planet | `0.0` | `r_out_user` | Required (> 1) in spherical mode; in Roche mode extends `r_max` past L1. |
| K3 | `Stellar Teff` | word 4 | real | K | `0.0` | `T_star_eff` | With K4 enables the diluted-blackbody Balmer continuum (`use_excited_H`). **Mandatory, and > 0, with `Spectrum type: Planck`** (line 11d): it is then the source of the whole photon grid as well. |
| K4 | `Stellar radius` | word 4 | real | R_sun | `0.0` | `R_star` | Converted to cm (`R_star*Rsun`). **Mandatory, and > 0, with `Spectrum type: Planck`** (line 11d). |
| K5 | `Deexc heat` | word 3 == `True`/`False` | flag | - | `.true.` | `incl_deexc_heat` | Collisional de-excitation heating of H(n=2); active only with the excited-H model (K3+K4). `False` reverts to the one-way coronal ledger. |
| K6 | `Wind-AE seed out` | word 4 | string (path) | - | `''` | `windae_seed_out` | Must be tested before K7 (substring). |
| K7 | `Wind-AE seed` | word 3 | string (path) | - | `inputdata/windae_seed.csv` | `windae_seed_file` | Path of the CSV warm-start seed the `IC mode: windae` generator loads. The keywords `grid` and `auto` instead pick the nearest solution in `inputdata/windae_grid/` to this planet's (Mp, Rp, XUV flux) (`wae_exhale_bridge.f90`). |
| K8 | `Jlya RT file` | word 4 | string (path) | - | `jlya_rt.txt` | `jlya_rt_file`, `jlya_mode=1` | Read a J_Lya(r) profile. |
| K9 | `Jlya escape-prob` | word 3 == `True` | flag | - | `jlya_mode=0` | `jlya_mode=2` | Requires `Stellar Lya flux > 0` or fatal `error stop`. |
| K10 | `Stellar Lya flux` | word 5 | real | erg/cm^2/s | `0.0` | `F_Lya_star` | With K15c on it is also band `B2` of the oxygen photolysis (the 1215.67 A line), so a value set for the H(n=2) pumping drives the photolysis as well. |
| K11 | `Lya stellar halfwidth` | word 5 | real | km/s | `70.0` | `dv_star_lya` | Half-width of the stellar Ly-alpha line, which sets how far into the line wings the beam penetrates (broad plateau of about +-70 km/s, Huang et al. 2017). |
| K12 | `Lya stellar boost` | word 5 | real | - | `5.0` | `lya_star_boost` | Bounded buildup of the trapped stellar beam: the mean intensity reaches `1 + (boost - 1)(1 - beta) T_star` times free streaming, so this value where the beam is fully trapped and 1 in the thin outer wind. Tuned to Huang et al. (2023) Fig. 11. |
| K12b | `Lya absorbing bottom` | word 4 == `True`/`False` | flag | - | `.false.` (reflecting bottom) | `lya_bottom_absorber` | Closes the Ly-alpha domain from below with a pure absorber instead of a reflector, as Huang et al. (2017) do: the H2 layer beneath the base really does absorb Ly-alpha through accidental resonances, and the default reflecting bottom keeps every downward photon and over-fills `Jbar` near the base. Puts the planet-ward face of the trapping slab at the bottom of the domain, at line-centre depth `tau(base) - tau(r)`, instead of mirroring the star-ward one: it shortens the frequency random walk of the published slab solution rather than adding a factor, since the two faces share one escaping population (Neufeld 1990 eq. 2.25). Acts on the in-line escape-probability RT only (`jlya_mode = 2`, i.e. `Jlya escape-prob: True`); `input_read` warns and ignores it for `jlya_mode` 0 and 1. See `src/modules/radiation/lya_rt.f90`, `docs/lya_destruction_channels.md`. |
| K13 | `du_th` | word 3 (+ optional word 4) | real(s) | - | `du_th=1.0e-3`, `du_th_plm=-1` | `du_th`, `du_th_plm` | Two values only if `Reconstruction scheme: PLM+WENO3` (order dependency on line 18). |
| K13b | `Reconstruction continuation` | word 3 (+ optional words 4, 5) | real (+ real and/or keyword) | - | `recon_lambda_step0=0` (off: the one-step hand-off) | `recon_lambda_step0`, `recon_lambda_dtu_tol`, `recon_lambda_adaptive` | `Reconstruction continuation: <dlambda> [<dtu_tol>] [fixed\ or adaptive]`. Walks the PLM -> WENO3 hand-off of a two-stage run along the homotopy `R_lambda = (1-lambda) R_PLM + lambda R_WENO3` instead of changing the discrete operator in one step. `dlambda` is the step in lambda per marching step (`<= 0` off; `1.0` reproduces the one-step switch, byte-identical). `dtu_tol` (default `1.2`) is the step control: lambda advances while `dtu <= dtu_tol * dtu` at the start of the ramp, and a step that exceeds it halves the lambda step instead of advancing. `fixed` disables that control; `adaptive` is the default. The two optional fields are told apart by content, not position. Requires `Reconstruction scheme: PLM+WENO3` with a non-empty PLM stage (order dependency on core line 18 and on K13); otherwise a WARNING is printed and the key is ignored. Once lambda reaches 1 the continuation disarms and the run continues in pure WENO3, so the JFNK finish and all acceptance gates see the production operator. Physics and measurements: `docs/f_plm_weno_continuation.md`. |
| K14 | `ATES_photoionization_rate` | word 2 == `True`/`true` | flag | - | `.false.` (Verner 1996) | `ates_photoion_rate` | Reverts He I (1^1S) photoionization to the legacy ATES fit. |
| K14b | `Legacy_HHe_rates` | word 2 == `True`/`true` | flag | - | `.false.` (Badnell/Mao + Voronov) | `legacy_hhe_rates` | Reverts H/He case-B recombination and collisional ionization to the legacy ATES fits (Hui & Gnedin 1997 recombination; Abel+1997/HG97 collisional ionization). Default uses Badnell RR (+ He II DR) minus Mao & Kaastra 2016 alpha_1 for case B, and Voronov 1997 collisional ionization. Free-free always uses the van Hoof et al. 2014 Gaunt table. |
| K14c | `Secondary_ionization` | word 2 == `False`/`True`/`Immediate` | flag | - | `.true.`, STAGED (SvS85 on after first convergence) | `use_sec_ion`, `sec_ion_immediate` | Shull & van Steenberg (1985) secondary ionization by fast photoelectrons (E0 > 40 eV): heating is scaled by f_heat(x) and H I / He I gain secondary ionizations; x is the ionized fraction of the H+He nuclei. STAGED activation (Update §38): the coupling is applied only after the wind first converges without it, then the run re-converges (stops are held N_stall steps after the flip). From a cold IC the immediate coupling amplifies the base startup transient into a NaN runaway on high-gravity cases. `False` = full photoelectron thermalization (bit-identical to the legacy path); `Immediate` = apply from step 0 (pre-staging behavior, A/B tests only). |
| K14d | `He_rec_coupling` | word 2 == `True`/`False` | flag | - | `.true.` (photons ionize/heat H; `False` = legacy lost-photon path) | `use_he_rec_coupling` | Couples He II -> He I recombination radiation to H I ionization (Draine 2011 on-the-spot y/z; see docs/QUESTIONS_2026-07-17.md). Adds an extra H I photoionization rate `n_HeII n_e [z alpha_B + y alpha_1]` with its photoelectron heating, and corrects the He II recombination to `alpha_B + y alpha_1`. `y` (Eq. 14.16) is the local `>= 24.6 eV` ground-capture share ionizing H; `z` (Sec. 15.5) is the density-dependent cascade share. `False` = pure case B (legacy; in TR mode the singlet recombination is then neither case A nor case B). |
| K14e | `He_H_charge_exchange` | word 2 == `True`/`False` | flag | - | `.true.` (He<->H pair active in every He system) | `he_h_charge_exchange` (charge_exchange module) | The He<->H charge-exchange pair (Huang 2023 Table 4 group B = Koskinen 2013 rates): B1 `He0+H+ -> He+ +H0` (endothermic, exp(-128000/T)) and B2 `He+ +H0 -> He0+H+` (the He+ loss channel where neutral H dominates, ~40% of He II at r~1.05 on WASP-121b, <1% by r>=1.2). Applied by dedicated routines in EVERY ionization system with He (including the advection pair and the analytic Newton Jacobians), independent of `cx_full` (which still gates the metal+He / metal+metal groups C/D). `False` = legacy no-He-CX path (bit-identical). |
| K14f | `Atomic rate set` | word 4 == `Koskinen2022`/`default`; any other word stops the run | flag | - | `default` (EXHALE's own atomic H/He rates) | `atomic_rate_set_k22` | Swaps the four atomic H/He rate coefficients for the fits Koskinen et al. (2022, ApJ 929, 52) list in their Table 1 as R1-R4, for a like-for-like comparison with their Model A: R1 `H+ + e -> H + hv` = `4.0e-12 (300/T)^0.64`, R2 `He+ + e -> He + hv` = `4.6e-12 (300/T)^0.64` (both Storey & Hummer 1995), R3 `H + e -> H+ + 2e` = `2.91e-8 U^0.39 exp(-U)/(0.232+U)` with `U = 13.6 eV / kT`, R4 `He + e -> He+ + 2e` = `1.75e-8 U^0.35 exp(-U)/(0.180+U)` with `U = 24.6 eV / kT` (both Voronov 1997, whose general form carries a `(1 + P sqrt(U))` factor with `P = 0` for these two). R3/R4 are the same Voronov fit, with the same parameters, that EXHALE already uses by default, so only the two recombination coefficients actually change value; the branch does pin the Voronov fit against `Legacy_HHe_rates`. The default set is the physically preferred one here -- case B recombination, the Lyman continuum being optically thick in this gas -- so the key exists to reproduce their choice, not to replace ours. It swaps RATE coefficients only: the recombination COOLING coefficients are untouched, with one deliberate exception: `lambda_rec_HeII` IS `kB T` times the He II recombination coefficient itself, so it follows the switch -- decoupling them would remove electrons at one rate and charge the gas at another. `lambda_rec_HII` is an independent Hui & Gnedin (1997) fit and does not move. Both sets are read by every consumer of the four (the coupled ionization systems, the molecular systems, `T_equation`, the steady residual, the cooling that uses alpha), because the switch sits inside the single definition of each coefficient in `Cool_coeff.f90`. With `He_rec_coupling` on the He II ground-capture `alpha_1` stays the Mao & Kaastra fit (Table 1 has no such split), so the closest comparison run also sets `He_rec_coupling: False`. |
| K15 | `Molecular chemistry` | word 3 == `True`/`true` | flag | - | `.false.` | `thereis_mol` | Requires He. Metals are allowed: they are solved in the same system as the molecular network. `He_diffusion` is allowed too since 2026-08-26 (milestone M4 of docs/binary_diffusion_design.md): the element transport closes over the molecular carriers. |
| K15b | `Stellar LW flux` | word 5 | real | erg/cm^2/s | `0.0` | `F_LW_star` | Band-integrated stellar flux over 912-1201 A at the planet's orbit. Drives H2 photodissociation in the molecular network (`lyman_werner.f90`): unattenuated rate `1.757e-7 * F_LW` s^-1 times the Draine & Bertoldi (1996) eq. (37) self-shielding factor of the star-ward H2 column, plus 0.4 eV of heating per dissociation. With K15c on it is **also** the first oxygen photolysis band (`LW`): 912-1201 A is one wavelength interval with one incident flux, and H2, H2O and OH absorb out of one beam, so the H2 rate carries the H2O + OH continuum factor of DB96 eq. (40) (identically 1 without K15c) and the H2O/OH rates there carry the fraction of the band the H2 lines have removed. 0 = off (bit-identical to the network without it); a warning is printed if the key is set without `Molecular chemistry`, and another if K15c is on while this is 0 (no H2O/OH photolysis over the interval where their cross sections peak). **Resolved rather than merely read (`Update_EXHALE_stage1.md` section 134).** With `Molecular chemistry` on, this key absent and a numerical spectrum in use (`Spectrum type: File`), the value is computed from that spectrum by `lyman_werner_band_flux_from_sed`: a trapezoid of the SED over 912-1201 A. No `(R_star/a)^2` dilution is applied because EXHALE's SED file is already AT THE PLANET (`read_sed`); the dilution belongs to the stellar-surface files the value used to be produced from by hand, and the two routes agreed to 4% on the solar proxy when the band was 912-1110 A (328.96 against the 343.0 the molecular cases then carried; they now carry 480.9 over 912-1201 A). A stated key always wins over the integral. With neither a key nor a spectrum the flux stays 0 and `write_setup_report` prints a WARNING naming it a missing channel rather than a modelling choice: the band exists whenever the star does. Whether the number came from the key or from the spectrum is recorded by `lw_from_spectrum` and echoed in the setup report. |
| K15c | `Oxygen chemistry` | word 3 == `True`/`true` | flag | - | `.false.` | `thereis_oxychem` | The A2 option: OH, H2O and CO solved in the coupled molecular ionization equilibrium, with the FUV photolysis of H2O and OH, so the base H2/H partition is computed rather than imported. Requires K15, He/H > 0 and a non-zero oxygen abundance (`metals.inp` `X_O` or `base.inp` `O_H_base`); each is a fatal `error stop`. Refuses `base.inp` `q_H2_base` (it computes that partition) unless a lower-atmosphere profile owns the region below the matching level, and refuses `He_metal_diffusion` (that arm counts an element over its ion stages alone, and CO carries an oxygen *and* a carbon nucleus, so it cannot follow two element factors). `He_diffusion` alone is accepted. Appends `OH H2O CO` to `Ion_species.txt` and writes `output/Oxygen_chemistry.txt` and `output/FUV_bands.txt`. See `docs/a2_oxygen_option_design.md`, `docs/a2_reaction_audit.md`. |
| K15d | `Molecular carrier transport` | word 4 == `True`/`true` | flag | - | `thereis_oxychem` (on whenever K15c is on; off for a molecular-only run) | `carrier_transport` | Vertical transport of the molecular carriers -- H2 always, plus OH, H2O and CO when the oxygen cycle is on: an implicit backward-Euler diffusion-advection step solved together with the same chemistry rows the local equilibrium uses (`diffusive_photochemistry.f90`). Molecular diffusion by Blanc's law over the background carriers, the eddy coefficient from `kzz_cell` (K19 or a profile, and nothing else), settling, no thermal diffusion. `False` restores the local steady state of milestone M2, a test of the chemistry alone, not a model of a base, since at a cool base the H2 chemical time and the flow time are comparable. With `kzz_cell` zero everywhere the transport is pure molecular diffusion and `write_setup_report` warns. Recorded in `EXHALE_resolved.out` as `carrier_transport`, beside `oxygen_chemistry`, `oxygen_reaction_set`, the five `fuv_band_*_flux` values and `oxygen_base_partition`. Acts only with K15 (`Molecular chemistry`), not K15c: every consumer of `carrier_transport` is guarded by `thereis_mol`, so with the molecular network off there are no carriers and the key transports nothing. That is not refused but reported: `input_read` warns that the key has no effect, as it does for `Molecular IR bands` and `Stellar LW flux` without the network. The retired name `Oxygen transport` is refused at startup rather than aliased (`Update_EXHALE_stage1.md` section 128). |
| K15e | `Coupled carrier solve` | word 4 == `True`/`true` | flag | - | `False` | `carrier_in_newton` | Solve `n(H2)` as a FOURTH Newton unknown per cell -- `(rho, rho v, E, n(H2)/n0)` -- instead of alternating a wind solve with a fixed-wind carrier relaxation. The alternation was measured not to converge and not to be able to: its drift gate stands below the smallest move the transport operator can make, and because the two halves are coupled through the particle count (which sets the temperature) a splitting that evaluates each at the other's old state gets the SIGN of the H2 front's motion wrong. With this on, the outer Picard loop runs once, the carrier row joins ` or norm(R) or ` on the scale of its own largest terms, the element budget becomes an admissibility constraint rather than a clamp, and the band geometry widens to `kl = ku = 11` with 23 Jacobian colors. On a MOLECULAR configuration (K15) it needs K15d (`Molecular carrier transport`) and is refused at startup without it (item B5i, 2026-09-08: with H2 eliminated the stationary residual is not a function of its unknowns, MEASURED seed dependence 3.2e5 times `Resid tol` per unit relative seed change against 9.7e-10 with the H2 row carried); in an atomic gas with `He_diffusion: True` it registers an element row and needs no carrier transport. The carrier densities of that system are carried as `ln n` (`EXHALE_CARRIER_LOG_UNKNOWN=0` restores the density unknown for measurement): the positivity of a carrier is then a property of the unknown space rather than a bound the step control has to enforce, with a floor at 1e-20 of the cell's own element budget, the fraction the carrier row scale already calls unobservable (item B5k, 2026-09-08; `docs/steady_solver_design.md` section 20). Default off: the coupled route changes the size of the steady system, so a run that does not ask for it does not pay for it (`Update_EXHALE_stage1.md` section 139). |
| K15f | `Ionization transport` | word 3 == `True`/`true` | flag | - | `False` | `ionization_transport` | Carry the HYDROGEN IONIZATION STATE with the flow: H+ becomes a fifth transported carrier of the same operator that carries H2 (`ic_Hp`, `diffusive_photochemistry.f90`), and the ionization sweep is handed the transported fraction instead of solving the H/H+ partition as a local root (`x_hp_fix` replaces row 1 of the molecular system, the way `x_h2_fix` replaces row 4). WHY: the local partition is the right answer only where a parcel is ionized faster than it leaves its shell, `P r/ or v or  >> 1`. On the Koskinen 2022 Model A comparison that ratio is 0.15-0.35 above 1.5 r_base, and the local root then gives x(H+) = 0.129/0.327/0.585 at 2.0/2.4/3.0 r_base where integrating along the same flow gives 0.057/0.078/0.106 and Model A, which advects every species, has 0.020/0.055/0.080 (`docs/k22_electron_density_excess.md` sec. 7). The proton's MOLECULAR diffusion is set to zero -- an ion in a neutral gas is ambipolar/Coulomb-limited, not hard-sphere, and the term this option exists for is the bulk advection -- while the eddy coefficient `kzz_cell` still applies, it being a bulk mixing coefficient. No base Dirichlet value: no handoff states an ionization fraction, so the lower ghosts stay with the sweep. Recorded in `EXHALE_resolved.out` as `ionization_transport`, and in a restart file's `# coupling:` line as `iontrans=T`. Needs K15d (`Molecular carrier transport`) and K15 (`Molecular chemistry`), and with `Solver: Newton` it needs K15e (`Coupled carrier solve`) as well: the proton is a transported carrier, so with K15e the stationary system carries its row and the sweep is handed the Newton unknown (`x_hp_fixed`), while without it the proton is outside the unknown space and the last equilibrium sweep would put back the local ionization state this option exists to leave. Each of the three is a fatal `error stop` in `input_read`. Marching needs none of them. The direct steady route (`EXHALE_PTC=1`) is refused outright, for the same reason as the sweep. Default off: it changes every ionization-dependent number of a run. |
| K15g | `Stellar FUV B1 flux` | word 6 | real | erg/cm^2/s | - | - | **RETIRED 2026-09-06 (item LW-NORM-B).** Band `B1` was 1110-1201 A, between the Lyman-Werner interval and Ly-alpha. The H2 Lyman and Werner lines pump on both sides of 1110 A at the 700-3200 K of a planetary base, so a band edge there normalized the H2 pumping per photon of a band narrower than the one the lines absorb from. B1 is now part of the Lyman-Werner band, which is 912-1201 A. A file that still carries this key STOPS the run with a message saying to delete the line and set K15b to the sum of the two numbers. |
| K15h | `Stellar FUV B3 flux` | word 6 | real | erg/cm^2/s | `0.0` | `F_FUV_B3` | The same over 1231-1450 A (band `B3`, branching 0.89/0.11/0.00). Band `B2` is the Ly-alpha line and is supplied by K10. |
| K15i | `Stellar FUV B4 flux` | word 6 | real | erg/cm^2/s | `0.0` | `F_FUV_B4` | The same over 1451-2304 A (band `B4`, branching 1.00/0.00/0.00). The H2O cross section falls three decades across this band, so the flat-`F_lambda` band average is off by a factor 4.6-6.1 against a real stellar spectrum; a run in which B4 matters is outside the treatment (`water_photolysis.f90` section 2). |
| K16 | `Molecular base` | word 3 == `True`/`true` | flag | - | `.false.` | `molecular_base` | EOS-only molecular base correction to `ntot_bc`. |
| K17 | `Lower atmosphere` | word 3 (+ optional word 4) | string + real | - / R_J | `lower_atm_mode=0` | `lower_atm_mode`, `lower_atm_r1bar` | `none`/`analytic`/`vulcan`. Triggers `run_lower_atm_prestep` (needs `EXHALE_ROOT`). |
| K17b | `Lower atmosphere profile` | word 4 | string | file name | `''` (off) | `lap_file` (`lower_atmosphere_profile.f90`) | The lower atmosphere handed over as a **table over an interval of pressure** instead of the single-level scalars of `base.inp` (section 2d). Matched BEFORE K17, whose label is a prefix of this one. With the key set the named file must exist. The profile then owns the base state, the elemental reservoirs and `K_zz`, and the `base.inp` keys of those three categories are **refused** (section 2c). Its `p_match_bar` is the level of the run: `n0 = p_match/(k_B T0 ntot_bc)` exactly as `p_base` of `base.inp` sets it, and a `Log10 lower boundary number density` beside the profile is refused unless it agrees within 1%. Absent key = present behavior, bit for bit. |
| K18 | `Lower column` | word 3 | real | R_J | `lower_col_r1bar=-1` | `lower_col_r1bar` | Analytic lower column 1-bar radius. |
| K19 | `He_Kzz` | word 2 | real | cm^2/s | `0.0` | `he_kzz` | The constant eddy diffusion coefficient a run states when it has **no** profile; it fills every entry of `kzz_cell`, which is what the element-diffusion operator (K23) and the molecular carrier transport (K15d) both read. Default 0 = pure molecular diffusion (it was `1.0e9` before 2026-08-25). Also overridable by `base.inp` (`Kzz_base`). With a lower-atmosphere profile in use the key is inert -- `K_zz` is then a profile -- and a warning says so. |
| K20 | `He_alphaT` | word 2 | real | - | `0.0` | `he_alphaT` | Thermal-diffusion factor. |
| K21 | `He_ambipolar` | word 2 == `False`/`false` | flag | - | `.true.` | `he_ambipolar` | Only `False` changes it (default on). |
| K22 | `He_metal_diffusion` | word 2 == `True`/`true` | flag | - | `.false.` | `he_metal_diffusion` | Tested before K23. |
| K23 | `He_diffusion` | word 2 == `True`/`true` | flag | - | `.false.` | `he_diffusion` | He/H diffusive separation. |
| K24 | `Stall` | words 3, 4 | real, int | - , steps | `stall_tol=1e-6`, `N_stall=2000` | `stall_tol`, `N_stall` | Stall-detector override. |
| K25 | `Energy solver` | word 3 == `Explicit` | flag | - | semi-implicit (`.true.`) | `use_semi_implicit_energy=.false.` | `Explicit` reverts the source update to forward Euler; any other word leaves the semi-implicit update in place. |
| K26 | `Time stepping` | word 3 == `Local` | flag | - | global (`.false.`) | `use_local_dt=.true.` | Cell-by-cell pseudo-time. |
| K27 | `Level tol` | word 3 | real | - | `lev_th=-1` | `lev_th` | Mass-flux level-stability tolerance. |
| K28 | `Solver` | word 2 == `Newton` (+ optional word 3) | flag + real | - | `use_newton_solver=.false.`, `newton_du_switch=1e-2` | `use_newton_solver`, `newton_du_switch` | JFNK hand-off. Distinct from K41/K42 (see §6.9). |
| K29 | `Valve eps` | word 3 | - | - | - | - | **Retired 2026-09-03 (section 152).** The label is still matched, and a file that carries it is REFUSED at startup with a message naming the replacement. There is no ghost velocity closure to select: the lower boundary takes the velocity from the outgoing acoustic characteristic of the first interior cell. |
| K30 | `Hydrostatic base` | word 3 | - | - | - | - | **Retired 2026-09-03, refused as K29.** There is no ghost pressure closure to select: the face pressure is the reservoir's, carried to the face along its own hydrostatic isentrope, and the ghost cells are the volume averages of that isentrope continued below the face. |
| K31 | `Shapiro filter` | word 3 (+ optional word 4) | real, int | - , steps | `shapiro_eps=-1`, `shapiro_every=4` | `shapiro_eps`, `shapiro_every` | Word 3 is the filter amplitude (`<= 0` disables it), optional word 4 the period in steps. Off by default; opt-in for a breathing base. |
| K31b | `Low-Mach damping` | word 3 (+ optional word 4) | real, real | - , Mach | `lowmach_damp_eps=-1` (off), `lowmach_damp_mach_th=1e-3` | `lowmach_damp_eps`, `lowmach_damp_mach_th` | Gated fourth-difference (Jameson-Schmidt-Turkel) dissipation of the `2 dr` contact/entropy mode that the contact-resolving HLLC flux stops damping as `v -> 0`. Added to the numerical momentum flux (and its work term to the energy flux) inside `RK_rhs`, so the marching loop and the JFNK steady residual see the same equation -- unlike K31, which touches the marching state only. The gate `[max(0, 1 - M^2/M_th^2)]^2` is exactly zero for `M >= M_th`. Explicit stability needs `eps4 < 1/(16 CFL)`; `input_read` warns otherwise. `eps4 <= 0` = off (default, bit-identical). See `src/modules/flux/low_mach_dissipation.f90`, `docs/hd209_metal_stagnation.md`. |
| K32 | `Base BC` | word 3 (+ optional word 4 if `pressure`) | string + real | - / microbar | `base_bc_mode=0` (density), `base_p_ubar=1.0` | `base_bc_mode`, `base_p_ubar` | `density` or `pressure`. It sets the base LEVEL, not the boundary closure, and it is the only base key left. The pressure it names is the pressure of the lower-atmosphere reservoir AT `r = 1` (the planet radius, the level `base.inp`'s `p_base` also refers to), not at the first cell face: `base_boundary` carries the reservoir from that level to `r_edg(0)` along its own hydrostatic isentrope, so refining the grid does not move the level the user stated. Pressure mode derives `n0`. |
| K32b | `Base ghost temperature` | word 4 | - | - | - | - | **Retired 2026-09-03, refused as K29.** There is no ghost temperature closure to select: the face temperature follows from the reservoir pressure and entropy at the base composition. |
| K32c | `Max steps` | word 3 | int | steps | `count_max=1000000` | `count_max` | Hard cap on marching iterations. The env variable `EXHALE_MAXSTEPS` is separate and only exits earlier. |
| K32d | `Coronal cutoff width` | word 4 | real | - | `coronal_cutoff_width=0.1` | `coronal_cutoff_width` | Roll-off width of the coronal-excitation guard below the 1e3 K CHIANTI fit floor (`Cool_coeff.f90`). Must be > 0; input_read aborts otherwise. Since the ground-term fine-structure statistical equilibrium (2026-08-12) the base temperature is insensitive to `w` over 0.02-1.2, so the value is no longer a tuning knob; the earlier 0.08-0.13 justification window is superseded (`docs/coronal_cutoff_width.md` section 7.2). |
| K32e | `Base IR field` | word 4 == `True`/`true` | flag | - | `.false.` | `base_ir_field` | Lets the infrared coolants of a molecular layer see the thermal radiation of the atmosphere below the base instead of emitting into vacuum: the lower atmosphere is taken to be black at those wavelengths and to radiate `B_nu(T0)` over the sky fraction `1 - sqrt(1 - (R_p/r)^2)`. Applies to the eight ground-term fine-structure lines of C I, C II, N II, O I (the incident field enters the statistical equilibrium as a photon occupation number, and the escape probability becomes two-sided) and to the H3+ bands (whose absorption is the Miller et al. 2013 total emission fit evaluated at `T0`, by Kirchhoff's law, since that fit is not a line list). Returns the net rate, emission minus absorption, so each channel stops cooling at its own radiative-equilibrium temperature. Every other channel keeps the optically thin, no-incident-field limit. Off = bit-identical to the previous single-face form; no effect on an atomic run. See `fine_structure_line_transfer` (`Cool_coeff.f90`), `h3p_net_cooling_rate` (`h3p_cooling.f90`), `docs/lower_atmosphere_coupling.md`. |
| K32f | `Molecular IR bands` | word 4 == `True`/`true` | flag | - | `.false.` | `mol_ir_bands` | Adds the infrared coolants a real H2 atmosphere carries below the H2 -> H front: the H2 quadrupole and magnetic dipole line spectrum (Roueff et al. 2019, A&A 630, A58), and the H2O and CO vibration-rotation bands (HITEMP through the Photochem correlated-k coefficients). Each emits in LTE and absorbs the same diluted `B_nu(T0)` K32e supplies, so each stops cooling at its own radiative equilibrium temperature; the closure keeps the stimulated-emission term, so the net rate is exactly zero for gas at `T0` with the whole sky black. Three new NET columns in `output/Cooling_breakdown.txt` (`H2_IR`, `H2O_IR`, `CO_IR`, negative where the band heats) and a column-integrated infrared block in `output/FUV_bands.txt`. H2 needs `Molecular chemistry: True`; H2O and CO need `Oxygen chemistry: True`, and the key is inert on species the run does not carry. With K32e off the bands emit into vacuum, which deepens the collapse rather than holding the layer, and `input_read` warns. Optically thin, LTE, cross sections tabulated 50-2000 K; the Planck-mean optical depths are written to `output/Cooling_breakdown.txt` so the thin assumption is measured. See `src/modules/lower_atmosphere/molecular_infrared_cooling.f90`, `cooling_data/molecular_infrared_bands.py`, `docs/lower_atmosphere_coupling.md`. |
| K32g | `Molecular reaction heat` | word 4 == `True`/`true` | flag | - | `.true.` | `mol_reaction_heat` | **Default ON; the key exists to turn the term OFF.** Deposits the energy the COLLISIONAL reactions of the H2/He network release into the gas. Why it cannot be left out: in an atomic gas the ionization energy a photon spends comes back as the Lyman photon of a radiative recombination and leaves, so it is correctly absent from the heating; in a molecular gas it does not. The cycle `H2 + hv -> H2+ + e` (the photon pays `I(H2) = 15.43 eV` and the photoelectron keeps the rest, which `PH_heat_HHe` already deposits), `H2+ + H2 -> H3+ + H` (+1.70 eV), `H3+ + e -> H2 + H` (+9.25 eV, DISSOCIATIVE) returns to H2 having turned one H2 into H + H and left `I(H2) - D0(H2) = 10.95 eV` in the gas with no photon to carry it away. Measured at 81 percent of the total heating rate at 1.02 `r_base` on the converged He/H = 0.0793 rung. The heat is built from ONE table of species enthalpies, so a closed chemical cycle releases exactly zero; the photon-driven reactions, the radiative recombinations, the collisional ionizations, the Penning channels and the H/He charge exchange are excluded because their energy is already in the ledger. The three-body association carries the Hollenbach & McKee (1979) collisional branching `(1 + n_cr/n)^-1`. Gated on K15 (`Molecular chemistry`), so an atomic run is byte-identical either way; `False` reproduces the state the code was in before it. New column `heat_mol_chem` in `output/Heating_breakdown.txt`, and one line in `EXHALE_setup.out`. See `src/modules/lower_atmosphere/molecular_reaction_heat.f90`, `Update_EXHALE_stage1.md` section 141. |
| K33 | `Base velocity` | word 3 | - | - | - | - | **Retired 2026-09-03, refused as K29.** Prescribing rho, v and p together is one condition too many for a subsonic inflow face. |
| K33b | `Base grid` | words 4 and 5 | real, int | R_planet, cells | `dr_base=2.0e-4`, `N_low_cells=50` | `dr_base`, `N_low_cells` | `Base grid [dr,cells]: <dr_base> [<N_low_cells>]` sets the uniform region of the `Mixed` grid: `N_low_cells` cells of size `dr_base` stacked on the base. The two are not independent: their product is the radial extent of that region (0.01 R_p by default), so they share one line, as `du_th [PLM,WENO3]` does; refining at fixed extent means dividing the first and multiplying the second (`5.0e-5 200` is the 4x refinement). `dr_base` must resolve the base scale height `H = kT/(mu g)`, which `write_setup_report` echoes as cells per `H`. Ignored by the `Uniform` and `Stretched` grid types. The default is written as a default-real literal, so spelling it out in `input.inp` does **not** reproduce a no-key run bit-for-bit. See `docs/hd189_base_checkerboard.md` §11. |
| K33c | `Grid cells` | word 3 | int | cells | `N=500` | `N` | Number of computational cells of the radial domain; ghost cells are added on top and are not counted. Fewer than 10 is a fatal `error stop 1`. Omitting the key keeps the 500 that used to be a compile-time constant, so an existing `input.inp` is unaffected. Every grid-sized array is allocated by `allocate_grid_arrays` once `N` is known. For the `Mixed` grid the split between the uniform base region and the stretched region is set separately by K33b, and `define_grid` checks the two are compatible. |
| K34 | `Viscosity` | word 2 (+ optional word 3) | `True`/`False` or real, real | - | `visc_on=.false.`, `visc_mu0=0` (off), `visc_s=0.7` | `visc_on`, `visc_mu0`, `visc_s` | `True` = calibrated `mu(T)` + dissipation `q_mu`; a number = diagnostic power law `mu0*T^s` in code units. One-word key, so the value is word 2 (was word 3, which no input file used). See `docs/viscosity_conduction.md`. |
| K34b | `Conduction` | word 2 | `True`/`False` | - | `cond_on=.false.` | `cond_on` | Heat conduction with `kappa(T) = 4.45e4 (T/1000 K)^0.7` (Watson+1981). Independent of K34. |
| K34c | `Caloric EOS` | word 3 == `ladder`/`monatomic` | string | - | `ladder` | `caloric_eos_monatomic` | The caloric equation of state. `ladder` (default): every atom, ion and electron stores 3/2 kT and H2 its rovibrational energy (Roueff et al. 2019 ladder), so gamma_eff depends on composition and T (section P53). `monatomic`: H2 too stores 3/2 kT, gamma = 5/3 everywhere -- the state of the code before P53, kept as a COMPARISON option for a published model whose energy equation is u = c_v T with a monatomic c_v (the Koskinen et al. 2022 Model A comparison, 2026-09-05). Any other word is refused. |
| K34d | `Photoelectron heating` | word 3 == `excess`/`full`/a fraction in (0,1] | string or real | - | `excess` | `photoheat_photon_fraction` (< 0 = excess; `photoheat_full_photon_energy` = fraction 1) | What a photoionization deposits as heat: `excess` (default, physical): the photoelectron's h nu - I. `full` (= 1): the whole photon energy h nu; a number f: the fraction f of h nu. Both are COMPARISON options (2026-09-05/06): Koskinen et al. (2022) Figure 9 shows a stellar heating rate 2.5-3 times what h nu - I gives for their ionization rate and spectrum, between the two accountings; a fraction reproduces their profile without claiming a physics, and a run with it must say so. Anything else is refused. |
| K35 | `Resid tol` | word 3 | real | - | `resid_th=-1` (the solver then uses 1e-5) | `resid_th` | Residual-norm convergence instead of du. **Renormalized 2026-09-03 (section 133)**: the rows are divided by a bound on their own largest term, not by ` or u or `, so every value is 1e-2 times its pre-133 equivalent (the old default 1e-3 is now 1e-5). **Section 143 changes what that bound is, on all three rows**: each row is divided by the largest term the row itself contains -- mass `max( or F or  r^2 at the two faces)/dV` (`mass_flux_row_scale`), momentum `max( or dF_2 or , or S_2 or )` (`momentum_row_scale`), energy `max( or dF_3 or , or S_3 or ,heat,cool)` (`energy_row_scale`), plus the operator-split viscous and conduction sources where those are on. The previous scale was the cell's state scale times its signal-crossing rate, `D_k( or v or +c_s)/dr`, and nothing in any of these rows moves at `c_s`: on the mass row it made the residual read the flux error times the local Mach number, a factor 2e4 too small in a quasi-hydrostatic layer, and the momentum and energy rows were mis-scaled the same way by 1e3 to 1e5. The section 62.3 gravitational bound is subsumed -- `S_2` is that term, taken from the expression the residual subtracts. **No value in the repository is renormalized by section 143** (the states the new measure accepts pass at the existing tolerance), but a run that reaches its answer through a steady solve or this gate stops at a different state, so its golden moves; a run that stops on `du` alone is byte-identical. **Section 145 changes how the scaled rows are combined into the one number this key is compared against**: it is now the maximum over cells of that cell's own scaled residual, `max_j  or R_kj or  / s_kj`, so a single cell out of tolerance keeps the state out of tolerance. The two forms the code carried before -- a ratio of volume-weighted sums and a ratio of maxima -- both averaged the cell scale away and neither was the local relative norm; measured on the two accepted roots of section 144, the volume form hid 7 cells on the hot Uranus and 70 on WASP-121b, 50 of them in the wind. Both are still REPORTED beside the gate value (section 145.4); only the acceptance changed. The tolerance itself is unchanged. |
| K35b | `Flux spread tol` | word 4 (+ optional word 5) | real, real | -, R_p | `flux_spread_th=2.0e-5`, `r_flux=1.2` | `flux_spread_th`, `r_flux` | The FLUX gate (section 133; redefined in section 145): a steady solve is accepted only when the radial spread of the RIEMANN FACE mass flux `F_{j+1/2} r^2` over `r >= r_flux` is below the tolerance AND the residual is below `Resid tol`. `<= 0` disables the flux gate. The defaults are unchanged by section 143, and the window still starts at 1.2 R_p for the reason `flux_spread_of_state` gives -- the base deliberately admits a small inflow, so the launch region measures the boundary condition rather than the wind. What section 143 changes is that the region the window excludes is now held by the residual gate instead of by nothing: on the states it produces the mass flux is flat to 1e-4 from 1.01 R_p outward. The spreads over `r >= 1.03` and `r >= 1.10` are reported beside the gate value (section 142). The marching `du` stop is **not** this gate -- different window, different functional, different threshold -- and section 143.5 records a matrix golden that stops on `du` at six times this tolerance. **Section 145 changes both the functional and the default.** The gate used to measure the CELL-CENTRED product `rho v r^2`, which the scheme does not conserve: on an accepted root whose face flux is a single double-precision value over 300 cells, that product spreads by `3.4e-3`, five orders larger, because it reconstructs the flux from cell averages. The gate now reads the face flux the HLLC solver actually returns, `face_mass_flux_r2`, and the threshold is re-derived for it: over the converged runs in the repository the face spread of an accepted root is `4.6e-13` to `1.6e-11`, and `2.0e-5` sits far above that and far below anything not converged. **An `input.inp` written before section 145 that still carries the old `5.0e-3` gets a warning from `input_read`, not a refusal**: on the new functional that value is six orders loose. |
| K36 | `Resid norm` | word 3 | string | - | RETIRED (section 145) | - | Was `vol`/`volume` or `Linf`/`linf`/`LINF`, selecting between two ways of combining the scaled rows. Section 145 makes the acceptance norm the cellwise maximum and reports the integrated norm beside it every time, so there is nothing left to select. The key is still READ and produces a warning saying so; `resid_vol` is gone from `parameters.f90`. |
| K37 | `CFL` | word 2 | real | - | `CFL=0.6` | `CFL` | Overrides the CFL number; a lower value gives a smaller `dt`. |
| K38 | `Transonic IC` | word 3 == `True` | flag | - | `.false.` | `transonic_ic` | Transonic isothermal-wind initial condition (`set_IC.f90`). Also sets `force_start`, so the first iterations run before the convergence test resumes. |
| K39 | `Hot Parker IC` | word 4 | real | K | `hot_parker_ic=.false.`, `T_wind_ic=1e4` | `T_wind_ic`, `hot_parker_ic` | Reads word 4 unconditionally when present. |
| K40 | `IC mode` | word 3 | string | - | `ic_mode=0` (cold) | `ic_mode` (+ `transonic_ic`/`hot_parker_ic`) | `cold`/`transonic`/`hot_parker`/`auto`/`windae`. Synonyms for K38/K39; explicit legacy keys take precedence. Unknown value warns and falls back to cold. |
| K41 | `Newton solver` | word 3 == `False` | flag | - | `use_newton_ieq=.true.` | `use_newton_ieq` | Ionization-equilibrium Newton solver toggle. |
| K42 | `Brent solver` | word 3 == `False` | flag | - | `use_brent_tsolve=.true.` | `use_brent_tsolve` | Temperature Brent solver toggle. |
| K43 | `Restart intent` | word 3, and word 4 with `stationary` | string | - | `trajectory` with `Run mode: phys`, else `relaxation` | `restart_intent`, `stationary_evaluate_only`, `stationary_equilibrate_loaded` | What a loaded state is continued as (`docs/restart_contract_design_20260909.md` section 2). `trajectory`: the clock continues from the state file's `t_phys`. `relaxation`: march toward stationarity, then the solver hand-off (what every restart did before the key existed). `stationary`: rebuild the derived quantities without a step, measure the stationary residual and the certification of the state AS LOADED, then enter the stationary solve at once; the second word `evaluate` writes the measured state back and stops instead, and `equilibrate` puts the loaded composition on its own fixed point before the measurement. Three refusals: the key without `Load IC? True`, `trajectory` with `Run mode: init`, `stationary` without `Solver: Newton`. |
| K44 | `Restart option change` | the rest of the line, tokens separated by commas, blanks or tabs | string list | - | no option may differ | `restart_option_change_named` (`IC_load`) | Which tokens of the state file's `# options` line are ALLOWED to differ between the loaded state and this run (decision 21 of `docs/To_be_determined_by_user_20260906.md`, option a). The restart contract otherwise refuses any option difference, which forbids the arm ladder this project converges with: converge without an option, restart with it on, converge again. Naming a token permits exactly that token to differ; every other difference still refuses the load and states the token. A named token that does not differ is reported and nothing else. The change is written into the state this run produces as one `# option_change` line (appendix D.2) and inherited by the rungs that follow. Refusals: an unknown token; a token that decides how many unknowns the state has (`metals`, `mol`, `oxychem`, `carrier`, `carrier_newton`, `iontrans`) -- the rows of the file are then not the rows of this run, which is a cold start and not a restart; a word naming the grid, the reservoir or the constant set instead of an option; the key without `Load IC? True`. |
| K45 | `Run mode` | word 3 | string | - | `init` | `run_mode`, `run_mode_given` | What this run is doing, stated rather than inferred (`docs/a0_run_mode_contract_20260906.md` section 6). `init`: initialization or continuation, which claims no elapsed time and permits local pseudo-time, pseudo-transient continuation and a chemistry that is not at its root, as numerical devices; the physical clock `t_phys` stays where it started and no output reports it. `phys`: physical integration, one global `dt` per step, the clock advanced only by a step accepted in full, the temporal error estimate of the step sampled every 20 accepted steps (`attempted_step.f90`), and the certification refusing a state whose chemistry has no root. `init` is the default in every configuration, and a value other than `init` or `phys` is a fatal `error stop 1`. Refusals: `phys` with `Time stepping: Local`, since cells advanced by different intervals do not form one trajectory. The mode also sets the default of K43 (`trajectory` under `phys`, `relaxation` otherwise), is written into the state file's metadata block as `mode` and read back on a restart (a `phys` run loading a file that states `mode=phys` without a finite non-negative `t_phys` is refused), and is stated in `EXHALE_setup.out` together with whether it was given or defaulted. |

### 2c. Additional startup file: `base.inp`

Not part of `input.inp`, but read at the same point (`read_base_inp` in
`input_read.f90`) and able to override core values after
`input.inp` is parsed. Optional; a missing file is a no-op. Format: keyword
lines, `#` comments ignored. Keys are matched as **labels**, by the same
`lbl_match` the `input.inp` reader uses (anchored at the start of the
left-trimmed line, terminated by `:`, `?`, `=`, whitespace or end-of-line);
before P2 they were matched as bare substrings anywhere in the line.

Every key carries a **category**, and the category is the contract: it says
what the value is allowed to do to the wind. The five categories are the ones
`docs/oxygen_chemistry_new_plan.md` P2 requires.

| Key | Category | Value word | Sets | Notes |
|---|---|---|---|---|
| (comments only) | provenance | - | - | Which code, network and profile produced the file. No key is parsed today; the machine-readable provenance keys are A1a/P0 of `oxygen_chemistry_new_plan.md`. |
| `T_base` | EOS boundary | word 2 | `T0` [K] | Gas temperature at the handoff level; overrides `Equilibrium temperature` (core line 5). |
| `r_base` | EOS boundary | word 2 | `R0` [R_J] | Radius of the handoff level; overrides `Planet radius` (core line 3), so the whole domain is measured from this level. |
| `p_base` | EOS boundary | word 2 | `p_base_bar` [bar] | The level every other value in the file refers to, and **since 2026-09-03 the level of the run**: `n0 = p_base/(k_B T0 ntot_bc)`, so `Log10 lower boundary number density` is no longer needed beside it. If it is given anyway the two must agree to 1%, and `Base BC: pressure` beside it must name the same level to round-off; either disagreement is refused at startup with both numbers and the two ways to fix it. A `Lower atmosphere profile:` file states the level the same way through its `p_match_bar`, and the same 1% rule applies to a density key beside it. Default `1e-6`. |
| `q_H2_base` | EOS boundary | word 2 | `q_h2_base` | H2 volume mixing ratio at the base. With `Molecular base: True` it replaces the chemical-equilibrium fit in `comp_ntot_bc` (`composition.f90`), i.e. it sets the base **particle count**. It is the **composition of the inflowing gas**, not only an EOS anchor (section 117 of `Update_EXHALE_stage1`): it is also imposed on the H2 partition of the lower ghost species, and the base particle count is then taken from that same species state, so the equation of state and the chemistry cannot describe different gas. Refused at startup above the attainable ceiling `0.5/(0.5 + He/H)` rather than capped. Absent (default `-1`) = fit used and the species partition left free, i.e. the historical behavior. |
| `HeH_base` | elemental reservoir | word 2 | `HeH` (sets `thereis_He` if > 0) | He/H nuclei ratio of the inflowing gas; overrides `He/H number ratio` (core line 8). |
| `<El>_H_base` | elemental reservoir | word 2 | `X_<El>`, hence `melem_ab(iel_<El>)` | `El` is any of the ten element symbols of `species_table` (`C N O Mg Si Ca Na K S Fe`), e.g. `O_H_base 4.90e-4`. Nuclei ratio El/H at the handoff level. **Overrides `metals.inp`** for that element, activates the metal system when it is the only nonzero abundance, and, on a restart, renormalizes that element's loaded column onto this ratio (key 20). Since P2, `thereis_metals`, `melem_ab` and `thereis_lowIP_metal` are all derived after the handoff, so a handoff element behaves exactly like a `metals.inp` element. |
| (none) | initial guess | - | - | No key today. A key here would seed a profile the solver may move away from -- e.g. the A0 improvement of using `q_H2_base` as the base-cell H2 seed in `set_IC`. |
| `Kzz_base` | boundary constraint | word 2 | `he_kzz` [cm^2/s] | Eddy diffusion coefficient at the base. It fills `kzz_cell`, which the element-diffusion operator (K23) and the molecular carrier transport (K15d) both read, so it is inert only when neither of those is on. |
| `q_H2O`, `q_CO`, ... | diagnostic | - | - | Written as `#` comments by `vulcan_to_base.py`. No consumer in the code (A1c); metadata, not physics. Unknown keys are ignored, so promoting one to a bare line changes nothing. |

**With a lower-atmosphere profile in use (key K17b) the EOS-boundary,
elemental-reservoir and boundary-constraint keys of this file are refused**
with an `error stop` naming the key and its category. The profile states all
of them at the matching level, and accepting a scalar beside it would
reintroduce the same-solution problem the profile removes by construction.
Provenance comments and diagnostic keys are still allowed, and the one
provenance comment that is parsed is `# solution_id <hash>`: if a `base.inp`
sits beside a profile and either lacks a `solution_id` or the two differ, the
run stops with both ids printed.

Generators: `src/utils/run_lower.py` writes `T_base`, `r_base`, `HeH_base`,
`Kzz_base` (analytic chemical-equilibrium column);
`src/utils/vulcan_to_base.py` writes those four plus `q_H2_base` and `p_base`,
and the molecular mixing ratios as comments. Neither writes `<El>_H_base` yet:
VULCAN's networks are H/C/N/O(/S), so the element abundances it could hand
over are its own input elemental ratios, and metal/alkali release is outside
its scope.

**Element-budget check.** `src/utils/element_budget.py <run_dir>` verifies
that the reservoirs actually held: for every element it compares
`n_El/n_H` per cell against the resolved abundance, hydrogen against
`rho/mass_per_H`. The resolved abundances, `mass_per_H` and `ntot_bc` are
written to `EXHALE_resolved.out` for that purpose.

### 2d. Additional startup file: `lower_atmosphere_profile.dat`

Named by the `Lower atmosphere profile:` key (K17b); the name above is only
the conventional one. It is the lower atmosphere's solution over an interval
of pressure, and it replaces the single-level scalars of `base.inp`. Read by
`read_lower_atmosphere_profile` (`src/modules/files_IO/lower_atmosphere_profile.f90`)
before `read_base_inp`, and applied by `apply_lower_atmosphere_profile` after
it, so nothing scalar can overwrite a profile value and the elemental
reservoirs (`melem_ab`, `thereis_metals`, `thereis_lowIP_metal`) are still
derived after the handoff.

Format: a block of `#` header lines, then a fixed-column table running **deep
to shallow** (strictly decreasing pressure), full double precision. Header
lines are `# key value`; the `# columns:` line names the columns and is what
the reader and `examples/exhale_io.py` index by.

| Header key | Required | Meaning |
|---|---|---|
| `solution_id` | yes | sha256 over the lower model's configuration, mechanism, thermodynamic data, stellar flux and elemental abundances. The fingerprint that makes "same solution" checkable. |
| `source_code` | no | `photochem` / `vulcan` / `analytic` |
| `source_version` | no | which build of the producer wrote the file, e.g. `photochem 0.9.0`. Not decoration: `photochem 0.9.0` is the corrected build this repository carries in `photochem/` and runs the handoff on (`README_photochem.md`), while `photochem 0.8.4` is the conda package, whose equilibrium solver leaves up to `2e-4` in a trace elemental ratio. The two are told apart per file by this line and by nothing else. |
| `mechanism` | no | mechanism file name and its own sha256 |
| `stellar_flux` | no | flux file name, sha256 and the dilution applied |
| `p_match_bar` | yes | the matching pressure: where EXHALE places its base |
| `p_top_bar` | yes | the shallowest level carried; must be `< p_match_bar`, else there is no overlap and the file is refused |
| `p_deep_bar` | no | the deepest level carried (diagnostic) |
| `trial_flux_H`, `trial_flux_He` | no | the elemental fluxes imposed at the lower model's upper boundary for this solution, g/s outward positive |
| `iteration` | no | closure iteration index. Its absence does not stop the run; it makes the profile non-iterable (`lower_profile_iterable F`), which is the right answer for a hand-written one-shot file. |
| `reached_steady_state` | no | `T`/`F` from the chemistry solver |
| `notes` | no | free text |

| Column | Required | Unit | Consumed as |
|---|---|---|---|
| `p` | yes | bar | the interpolation abscissa; strictly decreasing |
| `r` | yes | R_J | `R0` at the match; the radius axis of the `K_zz` interpolation |
| `T` | yes | K | `T0` at the match |
| `n_tot` | yes | cm^-3 | carried, not imposed: the base density is `n0` of `input.inp` and the base particle count follows from `T0`, `HeH`, `q_H2` through `comp_ntot_bc` exactly as on the scalar path |
| `rho` | yes | g cm^-3 | as `n_tot` |
| `Kzz` | yes | cm^2 s^-1 | `kzz_cell`, interpolated onto the grid |
| `q_H2` | yes | - | `q_h2_base` at the match |
| `q_H` | yes | - | carried |
| `X_He` | yes | - | `HeH` at the match (He/H **nuclei**, summed over every carrier) |
| `X_<El>` | no | - | `X_<El>`, hence `melem_ab(iel_<El>)`, for any of the ten element symbols of `species_table`. Overrides `metals.inp`. |
| `F_H`, `F_<El>` | no | g s^-1 | carried; the closure variable of `docs/phase_e_flux_closure_design.md` section 6 |
| `q_H2O`, `q_CO`, ... | no | - | carried; diagnostic, no consumer in the code |

Columns the reader has no consumer for are **kept, not dropped**, and columns
are found by name everywhere: the producers order their element list
differently from run to run.

Interpolation is linear in `log p`, never extrapolated. A target that falls
exactly on a level returns that level's value bit for bit. `K_zz` is
interpolated onto the EXHALE grid through the profile's own radius column
(the cell radius selects the bracketing levels); above the shallowest level
the file carries, a cell takes that level's value, because the eddy
coefficient is a lower-atmosphere property and the file makes no statement
above its top.

How the producers fill the `Kzz` column (`eddy_diffusion_coefficient`,
`src/utils/lower_profile_schema.py`): with no option stated, whatever
`K_zz(p)` the underlying solution carried, the default, in which nothing
changes and the `solution_id` is the one it always was. `--kzz-const K` puts
one value at every level. `--kzz-power ALPHA`, with the amplitude `--kzz-ref
K_ref` and the pressure it is stated at `--kzz-ref-bar p_ref` (default 1 bar),
puts `K_zz(p) = K_ref (p/p_ref)^-ALPHA`, the saturated gravity-wave form;
`ALPHA = 1/2` is the Lindzen (1981) slope Parmentier et al. (2013) fit and 0.4
is the Charnay et al. (2015) one. The two forms are mutually exclusive and
`--kzz-power` without `--kzz-ref` is refused; the power-law settings enter the
`solution_id` fingerprint only when they are in use. On the Photochem arm a
stated coefficient is the one the chemistry is **solved** on and not only the
one written out, so the file never states a composition that no single
`K_zz(p)` produced.

Refusals (all `error stop`): a missing named file; a missing `# columns:`
line; a missing required column; fewer than two levels in the table; a
non-positive or non-decreasing `p`; a missing `solution_id`, `p_match_bar` or
`p_top_bar`; `p_top_bar >= p_match_bar`; a `p_match_bar` outside the table's
coverage; a `p_top_bar` above the shallowest tabulated level (its deep end
needs no separate test, `p_top_bar >= p_match_bar` having already been
refused); a data row that does not carry the number of columns the schema
line names.

Header keys are label-matched with a trailing colon stripped, so `# columns`
and `# columns:` are both accepted; the form written by the adapters, and the
one to write by hand, is `# columns:`.

What the run records: `EXHALE_resolved.out` gains
`lower_profile_present`, and with a profile in use the file name,
`solution_id`, source, `p_match_bar`, `p_top_bar`, the trial fluxes, the
iteration index, `lower_profile_iterable`, `lower_profile_steady`, and -- once
a diffusion step has measured them -- the elemental fluxes over the overlap
window (`lower_profile_F_H_median` / `_spread`, `lower_profile_F_He_median` /
`_spread`, `lower_profile_flux_nface`). `lower_profile_flux_state` is
`unmeasured`, `window_empty` or `measured`; `window_empty` means the escape
window `[j_min:N]` and the profile's coverage do not intersect, so the closure
cannot be measured on that configuration. Alongside the overlap window, and
under the same `lap_flux_measured` guard, the run also reports the same
statistics over the steady-flux window `r >= r_esc`, where the elemental flux
is flat at a steady state: `steady_flux_window_r_lo_Rp`,
`steady_flux_window_nface`, `steady_F_H_median` / `_spread`,
`steady_F_He_median` / `_spread` and `steady_Mdot_median` / `_spread`. That is
the window `src/utils/element_flux_closure.py` falls back on when the overlap
is unmeasurable or too ragged.

With a profile in use the elemental face fluxes are written unconditionally to
`output/element_flux_profile.txt` (the same file the `EXHALE_DIFFUSION_CHECK=1`
diagnostic writes), carrying `F_He` and `F_H` on the same faces from the same
`X` and `J` the step used.

Example: `examples/17_lower_profile/`, whose
`make_example_profile.py` regenerates the synthetic column it ships.

## 3. Parsing semantics

### 3.0 Update (§5.6 Inc 1): the core block is now label-matched

As of §5.6 Inc 1, the Fortran core block is **no longer read by line order**.
Every core key is matched as a **label**: anchored at the start of the
left-trimmed line and terminated by a value separator (`:`, `?`, whitespace, or
`=`). The optional keyword-extension block uses the same anchored matching.
Practical consequences:

- **Line order is irrelevant** and blank / `#`-comment lines are skipped, so a
  reordered header or an inserted comment parses identically. Legacy positional
  files parse unchanged because their lines are self-labeling (e.g. `Planet
  radius [R_J]: 1.401`).
- **The value is still taken by word position** with `get_word(line, n)` on the
  matched line, so every "word N" entry in the §2 table is still accurate; only
  *which* line supplies the value changed (label lookup instead of sequence).
- **Collisions are resolved deterministically.** Anchoring removes the old
  capitalization-only hazard (`Newton solver:` no longer false-matches
  `Solver`). The one prefix pair that anchoring alone cannot separate,
  `Wind-AE seed out` vs `Wind-AE seed` (a whitespace-separated prefix), is
  resolved by testing the longer/most-specific key first in the keyword loop,
  as is `He_metal_diffusion` before `He_diffusion`.
- **A duplicated key is refused** (`refuse_duplicate_keys`, called on the loaded
  lines before any of them is consumed). A key states one quantity, so two lines
  carrying it are two answers to one question and the file does not say which is
  meant. Until 2026-09-03 such a pair resolved to its LAST occurrence -- a rule
  about file order, not about intent -- and a run then proceeded on a value
  nobody had chosen with nothing in its log to say so
  (`backup/regression/valve_sens/eps5` carried `Resid tol` twice and ran on the
  second). The run now stops with `error stop 1`, printing both line numbers and
  both values:

  ```
   (input_read.f90) ERROR: "Resid tol" appears more than once in input.inp:
       line 33: Resid tol: 1.0e-5
       line 36: Resid tol: 1.0e-6
      A key states one quantity, so it may appear only once. Delete the line that is not meant
      (a superseded value belongs in a "#" comment, which is not parsed).
  ```

  Each non-blank, non-`#` line is attributed to the LONGEST key it matches,
  because the key set is not prefix-free (`Lower atmosphere` is a prefix of
  `Lower atmosphere profile`, `Wind-AE seed` of `Wind-AE seed out`); the longest
  match is the key the parser itself consumes. A line matching no key is left to
  the unrecognized-line warning. The same rule and the same routine cover
  `base.inp`; `metals.inp` has its own check keyed on the QUANTITY rather than
  on the token, because there `C` and `CI` both set `X_C`.
- **A missing mandatory core key aborts** with `error stop 1` and a message
  naming the key, instead of silently misreading a neighbor.
- The conditional lines are unchanged in meaning but are now driven by
  **content, not position**: `Spectrum type` selects which property line
  (`Spectrum file` / `Power-law index` / `Photon energy`) is consumed, and
  `Planck` consumes none of them, reading the two optional stellar keys K3
  and K4 instead; the
  energy-band line is read unless the spectrum is monochromatic; the X-ray
  luminosity line is read only when X-rays are included.

The rest of §3 documents the original positional design that Inc 1 replaced.

### 3.1 Core reads (formerly positional; now label-matched)

Before §5.6 Inc 1 the core block was read by a fixed sequence of
`read(11, '(A)') line` statements, each followed by `get_word(line, n)`; the
parser never inspected the label and trusted that line N held the intended
parameter with its value at word position N. That positional design (described
here for reference) had these consequences, now removed by the label matching in
§3.0:

- Inserting or deleting any core line shifted every later core read.
- The number of core lines is not constant: it grows or shrinks with
  `Spectrum type` (line 11a/11b/11c), with the monochromatic flag (line 13
  skipped), and with the X-ray choice (lines 13 word 8 and 14). The Fortran
  parser tracked these with matching `if` branches, so it stayed consistent, but
  any external reader that assumes fixed line numbers will not (see §6).

### 3.2 Keyword matching (historical: the substring design)

*Superseded by §3.0. Kept because the clause ordering it describes is still what
the code does.* The loop tested each line against a chain of
`else if (index(line, 'KEY') > 0)` clauses. Matching was:

- **Substring, not exact.** A line matched a key if the key text appeared
  anywhere in it. This made clause ordering significant when one key is a
  substring of another. The two cases that need it are still ordered the same
  way under anchored matching: `Wind-AE seed out` (K6) before `Wind-AE seed`
  (K7), and `He_metal_diffusion` (K22) before `He_diffusion` (K23).
- **Case-sensitive** for the key text. The `Solver` clause (K28) matches the
  literal `Solver`, while `Newton solver`/`Brent solver`/`Energy solver` use a
  lowercase `solver`. Under the old substring rule the branches stayed separate
  only by that capitalization; anchoring now separates them structurally, and a
  mistyped `Newton Solver:` draws an unknown-line warning instead of
  false-matching.
- **First matching clause wins.** Each line falls into at most one branch (the
  `else if` chain), so a line never fires two keys. This is unchanged.
- **Absent key = compiled default.** Every keyword variable is initialized
  before the loop or in `parameters.f90`, so omitting a keyword line leaves the
  default in place. There is no required *keyword*. The core block, by
  contrast, is mandatory: 19 of its keys are read through `req()`
  unconditionally, and three more (`Spectrum file` / `Power-law index` /
  `Photon energy`, whichever `Spectrum type` selects, none of them for
  `Planck`; the `[E_low` energy-band
  line via `req_eband()`; `Log10 of X-ray luminosity`) are mandatory when their
  condition applies. `Include He23S` is the one core-block line that is
  optional (row 19).

### 3.3 Value-word conventions differ by key

Value word position is not uniform: the `He_*` keys (`He_Kzz`, `He_alphaT`,
`He_ambipolar`, `He_metal_diffusion`, `He_diffusion`, `He_rec_coupling`,
`He_H_charge_exchange`), `CFL`, `Solver`, `Viscosity`, `Conduction`,
`ATES_photoionization_rate`, `Legacy_HHe_rates` and `Secondary_ionization` take
their value at word 2 (short one-word labels), while multi-word labels (`Stellar Lya flux`, `Base BC`, `Domain mode`,
etc.) place the value at word 3, 4, or 5. The table above records each
position. This is a direct consequence of `get_word` counting whitespace tokens
rather than parsing `key: value`.

### 3.4 Order dependencies

- Placement of a keyword line relative to the core block no longer matters
  (with the positional reads it did: a keyword line before `Force start:` was
  consumed as a core line and misread).
- `du_th` (K13) behavior depends on `Reconstruction scheme:` (core line 18):
  two thresholds are honored only when `recon_two_stage` is set by `PLM+WENO3`;
  otherwise the second value is ignored.
- `Reconstruction continuation` (K13b) depends on both `Reconstruction scheme:`
  (core line 18) and `du_th` (K13): it is honored only when `recon_two_stage`
  is set AND the two thresholds leave a non-empty PLM stage
  (`du_th_plm > du_th`). Resolved at the end of `input_read`, after both are
  known; a key that cannot apply is cleared with a WARNING rather than
  silently ignored.
- `IC mode` (K40) and the legacy `Transonic IC`/`Hot Parker IC` (K38/K39)
  interact: the explicit legacy keys take precedence over `IC mode: auto`.

### 3.5 Error behavior

The parser aborts with `error stop 1` (a nonzero exit) on genuine input errors:

- Unknown `Spectrum type`.
- `Jlya escape-prob: True` with `Stellar Lya flux <= 0`.
- `Domain mode: Spherical` without `Outer radius > 1`.
- Molecular chemistry combined with no He.
- `Lower atmosphere` requested without a 1-bar radius or without `EXHALE_ROOT`,
  or a failed generator.
- A key stated more than once, in `input.inp` or in `base.inp`
  (`refuse_duplicate_keys`), or a quantity stated more than once in
  `metals.inp` (`refuse_second_statement`). See the duplicate-key entry in
  §3.0.

Malformed reads inside `read(str, *)` propagate a Fortran I/O error at the point
of the read (the core block does not recover from a malformed line; the companion
files `metals.inp`/`opacity.inp` do skip malformed lines with a warning).

## 4. Appendix A: `metals.inp` schema

Reader: `src/modules/files_IO/metals_input_read.f90` (`read_metals_input`,
called from `input_read`). Any positive element abundance turns metals
on (`thereis_metals`); there is no metals switch in `input.inp`. A missing file
leaves all `X_*` at zero. The abundances are read here but **consumed later**:
since P2, `thereis_metals`, `melem_ab` and `thereis_lowIP_metal` are derived in
the composition block after `read_base_inp`, so a `<El>_H_base` handoff key
(section 2c) overrides this file for that element.

Each quantity may be stated only once (`refuse_second_statement`). The check is
on the QUANTITY, not on the token, because the labels alias: `C` and `CI` both
set `X_C`, so a file carrying both states one abundance twice and the run stops
naming both lines.

Format: one entry on each line, `<token> <value>`, split on the first
whitespace. Blank lines and lines starting with `#` are ignored; malformed lines
warn and are skipped. **Element labels are matched case-sensitively** (so that
`Si` vs `S`, `Na` vs `N`, `Ca` vs `C` stay distinct); each element accepts its
neutral-ion label or its bare symbol.

| Token(s) | Type | Units | Default | Sets | Notes |
|---|---|---|---|---|---|
| `CI` or `C` | real | n_X/n_H | 0.0 | `X_C` | Carbon abundance by number |
| `NI` or `N` | real | n_X/n_H | 0.0 | `X_N` | Nitrogen |
| `OI` or `O` | real | n_X/n_H | 0.0 | `X_O` | Oxygen |
| `MgI` or `Mg` | real | n_X/n_H | 0.0 | `X_Mg` | Magnesium |
| `SiI` or `Si` | real | n_X/n_H | 0.0 | `X_Si` | Silicon |
| `CaI` or `Ca` | real | n_X/n_H | 0.0 | `X_Ca` | Calcium |
| `NaI` or `Na` | real | n_X/n_H | 0.0 | `X_Na` | Sodium |
| `KI` or `K` | real | n_X/n_H | 0.0 | `X_K` | Potassium |
| `SI` or `S` | real | n_X/n_H | 0.0 | `X_S` | Sulfur |
| `FeI` or `Fe` | real | n_X/n_H | 0.0 | `X_Fe` | Iron |
| `cx_full` or `CX_FULL` | 0/1 | - | `.false.` | `cx_full` | Value > 0.5 adds the metal+He (group C) and metal+metal (group D) reactions of Huang 2023 Table 4 to the metal+H set (group A) that is always active. It does **not** gate the He<->H pair (group B): that is on by default under its own `input.inp` key `He_H_charge_exchange` (K14e) and is applied independently of `cx_full`. |
| `cno_cool` or `CNO_COOL` | 0/1 | - | `.true.` | `cno_chianti` | 1 = CHIANTI v11 C/N/O fits (default); 0 = legacy AIOLOS fits (no N cooling). |
| `eos_metals` or `EOS_METALS` | 0/1 | - | `.true.` | `eos_include_metals` | 1 = metals contribute to the mass/electron/particle budget (default); 0 = legacy trace approximation. |
| `pp_metals` or `pp_metal_mode` | 0/1/2 | - | `1` | `pp_metal_mode` | Metal treatment in the advection post-process: 0 metal-free, 1 frozen, 2 re-solve. Out-of-range values are clamped with a warning. |

Any other token warns ("unrecognized ion") and is skipped.

## 5. Appendix B: `opacity.inp` schema

Reader: `src/modules/files_IO/opacity_input_read.f90` (`read_opacity_input`). A
missing file leaves the `parameters.f90` defaults (the legacy analytic cross
sections, model `A`). Format: `KEY = VALUE`, one entry on each line, split on the
first `=`. Blank lines and `#` comments are ignored; malformed lines warn and are
skipped. **Keys are upper-cased before matching** (so key names are
case-insensitive, unlike `metals.inp`).

| Key | Type | Units | Default | Sets | Notes |
|---|---|---|---|---|---|
| `OPACITY_MODEL` | char (first char) | - | `A` | `opacity_model` | Must be one of `A`/`C`/`P`/`T`, else fatal `error stop 1`. The value's first character is upper-cased. |
| `OPA_CONST_HI` | real | factor | 1.0 | `opa_const_HI` | HI cross-section scale factor |
| `OPA_CONST_HEI` | real | factor | 1.0 | `opa_const_HeI` | He I |
| `OPA_CONST_HEII` | real | factor | 1.0 | `opa_const_HeII` | He II |
| `OPA_CONST_HEITR` | real | factor | 1.0 | `opa_const_HeITR` | He I triplet |
| `OPA_PB_FACTOR` | real | - | 0.0 | `opa_pb_factor` | Robinson & Catling pressure-broadening `a` |
| `OPA_PB_EXPONENT` | real | - | 1.0 | `opa_pb_exponent` | Exponent `n` |
| `OPA_PB_PIVOT` | real | dyne/cm^2 | 1.0e5 | `opa_pb_pivot` | Pivot pressure (0.1 bar) |
| `OPA_FILE_HI` | string (path) | - | (unset) | `opa_file_HI` | `.opa` table for HI |
| `OPA_FILE_HEI` | string (path) | - | (unset) | `opa_file_HeI` | He I table |
| `OPA_FILE_HEII` | string (path) | - | (unset) | `opa_file_HeII` | He II table |
| `OPA_FILE_HEITR` | string (path) | - | (unset) | `opa_file_HeITR` | He I triplet table |

Any other key warns ("unknown key") and is skipped.

## 6. Discrepancies and fragilities

Framed tentatively. Items 6.1--6.4 and 6.8 were the observations that motivated
§5.6 Inc 1--3; all five have since been resolved and are kept with their
resolutions, because the resolution is the thing worth knowing.

### 6.1 Positional fragility (resolved)

The core block used to be read by absolute line order, so inserting or removing
any core line shifted every later read; the parser stayed self-consistent only
because it mirrored the conditional lines (spectrum property, energy band,
X-ray) with matching `if` branches. **Resolved by Inc 1:** every core key is
found by anchored label match, order is irrelevant, and a missing mandatory key
aborts with a message naming it (§3.0).

### 6.2 `EXHALE_transit.py` header reads (resolved)

`EXHALE_transit.py` used to read Rp, Mp, T0, a_orb and Mstar by a positional
`num` counter offset by one from the physical line, which would have grabbed the
wrong value had any line been inserted before `Parent star mass`. **Resolved by
Inc 2:** `read_input_params` in `exhale_transit_lib.py` matches every field by
label through `find_input_label`; there is no counter left in the script.

### 6.3 The readers surface different subsets (partly resolved)

- Fortran matches every key by anchored label.
- `exhale_io.py:read_input` is label-based (splits on `:`, matches by string
  prefix) and returns only n0, Rp, Mp, T0, a_orb, r_esc, HeH, Mstar, LEUV, and a
  `spherical` flag, plus the full raw dict. It never parses the energy-band line
  (no `:`), and the X-ray luminosity is present only in the raw dict, not as a
  named field.
- `EXHALE_transit.py` matches by label throughout.

The fragility difference is gone; what remains is that the two Python loaders
surface different *subsets* of the file.

### 6.4 `exhale_io.mdot_log10` and the four flux methods (resolved)

`mdot_log10` used to scale only `Rate/2 + Mdot/2` and `Mdot/4`, with no branch
for `Rate/4 + Mdot` or `alpha`. **Resolved:** `examples/exhale_io.py` now holds
an explicit `_MDOT_FACTOR` table with all four methods (`Rate/4` and `alpha`
take no output correction, by design), and `_mdot_factor()` raises `ValueError`
on anything else.

### 6.5 The GUI writes only the core block

`EXHALE_interface_functions.py:start_func` stops writing at `Force start:`. It
emits no keyword-block line, and its reconstruction list offers only
`PLM`/`WENO3` (not `PLM+WENO3`). Every extension key (du_th, Solver, Stellar
Teff/radius, Domain mode, Base BC, IC mode, and the rest) and any two-stage
reconstruction must be hand-appended. A GUI-produced file is therefore a bare
legacy run; the shipped examples are hand-extended.

### 6.6 GUI energy-band separator differs from the examples

The GUI writes the energy-band line with ` - ` separators (GUI-1007),
whereas every shipped example uses ` , ` separators. Both parse because Fortran
reads only word positions 4/6/8, which land on the numbers regardless of the
separator token. It is a cosmetic inconsistency, not a parse error, but it means
the "canonical" format is ambiguous.

### 6.7 The `Use only EUV?` double negative

`Use only EUV? False` sets `thereis_Xray = .true.` (X-rays INCLUDED), and
`True` leaves X-rays off. The stored flag is the logical negation of
the label, which is easy to misread. This is a semantic gotcha rather than a
bug.

### 6.8 Keyword substring collisions (resolved)

Under `index()` substring matching, safety rested on convention: `Wind-AE seed
out` tested before `Wind-AE seed`, and the `Solver` clause distinguished from
`Newton solver`/`Brent solver`/`Energy solver` purely by capitalization, so a
user writing `Newton Solver:` would false-match plain `Solver`. **Resolved:**
anchored `lbl_match` plus the `known_keys` list, which makes an unrecognized
non-blank line print a warning naming it, added for exactly this typo. The
`Wind-AE seed out` / `He_metal_diffusion` orderings are still needed and still
in place, since those are whitespace-separated prefixes that anchoring alone
does not separate.

### 6.9 Two unrelated "solver" concepts share similar labels

`Solver: Newton` (K28) selects the JFNK marching/steady solver
(`use_newton_solver`), whereas `Newton solver: False` (K41) toggles the
ionization-equilibrium Newton solver (`use_newton_ieq`), and `Brent solver`
(K42) toggles the temperature solver. Three distinct switches with confusingly
similar names; only the code disambiguates them.

### 6.10 Manual vs code

The user manual (`docs/EXHALE_user_manual.tex`, §"input.inp", around lines
466-566) describes the same two-region layout (fixed-order core block plus
optional keyword block) and the conditional spectrum-property line, so it
appears broadly consistent with the code. The manual's keyword table is a
curated subset (it documents `Domain mode`, `Outer radius`, and the headline
physics options) and does not enumerate every one of the 57 keyword keys; this
schema is the complete list. No outright contradiction was found; the manual is
simply less exhaustive than the parser.

### 6.11 Cross-file coupling and inconsistent case conventions

Whether metals are active is decided by the presence of `metals.inp` or of a
`<El>_H_base` handoff key, not by any `input.inp` key, and `base.inp` overrides
`T0`/`R0`/`HeH`/`he_kzz` and the elemental abundances after `input.inp` is read,
as well as setting the base H2 mixing ratio (`q_H2_base`) and the handoff level
(`p_base`). Every override is echoed line by line at startup and the resolved
result is written to `EXHALE_resolved.out`, so "silently" no longer applies;
what remains is that the file is picked up by presence alone. Additionally, `metals.inp` element labels are matched
case-sensitively while `opacity.inp` keys are upper-cased before matching, so the
two companion readers follow opposite case conventions.

---

## 7. Appendix D: the restart files and their `#` headers

`Load IC? True` reads `output/Hydro_ioniz_IC.txt` and
`output/Ion_species_IC.txt`, which are copies of the `Hydro_ioniz.txt` and
`Ion_species.txt` a previous run wrote. Both are plain column files whose `#`
lines are comments; `load_IC` skips them, except for two it reads:

| line | in | read by | meaning |
|---|---|---|---|
| `# columns r[Rp] HI HII ...` | `Ion_species_IC.txt` | `load_IC` | maps species columns by label, so a file with a different column set still restores every element it does carry |
| `# coupling: sec_ion=... sec_ion_step=... recon=... certified=... mode=...` | `Hydro_ioniz_IC.txt` | `parse_coupling_header` | the run state the profile was produced under |
| `# restart_schema 1` and the seven field lines under it | both files | `parse_restart_metadata_line`, `verify_restart_metadata` | the CONFIGURATION the state is a state of (D.2) |

The coupling line is written by `write_coupling_state_header`
(`src/modules/functions/utilities.f90`) and carries every switch that changes
*while* a run converges rather than being fixed by `input.inp`:

| field | value | acted on at restart |
|---|---|---|
| `sec_ion` | `T` / `F`: whether the SvS85/Dalgarno photoelectron secondary ionization was being applied when the file was written | **yes** -- the run starts the coupling there instead of re-staging it |
| `sec_ion_step` | the step it was armed at, `-1` if it was never staged in | reported |
| `recon` | the reconstruction stage in force (`PLM` or `WENO3`) | reported when it disagrees |
| `iontrans` | `T` when the state was produced with the hydrogen ionization state carried (written only then), so its H+ column is a transported quantity and not the local root of its cell | reported, both ways round |
| `certified` | `T` when the A2 stationary certification was made on this state and every active equation was within its tolerance. `cert_reason` stands beside it in two cases: on a refusal it says why (`no_stationary_claim`, `failing_entries`), and on a certified state carrying an element or carrier row it reads `certified_in_wind`, because those rows are gated at `r >= cert_regime_wind_r` and reported below it (decision 22 (a), `docs/certification_tolerance_anchoring_20260910.md`) | read: it is the claim a `Restart intent: stationary` evaluation is held to (its exit status) |
| `mode` | `init` or `phys`: which of the three run states produced the state | read: a `phys` continuation needs a `phys` file with a clock |
| `t_phys` | the elapsed time of a physical state, in seconds; absent for an initialization snapshot | read by a `phys` restart, which continues from it |

`valve` and `fluxconst` left the line with section 152: they were switches of
the OLD base boundary, and the characteristic boundary has neither, so a
header reporting them would report settings the run cannot have. The parser
reads the fields it recognizes one at a time, so an older file still restores
its `sec_ion` state.

`sec_ion` is restored because it is a property of the *state's composition*;
`certified`, `mode` and `t_phys` are properties of the *state*; the rest are
properties of the *run*, which `input.inp` states.
Unknown keys are skipped, so a file written by a later version that carries more
of them is still readable, and a file with **no** `# coupling:` line at all --
anything written before it existed -- restarts exactly as it did before, with
the staging re-run. `write_setup_report` prints which of the three cases the run
is in. Rationale and measurements: `docs/Update_EXHALE_stage1.md` section 144.

The line is a `#` comment, so no numeric parse changes and
`backup/regression/run_check.sh`, which compares with `grep -v '^ *#'`, is
unaffected.

### D.1 The `# provenance:` lines

Every profile file carries three further `#` lines, written by
`write_provenance_header` (same module). Nothing reads them -- they exist so that
a profile found on disk in a year can be tied to an executable, an input and a
set of physics options without asking anyone:

```
# provenance: git=6d07d48afd41 tree=dirty run=2026-09-03T20:19:09
# provenance: ck_input=1260237361560544917 ck_base=-1 ck_metals=4583919659595562347
# provenance: recon=WENO3 base_bc=isothermal carrier=none restart_schema=3 resid_def=145 N=500
```

| field | meaning |
|---|---|
| `git` | the revision the EXECUTABLE was built from, stamped in by the Makefile into a generated `build/build_stamp.f90`. The run never calls `git` |
| `tree` | `dirty` / `clean`: whether the working tree carried uncommitted changes at build time |
| `run` | when THIS run wrote the file, from `date_and_time`. Deliberately not a build time: a timestamp compiled into a constant makes two builds of one source different files, and the md5 of the binary is how a bitwise regression verdict is tied to a revision (section 145.8) |
| `ck_input`, `ck_base`, `ck_metals` | checksums of `input.inp`, `base.inp` and `metals.inp` AS READ, by `file_rolling_checksum` -- two modular rolling sums, exact in `integer(8)` on every compiler, and deliberately not called a digest. `-1` means the file was absent |
| `recon` | the reconstruction in force (`PLM` / `WENO3`) |
| `base_bc` | the lower boundary model. There is one, `characteristic` (section 152 replaced the three component-wise ghost closures the field used to distinguish and retired their keys); the field is kept, and kept in this position, so that a reader of an older profile can still tell which boundary produced it |
| `carrier` | the molecular carrier model: `none` (no molecules), `local` or `transported` |
| `restart_schema` | the version of the restart layout this file follows |
| `resid_def` | the version of the residual definition the acceptance gates used |
| `N` | the grid size |

The regression compares with `grep -v '^ *#'`, so none of this moves a golden.
Rationale: `docs/Update_EXHALE_stage1.md` section 145.7; review section 5, item 12.11.

### D.2 The restart metadata block

Both state files carry eight further `#` lines (plus the option-change
history described at the end of this section), written by
`write_restart_metadata_header` and read by `verify_restart_metadata` (both in
`src/modules/files_IO/load_IC.f90`, so the writer and the loader cannot
disagree about the fields). They state the CONFIGURATION the state is a state
of: a solution reloaded under another composition, grid, constant set or
equation set is not the same problem's state, and until this block existed
nothing in the files said so.

```
# restart_schema 1
# reservoir He/H 8.5099999999999995E-02 C/H 2.6899999999999998E-04 ...
# species_columns 34 r HI HII HeI HeII HeIII HeITR CI CII CIII ...
# grid N 500 R0[cm] 1.5781859000000000E+10 r_min[Rp] 1.0001973003694156E+00 r_max[Rp] 1.5666182019207202E+00 mode Mixed
# constants set IAU2015+CODATA2018 RJ[cm] 7.1492000000000000E+09 kB[erg/K] 1.3806490000000000E-16
# options He23S=T metals=T eos_metals=T mol=F ... jlya=2
# t_phys[s] 0.0000000000000000E+00
# source git=35d9dd5d3ca7 tree=dirty
```

| field | what it states | at restart |
|---|---|---|
| `restart_schema` | the version of the block | a version this build does not compare is refused, rather than compared in part |
| `reservoir` | the elemental composition the state was solved at: `He/H` and every trace element the run carries a reservoir for | **refused** when a ratio differs by more than `heh_dev_tol` (1e-6, the same round-off allowance the He/H check uses) or when an element is in one reservoir and not the other |
| `species_columns` | the column count and the species names, built from `species_table.f90` in the order the `# columns` line uses | reported. Species are restored by their label, and an element the file does not carry is built from the abundance and reported |
| `grid` | the cell count, `R0`, the first and last physical cell center and the construction | **refused** on any difference. Centers from one construction against faces, widths and window indices from another are not one discretization, and a changed Jupiter radius is a changed physical grid; a conservative remap is a separate workflow |
| `constants` | the constant set, with `RJ` and `k_B` written out | **refused** on any difference |
| `options` | the switches that decide which equations the state solves, one token each (`opt_name` in `load_IC.f90`) | compared TOKEN BY TOKEN. **Refused**, with the token named, on any difference the input did not permit with `Restart option change:` (K44); a permitted token may differ and the change is recorded (see below). Numeric inputs (fluxes, tolerances) are NOT here -- they change the coefficients of the same system, and a state carried across such a change is a legitimate starting point whose residual simply is not zero |
| `t_phys[s]` | the clock of a physical state | read only when the `# coupling:` line did not state it, so there is one authority for the clock |
| `source` | the revision the executable was built from | reported. A run whose own initial state came from a pair with no block writes `restart_input=provenance_unknown` here, and the mark is inherited: a state descended from an unknown configuration stays marked |

One caveat for a reader of these files: the field of a header line is its
FIRST token. `species_columns` contains the word `columns`, so a reader that
looks for the label line with a substring test (`'columns' in line`) picks up
the metadata line instead and shifts every column index by one. Every reader
in this tree matches the token; two that did not were fixed when the block
landed (`load_IC.f90`'s `parse_labels` call and
`src/tests/grid_and_gates/restart_round_trip.sh`).

A pair with no `restart_schema` line -- anything written before the block
existed -- is loaded exactly as it was before, and the run says
`provenance unknown` and carries that mark into every file it writes. A pair
whose two halves do not carry the same block was not written by one run and is
refused. `EXHALE_setup.out` states both: whether the loaded state's
configuration could be compared with this run's at all, and which options this
run was allowed to change.

**The option-change history.** A restart that `Restart option change:` (K44)
allowed to change a named option writes one further line into both state
files, after the eight fields:

```
# option_change sec_ion=T -> sec_ion=F at restart of git=35d9dd5d3ca7 tree=dirty
```

Its from and to sides carry only the tokens that actually differed, and the
tail names what produced the state that was loaded. These lines are
informational and never compared. A rung inherits the lines of the state it
was restarted from, oldest first, so the rungs of an arm ladder can be read
back from its last state; the history holds 32 lines, and a history that has
lost older lines carries one line saying how many.

Every number in the block is written with `ES23.16`, seventeen significant
digits, which round trips a double exactly, so a field compared as text is
compared as the number it stands for. The lines are `#` comments, so no
numeric parse and no golden verdict changes.
Design: `docs/restart_contract_design_20260909.md`.

# EXHALE `input.inp` schema (authoritative)

This document is the authoritative schema of EXHALE's `input.inp` file. It began
as §5.6 Inc 0 of `docs/refactor_plan_system_composition_parser.md` ("document the
schema"), the prerequisite for the increments that converted the Fortran
positional reads to anchored label matching and unified the Python loaders.
Those increments have since been carried out; §3.0 describes the parser as it
now behaves, and §3.1--3.2 are retained as a record of the positional design it
replaced.

**Line-number caveat.** The `F<n>` citations throughout point at
`input_read.f90` as it stood when each section was written and have drifted by
tens to a couple of hundred lines. Treat them as navigational hints, not
addresses; the `GUI<n>` citations into
`src/utils/EXHALE_interface_functions.py` are current.

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

- **Fortran** (`input_read.f90`) — the authoritative reader; drives the
  simulation.
- **`examples/exhale_io.py`** (`read_input`) — a label-based loader used by the
  analysis notebooks. It splits each line on the first `:` and matches keys by
  string prefix, so it is robust to line order, but it surfaces only a subset of
  the parameters.
- **`EXHALE_transit.py`** (top level of `EXHALE/`) via
  `exhale_transit_lib.read_input_params` — the transit post-processor. Every
  field it needs (Rp, Mp, T0, a_orb, Mstar, LEUV, `2D approximate method`) is
  matched by label through `find_input_label`; the positional counter it used
  before Inc 2 is gone.
- **`src/utils/EXHALE_interface_functions.py`** — the Tk GUI writer
  (`start_func`). It writes the core block only (through `Force start`); it does
  not emit any keyword-block line.

Two companion configuration files are also read at startup and are covered in
the appendices: `metals.inp` (its mere presence turns metals on) and
`opacity.inp`. A third optional file, `base.inp`, is read by
`read_base_inp` and can override a few core values after `input.inp` is parsed;
it is summarized at the end of the main table notes.

Word positions below refer to whitespace-separated tokens counted from 1, as
`get_word` counts them. "F<n>" cites a line in `input_read.f90`. "GUI<n>" cites
the writing line in `EXHALE_interface_functions.py`. "IO" means
`exhale_io.py:read_input`. "TR" means `EXHALE_transit.py` with its `num`
counter.

## 2. Main table

### 2a. Core block (fixed order; consumed first)

Rows are in the exact order `input_read.f90` reads them.

| # | Line / label | Match mode | Type | Units | Default (if absent) | Read by | Notes |
|---|---|---|---|---|---|---|---|
| 1 | `Planet name:` | positional, word 3 | string | - | mandatory | F49-50; GUI919; IO raw `Planet name` | Sets `p_name`. Value is a single token (word 3); a name with spaces would be truncated. |
| 2 | `Log10 lower boundary number density [cm^-3]:` | positional, word 7 | real | log10(cm^-3) | mandatory | F53-55; GUI923; IO `n0_log10` | `n0` is stored as the exponent, then raised to `10**n0` at F612. Not read by TR. |
| 3 | `Planet radius [R_J]:` | positional, word 4 | real | R_Jupiter | mandatory | F58-60; GUI927; IO `Rp_RJ`; TR num==2 word 4 | Converted to cm (`R0*RJ`) at F652. |
| 4 | `Planet mass [M_J]:` | positional, word 4 | real | M_Jupiter | mandatory | F63-65; GUI931; IO `Mp_MJ`; TR num==3 word 4 | Converted to g (`Mp*MJ`) at F653. |
| 5 | `Equilibrium temperature [K]:` | positional, word 4 | real | K | mandatory | F68-70; GUI935; IO `T0`; TR num==4 word 4 | Base temperature `T0`. Can be overridden later by `base.inp`. |
| 6 | `Orbital distance [AU]:` | positional, word 4 | real | AU | mandatory | F73-75; GUI939; IO `a_AU`; TR num==5 word 4 | Converted to cm (`a_orb*AU`) at F654. |
| 7 | `Escape radius [R_p]:` | positional, word 4 | real | R_planet | mandatory | F78-80; GUI943; IO `r_esc` | Inner escape radius. Not read by TR. |
| 8 | `He/H number ratio:` | positional, word 4 | real | - | mandatory | F83-85; GUI947; IO `HeH` | `HeH > 0` sets `thereis_He`. Can be overridden by `base.inp` and zeroed by a sub-He monochromatic spectrum (F179-185). Not read by TR. |
| 9 | `2D approximate method:` | positional, word 4 (+ word 6 if `alpha`) | keyword string | - | mandatory | F126-139; GUI951; IO raw `2D approximate method`; TR content match (F157) | Word 4 is one of `Mdot/4`, `Rate/2`, `Rate/4`, `alpha`. `Rate/2` is rewritten to `Rate/2 + Mdot/2`, `Rate/4` to `Rate/4 + Mdot` (F138-139). If word 4 is `alpha`, word 6 is read into `a_tau` (F130-132). Sets the day-night flux factor. |
| 10 | `Parent star mass [M_sun]:` | positional, word 5 | real | M_sun | mandatory | F142-144; GUI955; IO `Mstar_Msun`; TR num==9 word 5 | Converted to g (`Mstar*Msun`) at F655. |
| 11 | `Spectrum type:` | positional, word 3 | string | - | mandatory | F147-148; GUI959; IO raw `Spectrum type` | Word 3 must be `Load`, `Power-law`, or `Monochromatic`; anything else is a fatal `error stop 1` (F187-192). The GUI writes `Load from file..`, whose word 3 is `Load`. Selects the conditional line 11a/11b/11c below. |
| 11a | `Spectrum file:` (only if `Load`) | positional, word 3 | string (path) | - | conditional | F154-155; GUI964 | Sets `sed_file`, `do_read_sed = .true.`. |
| 11b | `Power-law index:` (only if `Power-law`) | positional, word 3 | real | - | conditional | F159-161; GUI968 | Sets `PLind`, `is_PL_sed = .true.`. |
| 11c | `Photon energy [eV]:` (only if `Monochromatic`) | positional, word 4 | real | eV | conditional | F170-172; GUI972 | Sets `e_low`, `is_monochr = .true.`. If below the He I threshold, He is removed and `HeH` is zeroed (F179-185). |
| 12 | `Use only EUV?` | positional, word 4 | bool-ish | - | mandatory | F196-198; GUI985/987 | Word 4 `== 'False'` sets `thereis_Xray = .true.` (X-rays INCLUDED). `True` leaves X-rays off. See the double-negative note in §6.7. |
| 13 | `[E_low,E_mid(,E_high)] = [ ... ]` (skipped if monochromatic) | positional, words 4, 6, 8 | reals | eV | conditional | F201-233; GUI999/1004 | If X-rays off: read `e_low` (word 4) and `e_mid` (word 6); `e_top` defaults to `1.24e3`. If X-rays on: also read `e_top` (word 8). The GUI writes `-`-separated values; examples use `,`; both parse because only word positions 4/6/8 matter. |
| 14 | `Log10 of X-ray luminosity [erg/s]:` (only if X-rays on) | positional, word 6 | real | log10(erg/s) | conditional | F237-239; GUI1013 | Sets `LX`; `LX = 0` when X-rays off. |
| 15 | `Log10 of EUV luminosity [erg/s]:` | positional, word 6 | real | log10(erg/s) | mandatory | F245-247; GUI1017; IO `LEUV`; TR content match (F156) | Sets `LEUV`. |
| 16 | `Grid type:` | positional, word 3 | string | - | mandatory | F250-251; GUI1021 | `Uniform`, `Mixed`, or `Stretched`. |
| 17 | `Numerical flux:` | positional, word 3 | string | - | mandatory | F254-255; GUI1025 | `HLLC`, `ROE`, or `LLF`. |
| 18 | `Reconstruction scheme:` | positional, word 3 | string | - | mandatory | F260-267; GUI1029 | `PLM`, `WENO3`, or `PLM+WENO3` (two-stage, `recon_two_stage`). The GUI list offers only `PLM`/`WENO3`; `PLM+WENO3` is hand-added. This choice governs how the later `du_th` values are used. |
| 19 | `Include He23S?` | positional, word 3 | bool-ish | - | **optional**, default on | F251-256; GUI1033/1035 | The one core-block line that may be omitted: `thereis_HeITR` defaults to `.true.` in `parameters.f90`, so a file without the line gets the triplet. Word 3 `True`/`true` sets it, `False`/`false` clears it (the deliberate opt-out, which the HeITR-off regression case uses); any other word leaves the default. Forced off later if He is absent. |
| 20 | `Load IC?` | positional, word 3 | bool-ish | - | mandatory | F275-277; GUI1040/1042 | `True` sets `do_load_IC`. A restart restores every species density the file carries, with one exception: an element whose reservoir a lower-atmosphere handoff states (`<El>_H_base`, or an elemental ratio of the `Lower atmosphere profile:` file) has its loaded column renormalized onto that `El/H` by a single factor common to all its ionization stages, reported in the run log. Elements the handoff does not state are loaded unchanged. |
| 21 | `Do only PP:` | positional, word 4 | bool-ish | - | mandatory | F280-285; GUI1056/1058 | `True` sets `do_only_pp` and clears `force_start`, and arms `sec_ion_active` before the loop so the single equilibrium solve uses the same physics as the post-process (`Update_EXHALE.md` sec. 94). |
| 22 | `Force start:` | positional, word 3 | bool-ish | - | mandatory | F288-293; GUI1062/1064 | `True` sets `force_start` and clears `do_only_pp`. Historically the last positional line; with anchored matching its position no longer matters. |

### 2b. Keyword-extension block (optional; scanned after the core block)

Every line is matched by `lbl_match(line, 'KEY')` in the keyword `do` loop, in
the order shown. Any line may be omitted. None of these are read
by the Python loaders except where noted.

| # | Key substring | Match / value word | Type | Units | Default (if absent) | Sets | Notes |
|---|---|---|---|---|---|---|---|
| K1 | `Domain mode` | word 3 == `Spherical` | flag | - | `.false.` (Roche) | `spherical_domain` | F329-331. Spherical potential out to `Outer radius`. |
| K2 | `Outer radius` | word 4 | real | R_planet | `0.0` | `r_out_user` | F332-334. Required (> 1) in spherical mode; in Roche mode extends `r_max` past L1. |
| K3 | `Stellar Teff` | word 4 | real | K | `0.0` | `T_star_eff` | F335-337. With K4 enables the diluted-blackbody Balmer continuum (`use_excited_H`, F588). |
| K4 | `Stellar radius` | word 4 | real | R_sun | `0.0` | `R_star` | F338-340. Converted to cm (`R_star*Rsun`, F587). |
| K5 | `Deexc heat` | word 3 == `True`/`False` | flag | - | `.true.` | `incl_deexc_heat` | F341-343. Collisional de-excitation heating of H(n=2); active only with the excited-H model (K3+K4). `False` reverts to the one-way coronal ledger. |
| K6 | `Wind-AE seed out` | word 4 | string (path) | - | `''` | `windae_seed_out` | F344-345. Must be tested before K7 (substring). |
| K7 | `Wind-AE seed` | word 3 | string (path) | - | `inputdata/windae_seed.csv` | `windae_seed_file` | F346-347. |
| K8 | `Jlya RT file` | word 4 | string (path) | - | `jlya_rt.txt` | `jlya_rt_file`, `jlya_mode=1` | F348-350. Read a J_Lya(r) profile. |
| K9 | `Jlya escape-prob` | word 3 == `True` | flag | - | `jlya_mode=0` | `jlya_mode=2` | F351-353. Requires `Stellar Lya flux > 0` or fatal `error stop` (F597-605). |
| K10 | `Stellar Lya flux` | word 5 | real | erg/cm^2/s | `0.0` | `F_Lya_star` | F354-356. With K15c on it is also band `B2` of the oxygen photolysis (the 1215.67 A line), so a value set for the H(n=2) pumping drives the photolysis as well. |
| K11 | `Lya stellar halfwidth` | word 5 | real | km/s | `70.0` | `dv_star_lya` | F357-359. |
| K12 | `Lya stellar boost` | word 5 | real | - | `5.0` | `lya_star_boost` | F360-362. |
| K12b | `Lya absorbing bottom` | word 4 == `True`/`False` | flag | - | `.false.` (reflecting bottom) | `lya_bottom_absorber` | Closes the Ly-alpha domain from below with a pure absorber instead of a reflector, as Huang et al. (2017) do: the H2 layer beneath the base really does absorb Ly-alpha through accidental resonances, and the default reflecting bottom keeps every downward photon and over-fills `Jbar` near the base. Adds the downward Neufeld escape `beta_bot`, computed on `tau(base) - tau(r)`, as a channel parallel to the upward escape, so the total escape probability gains a third factor. Acts on the in-line escape-probability RT only (`jlya_mode = 2`, i.e. `Jlya escape-prob: True`); `input_read` warns and ignores it for `jlya_mode` 0 and 1. See `src/modules/radiation/lya_rt.f90`, `docs/lya_destruction_channels.md`. |
| K13 | `du_th` | word 3 (+ optional word 4) | real(s) | - | `du_th=1.0e-3`, `du_th_plm=-1` | `du_th`, `du_th_plm` | F363-383. Two values only if `Reconstruction scheme: PLM+WENO3` (order dependency on line 18). |
| K14 | `ATES_photoionization_rate` | word 2 == `True`/`true` | flag | - | `.false.` (Verner 1996) | `ates_photoion_rate` | F384-388. Reverts He I (1^1S) photoionization to the legacy ATES fit. |
| K14b | `Legacy_HHe_rates` | word 2 == `True`/`true` | flag | - | `.false.` (Badnell/Mao + Voronov) | `legacy_hhe_rates` | Reverts H/He case-B recombination and collisional ionization to the legacy ATES fits (Hui & Gnedin 1997 recombination; Abel+1997/HG97 collisional ionization). Default uses Badnell RR (+ He II DR) minus Mao & Kaastra 2016 alpha_1 for case B, and Voronov 1997 collisional ionization. Free-free always uses the van Hoof et al. 2014 Gaunt table. |
| K14c | `Secondary_ionization` | word 2 == `False`/`True`/`Immediate` | flag | - | `.true.`, STAGED (SvS85 on after first convergence) | `use_sec_ion`, `sec_ion_immediate` | Shull & van Steenberg (1985) secondary ionization by fast photoelectrons (E0 > 40 eV): heating is scaled by f_heat(x) and H I / He I gain secondary ionizations; x is the ionized fraction of the H+He nuclei. STAGED activation (Update §38): the coupling is applied only after the wind first converges without it, then the run re-converges (stops are held N_stall steps after the flip) — from a cold IC the immediate coupling amplifies the base startup transient into a NaN runaway on high-gravity cases. `False` = full photoelectron thermalization (bit-identical to the legacy path); `Immediate` = apply from step 0 (pre-staging behavior, A/B tests only). |
| K14d | `He_rec_coupling` | word 2 == `True`/`False` | flag | - | `.true.` (photons ionize/heat H; `False` = legacy lost-photon path) | `use_he_rec_coupling` | Couples He II -> He I recombination radiation to H I ionization (Draine 2011 on-the-spot y/z; see docs/QUESTIONS_2026-07-17.md). Adds an extra H I photoionization rate `n_HeII n_e [z alpha_B + y alpha_1]` with its photoelectron heating, and corrects the He II recombination to `alpha_B + y alpha_1`. `y` (Eq. 14.16) is the local `>= 24.6 eV` ground-capture share ionizing H; `z` (Sec. 15.5) is the density-dependent cascade share. `False` = pure case B (legacy; in TR mode the singlet recombination is then neither case A nor case B). |
| K14e | `He_H_charge_exchange` | word 2 == `True`/`False` | flag | - | `.true.` (He<->H pair active in every He system) | `he_h_charge_exchange` (charge_exchange module) | The He<->H charge-exchange pair (Huang 2023 Table 4 group B = Koskinen 2013 rates): B1 `He0+H+ -> He+ +H0` (endothermic, exp(-128000/T)) and B2 `He+ +H0 -> He0+H+` (the He+ loss channel where neutral H dominates, ~40% of He II at r~1.05 on WASP-121b, <1% by r>=1.2). Applied by dedicated routines in EVERY ionization system with He (including the advection pair and the analytic Newton Jacobians), independent of `cx_full` (which still gates the metal+He / metal+metal groups C/D). `False` = legacy no-He-CX path (bit-identical). |
| K15 | `Molecular chemistry` | word 3 == `True`/`true` | flag | - | `.false.` | `thereis_mol` | F389-392. Requires He. Metals are allowed: they are solved in the same system as the molecular network. `He_diffusion` is allowed too since 2026-08-26 (milestone M4 of docs/binary_diffusion_design.md): the element transport closes over the molecular carriers. |
| K15b | `Stellar LW flux` | word 5 | real | erg/cm^2/s | `0.0` | `F_LW_star` | Band-integrated stellar flux over 912-1110 A at the planet's orbit. Drives H2 photodissociation in the molecular network (`lyman_werner.f90`): unattenuated rate `1.757e-7 * F_LW` s^-1 times the Draine & Bertoldi (1996) eq. (37) self-shielding factor of the star-ward H2 column, plus 0.4 eV of heating per dissociation. With K15c on it is **also** the first oxygen photolysis band (`LW`): 912-1110 A is one wavelength interval with one incident flux, and H2, H2O and OH absorb out of one beam, so the H2 rate carries the H2O + OH continuum factor of DB96 eq. (40) (identically 1 without K15c) and the H2O/OH rates there carry the fraction of the band the H2 lines have removed. 0 = off (bit-identical to the network without it); a warning is printed if the key is set without `Molecular chemistry`, and another if K15c is on while this is 0 (no H2O/OH photolysis over the interval where their cross sections peak). |
| K15c | `Oxygen chemistry` | word 3 == `True`/`true` | flag | - | `.false.` | `thereis_oxychem` | The A2 option: OH, H2O and CO solved in the coupled molecular ionization equilibrium, with the FUV photolysis of H2O and OH, so the base H2/H partition is computed rather than imported. Requires K15, He/H > 0 and a non-zero oxygen abundance (`metals.inp` `X_O` or `base.inp` `O_H_base`); each is a fatal `error stop`. Refuses `base.inp` `q_H2_base` (it computes that partition) unless a lower-atmosphere profile owns the region below the matching level, and refuses `He_metal_diffusion` (that arm counts an element over its ion stages alone, and CO carries an oxygen *and* a carbon nucleus, so it cannot follow two element factors). `He_diffusion` alone is accepted. Appends `OH H2O CO` to `Ion_species.txt` and writes `output/Oxygen_chemistry.txt` and `output/FUV_bands.txt`. See `docs/a2_oxygen_option_design.md`, `docs/a2_reaction_audit.md`. |
| K15d | `Oxygen transport` | word 3 == `False`/`false` | flag | - | `.true.` (on whenever K15c is on) | `oxygen_transport` | Vertical transport of the molecular carriers H2, OH, H2O and CO: an implicit backward-Euler diffusion-advection step solved together with the same chemistry rows the local equilibrium uses (`diffusive_photochemistry.f90`). Molecular diffusion by Blanc's law over the background carriers, the eddy coefficient from `kzz_cell` (K19 or a profile, and nothing else), settling, no thermal diffusion. Only `False` changes it, and that restores the local steady state of milestone M2 — a test of the chemistry alone, not a model of a base, since at a cool base the H2 chemical time and the flow time are comparable. With `kzz_cell` zero everywhere the transport is pure molecular diffusion and `write_setup_report` warns. Recorded in `EXHALE_resolved.out` as `oxygen_transport`, beside `oxygen_chemistry`, `oxygen_reaction_set`, the five `fuv_band_*_flux` values and `oxygen_base_partition`. Inert without K15c. |
| K15e | `Stellar FUV B1 flux` | word 6 | real | erg/cm^2/s | `0.0` | `F_FUV_B1` | Band-integrated stellar flux over **1110-1201 A** at the planet's orbit (band `B1`), one of the five oxygen photolysis bands. It starts at 1110 A because 912-1110 A is the Lyman-Werner interval and is entered through K15b instead; a value measured over 912-1201 A must be re-split between the two keys. H2O branching OH+H / H2+O(1D) / O+H+H = 0.89/0.11/0.00. Inert without K15c, which warns if a `Stellar FUV B*` key is set with the oxygen chemistry off. |
| K15f | `Stellar FUV B3 flux` | word 6 | real | erg/cm^2/s | `0.0` | `F_FUV_B3` | The same over 1231-1450 A (band `B3`, branching 0.89/0.11/0.00). Band `B2` is the Ly-alpha line and is supplied by K10. |
| K15g | `Stellar FUV B4 flux` | word 6 | real | erg/cm^2/s | `0.0` | `F_FUV_B4` | The same over 1451-2304 A (band `B4`, branching 1.00/0.00/0.00). The H2O cross section falls three decades across this band, so the flat-`F_lambda` band average is off by a factor 4.6-6.1 against a real stellar spectrum; a run in which B4 matters is outside the treatment (`water_photolysis.f90` section 2). |
| K16 | `Molecular base` | word 3 == `True`/`true` | flag | - | `.false.` | `molecular_base` | F393-396. EOS-only molecular base correction to `ntot_bc`. |
| K17 | `Lower atmosphere` | word 3 (+ optional word 4) | string + real | - / R_J | `lower_atm_mode=0` | `lower_atm_mode`, `lower_atm_r1bar` | F397-404. `none`/`analytic`/`vulcan`. Triggers `run_lower_atm_prestep` (needs `EXHALE_ROOT`). |
| K17b | `Lower atmosphere profile` | word 4 | string | file name | `''` (off) | `lap_file` (`lower_atmosphere_profile.f90`) | The lower atmosphere handed over as a **table over an interval of pressure** instead of the single-level scalars of `base.inp` (section 2d). Matched BEFORE K17, whose label is a prefix of this one. With the key set the named file must exist. The profile then owns the base state, the elemental reservoirs and `K_zz`, and the `base.inp` keys of those three categories are **refused** (section 2c). Absent key = present behavior, bit for bit. |
| K18 | `Lower column` | word 3 | real | R_J | `lower_col_r1bar=-1` | `lower_col_r1bar` | F405-407. Analytic lower column 1-bar radius. |
| K19 | `He_Kzz` | word 2 | real | cm^2/s | `0.0` | `he_kzz` | F408-411. The constant eddy diffusion coefficient a run states when it has **no** profile; it fills every entry of `kzz_cell`, which is what the element-diffusion operator reads. Default 0 = pure molecular diffusion (it was `1.0e9` before 2026-08-25). Also overridable by `base.inp` (`Kzz_base`). With a lower-atmosphere profile in use the key is inert -- `K_zz` is then a profile -- and a warning says so. |
| K20 | `He_alphaT` | word 2 | real | - | `0.0` | `he_alphaT` | F412-414. Thermal-diffusion factor. |
| K21 | `He_ambipolar` | word 2 == `False`/`false` | flag | - | `.true.` | `he_ambipolar` | F415-418. Only `False` changes it (default on). |
| K22 | `He_metal_diffusion` | word 2 == `True`/`true` | flag | - | `.false.` | `he_metal_diffusion` | F419-422. Tested before K23. |
| K23 | `He_diffusion` | word 2 == `True`/`true` | flag | - | `.false.` | `he_diffusion` | F423-427. He/H diffusive separation. |
| K24 | `Stall` | words 3, 4 | real, int | - , steps | `stall_tol=1e-6`, `N_stall=2000` | `stall_tol`, `N_stall` | F428-434. Stall-detector override. |
| K25 | `Energy solver` | word 3 == `Explicit` | flag | - | semi-implicit (`.true.`) | `use_semi_implicit_energy=.false.` | F435-442. |
| K26 | `Time stepping` | word 3 == `Local` | flag | - | global (`.false.`) | `use_local_dt=.true.` | F443-450. Cell-by-cell pseudo-time. |
| K27 | `Level tol` | word 3 | real | - | `lev_th=-1` | `lev_th` | F451-455. Mass-flux level-stability tolerance. |
| K28 | `Solver` | word 2 == `Newton` (+ optional word 3) | flag + real | - | `use_newton_solver=.false.`, `newton_du_switch=1e-2` | `use_newton_solver`, `newton_du_switch` | F456-467. JFNK hand-off. Distinct from K41/K42 (see §6.9). |
| K29 | `Valve eps` | word 3 | real | - | `valve_eps=-1` | `valve_eps` | F468-472. Softplus base valve. |
| K30 | `Hydrostatic base` | word 3 == `True` | flag | - | `.false.` | `hydrostatic_base` | F473-477. |
| K31 | `Shapiro filter` | word 3 (+ optional word 4) | real, int | - , steps | `shapiro_eps=-1`, `shapiro_every=4` | `shapiro_eps`, `shapiro_every` | F478-484. |
| K31b | `Low-Mach damping` | word 3 (+ optional word 4) | real, real | - , Mach | `lowmach_damp_eps=-1` (off), `lowmach_damp_mach_th=1e-3` | `lowmach_damp_eps`, `lowmach_damp_mach_th` | Gated fourth-difference (Jameson-Schmidt-Turkel) dissipation of the `2 dr` contact/entropy mode that the contact-resolving HLLC flux stops damping as `v -> 0`. Added to the numerical momentum flux (and its work term to the energy flux) inside `RK_rhs`, so the marching loop and the JFNK steady residual see the same equation -- unlike K31, which touches the marching state only. The gate `[max(0, 1 - M^2/M_th^2)]^2` is exactly zero for `M >= M_th`. Explicit stability needs `eps4 < 1/(16 CFL)`; `input_read` warns otherwise. `eps4 <= 0` = off (default, bit-identical). See `src/modules/flux/low_mach_dissipation.f90`, `docs/hd209_metal_stagnation.md`. |
| K32 | `Base BC` | word 3 (+ optional word 4 if `pressure`) | string + real | - / microbar | `base_bc_mode=0` (density), `base_p_ubar=1.0` | `base_bc_mode`, `base_p_ubar` | F485-498. `density` or `pressure`. Pressure mode derives `n0` (F720-724). |
| K32b | `Base ghost temperature` | word 4 | string | - | `base_ghost_T_continuous=.false.` (isothermal) | `base_ghost_T_continuous` | `isothermal` (legacy, ghost pressure `ntot_bc + dp_bc`, i.e. `T_ghost = T0`) or `continuous` (`dT/dr = 0`, ghost pressure `(ntot_bc + dp_bc)*T_1`). Unknown value warns and keeps isothermal. Ignored when K30 is set. |
| K32c | `Max steps` | word 3 | int | steps | `count_max=1000000` | `count_max` | Hard cap on marching iterations. The env variable `EXHALE_MAXSTEPS` is separate and only exits earlier. |
| K32d | `Coronal cutoff width` | word 4 | real | - | `coronal_cutoff_width=0.1` | `coronal_cutoff_width` | Roll-off width of the coronal-excitation guard below the 1e3 K CHIANTI fit floor (`Cool_coeff.f90`). Must be > 0; input_read aborts otherwise. Since the ground-term fine-structure statistical equilibrium (2026-08-12) the base temperature is insensitive to `w` over 0.02-1.2, so the value is no longer a tuning knob; the earlier 0.08-0.13 justification window is superseded (`docs/coronal_cutoff_width.md` section 7.2). |
| K32e | `Base IR field` | word 4 == `True`/`true` | flag | - | `.false.` | `base_ir_field` | Lets the infrared coolants of a molecular layer see the thermal radiation of the atmosphere below the base instead of emitting into vacuum: the lower atmosphere is taken to be black at those wavelengths and to radiate `B_nu(T0)` over the sky fraction `1 - sqrt(1 - (R_p/r)^2)`. Applies to the eight ground-term fine-structure lines of C I, C II, N II, O I (the incident field enters the statistical equilibrium as a photon occupation number, and the escape probability becomes two-sided) and to the H3+ bands (whose absorption is the Miller et al. 2013 total emission fit evaluated at `T0`, by Kirchhoff's law, since that fit is not a line list). Returns the net rate, emission minus absorption, so each channel stops cooling at its own radiative-equilibrium temperature. Every other channel keeps the optically thin, no-incident-field limit. Off = bit-identical to the previous single-face form; no effect on an atomic run. See `fine_structure_line_transfer` (`Cool_coeff.f90`), `h3p_net_cooling_rate` (`h3p_cooling.f90`), `docs/lower_atmosphere_coupling.md`. |
| K32f | `Molecular IR bands` | word 4 == `True`/`true` | flag | - | `.false.` | `mol_ir_bands` | Adds the infrared coolants a real H2 atmosphere carries below the H2 -> H front: the H2 quadrupole and magnetic dipole line spectrum (Roueff et al. 2019, A&A 630, A58), and the H2O and CO vibration-rotation bands (HITEMP through the Photochem correlated-k coefficients). Each emits in LTE and absorbs the same diluted `B_nu(T0)` K32e supplies, so each stops cooling at its own radiative equilibrium temperature; the closure keeps the stimulated-emission term, so the net rate is exactly zero for gas at `T0` with the whole sky black. Three new NET columns in `output/Cooling_breakdown.txt` (`H2_IR`, `H2O_IR`, `CO_IR`, negative where the band heats) and a column-integrated infrared block in `output/FUV_bands.txt`. H2 needs `Molecular chemistry: True`; H2O and CO need `Oxygen chemistry: True`, and the key is inert on species the run does not carry. With K32e off the bands emit into vacuum, which deepens the collapse rather than holding the layer, and `input_read` warns. Optically thin, LTE, cross sections tabulated 50-2000 K; the Planck-mean optical depths are written to `output/Cooling_breakdown.txt` so the thin assumption is measured. See `src/modules/lower_atmosphere/molecular_infrared_cooling.f90`, `cooling_data/molecular_infrared_bands.py`, `docs/lower_atmosphere_coupling.md`. |
| K33 | `Base velocity` | word 3 | string | - | `base_v_massflux=.false.` (valve) | `base_v_massflux` | F499-508. `valve` or `massflux`. |
| K33b | `Base grid` | words 4 and 5 | real, int | R_planet, cells | `dr_base=2.0e-4`, `N_low_cells=50` | `dr_base`, `N_low_cells` | `Base grid [dr,cells]: <dr_base> [<N_low_cells>]` sets the uniform region of the `Mixed` grid: `N_low_cells` cells of size `dr_base` stacked on the base. The two are not independent — their product is the radial extent of that region (0.01 R_p by default) — so they share one line, as `du_th [PLM,WENO3]` does; refining at fixed extent means dividing the first and multiplying the second (`5.0e-5 200` is the 4x refinement). `dr_base` must resolve the base scale height `H = kT/(mu g)`, which `write_setup_report` echoes as cells per `H`. Ignored by the `Uniform` and `Stretched` grid types. The default is written as a default-real literal, so spelling it out in `input.inp` does **not** reproduce a no-key run bit-for-bit. See `docs/hd189_base_checkerboard.md` §11. |
| K33c | `Grid cells` | word 3 | int | cells | `N=500` | `N` | Number of computational cells of the radial domain; ghost cells are added on top and are not counted. Fewer than 10 is a fatal `error stop 1`. Omitting the key keeps the 500 that used to be a compile-time constant, so an existing `input.inp` is unaffected. Every grid-sized array is allocated by `allocate_grid_arrays` once `N` is known. For the `Mixed` grid the split between the uniform base region and the stretched region is set separately by K33b, and `define_grid` checks the two are compatible. |
| K34 | `Viscosity` | word 2 (+ optional word 3) | `True`/`False` or real, real | - | `visc_on=.false.`, `visc_mu0=0` (off), `visc_s=0.7` | `visc_on`, `visc_mu0`, `visc_s` | `True` = calibrated `mu(T)` + dissipation `q_mu`; a number = diagnostic power law `mu0*T^s` in code units. One-word key, so the value is word 2 (was word 3, which no input file used). See `docs/viscosity_conduction.md`. |
| K34b | `Conduction` | word 2 | `True`/`False` | - | `cond_on=.false.` | `cond_on` | Heat conduction with `kappa(T) = 4.45e4 (T/1000 K)^0.7` (Watson+1981). Independent of K34. |
| K35 | `Resid tol` | word 3 | real | - | `resid_th=-1` | `resid_th` | F516-520. Residual-norm convergence instead of du. |
| K36 | `Resid norm` | word 3 | string | - | `resid_vol=.true.` (vol) | `resid_vol` | F521-528. `vol`/`volume` or `Linf`/`linf`/`LINF`. |
| K37 | `CFL` | word 2 | real | - | `CFL=0.6` | `CFL` | F529-531. |
| K38 | `Transonic IC` | word 3 == `True` | flag | - | `.false.` | `transonic_ic` | F532-534. Forces first iterations (F582-583). |
| K39 | `Hot Parker IC` | word 4 | real | K | `hot_parker_ic=.false.`, `T_wind_ic=1e4` | `T_wind_ic`, `hot_parker_ic` | F535-538. Reads word 4 unconditionally when present. |
| K40 | `IC mode` | word 3 | string | - | `ic_mode=0` (cold) | `ic_mode` (+ `transonic_ic`/`hot_parker_ic`) | F539-566. `cold`/`transonic`/`hot_parker`/`auto`/`windae`. Synonyms for K38/K39; explicit legacy keys take precedence. Unknown value warns and falls back to cold. |
| K41 | `Newton solver` | word 3 == `False` | flag | - | `use_newton_ieq=.true.` | `use_newton_ieq` | F567-569. Ionization-equilibrium Newton solver toggle. |
| K42 | `Brent solver` | word 3 == `False` | flag | - | `use_brent_tsolve=.true.` | `use_brent_tsolve` | F570-572. Temperature Brent solver toggle. |

### 2c. Additional startup file: `base.inp`

Not part of `input.inp`, but read at the same point (`read_base_inp`,
`input_read.f90` F1176-1290) and able to override core values after
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
| `T_base` | EOS boundary | word 2 | `T0` [K] | F1234-1236 |
| `r_base` | EOS boundary | word 2 | `R0` [R_J] | F1237-1239 |
| `p_base` | EOS boundary | word 2 | `p_base_bar` [bar] | F1250-1252. The level every other value refers to. Default `1e-6`; a value away from 1 microbar is echoed at startup, never rejected. |
| `q_H2_base` | EOS boundary | word 2 | `q_h2_base` | F1247-1249. H2 volume mixing ratio at the base. With `Molecular base: True` it replaces the chemical-equilibrium fit in `comp_ntot_bc` (`composition.f90`), i.e. it sets the base **particle count**. It is the **composition of the inflowing gas**, not only an EOS anchor (section 117 of `Update_EXHALE`): it is also imposed on the H2 partition of the lower ghost species, and the base particle count is then taken from that same species state, so the equation of state and the chemistry cannot describe different gas. Refused at startup above the attainable ceiling `0.5/(0.5 + He/H)` rather than capped. Absent (default `-1`) = fit used and the species partition left free, i.e. the historical behavior. |
| `HeH_base` | elemental reservoir | word 2 | `HeH` (sets `thereis_He` if > 0) | F1240-1243 |
| `<El>_H_base` | elemental reservoir | word 2 | `X_<El>`, hence `melem_ab(iel_<El>)` | F1254-1268. `El` is any of the ten element symbols of `species_table` (`C N O Mg Si Ca Na K S Fe`), e.g. `O_H_base 4.90e-4`. Nuclei ratio El/H at the handoff level. **Overrides `metals.inp`** for that element, activates the metal system when it is the only nonzero abundance, and, on a restart, renormalizes that element's loaded column onto this ratio (key 20). Since P2, `thereis_metals`, `melem_ab` and `thereis_lowIP_metal` are all derived after the handoff, so a handoff element behaves exactly like a `metals.inp` element. |
| (none) | initial guess | - | - | No key today. A key here would seed a profile the solver may move away from -- e.g. the A0 improvement of using `q_H2_base` as the base-cell H2 seed in `set_IC`. |
| `Kzz_base` | boundary constraint | word 2 | `he_kzz` [cm^2/s] | F1244-1246. Eddy diffusion coefficient imposed at the base by the element-diffusion operator; inert with `He_diffusion` off. |
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
`K_zz(p)` the underlying solution carried — the default, in which nothing
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
- **A duplicated key resolves to its LAST occurrence** (both the core lookup and
  the keyword loop keep the last match), matching the pre-change keyword-loop
  behavior.
- **A missing mandatory core key aborts** with `error stop 1` and a message
  naming the key, instead of silently misreading a neighbor.
- The conditional lines are unchanged in meaning but are now driven by
  **content, not position**: `Spectrum type` selects which property line
  (`Spectrum file` / `Power-law index` / `Photon energy`) is consumed; the
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
  `Photon energy`, whichever `Spectrum type` selects; the `[E_low` energy-band
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
- `IC mode` (K40) and the legacy `Transonic IC`/`Hot Parker IC` (K38/K39)
  interact: the explicit legacy keys take precedence over `IC mode: auto`.

### 3.5 Error behavior

The parser aborts with `error stop 1` (a nonzero exit) on genuine input errors:

- Unknown `Spectrum type` (F187-192).
- `Jlya escape-prob: True` with `Stellar Lya flux <= 0` (F597-605).
- `Domain mode: Spherical` without `Outer radius > 1` (F661-665).
- Molecular chemistry combined with no He.
- `Lower atmosphere` requested without a 1-bar radius or without `EXHALE_ROOT`,
  or a failed generator (F834-867).

Malformed reads inside `read(str, *)` propagate a Fortran I/O error at the point
of the read (the core block does not recover from a malformed line; the companion
files `metals.inp`/`opacity.inp` do skip malformed lines with a warning).

## 4. Appendix A: `metals.inp` schema

Reader: `src/modules/files_IO/metals_input_read.f90` (`read_metals_input`,
called from `input_read` at F81). Any positive element abundance turns metals
on (`thereis_metals`); there is no metals switch in `input.inp`. A missing file
leaves all `X_*` at zero. The abundances are read here but **consumed later**:
since P2, `thereis_metals`, `melem_ab` and `thereis_lowIP_metal` are derived in
the composition block after `read_base_inp`, so a `<El>_H_base` handoff key
(section 2c) overrides this file for that element.

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

The GUI writes the energy-band line with ` - ` separators (GUI999-1007),
whereas every shipped example uses ` , ` separators. Both parse because Fortran
reads only word positions 4/6/8, which land on the numbers regardless of the
separator token. It is a cosmetic inconsistency, not a parse error, but it means
the "canonical" format is ambiguous.

### 6.7 The `Use only EUV?` double negative

`Use only EUV? False` sets `thereis_Xray = .true.` (X-rays INCLUDED), and
`True` leaves X-rays off (F196-198). The stored flag is the logical negation of
the label, which is easy to misread. This is a semantic gotcha rather than a
bug.

### 6.8 Keyword substring collisions (resolved)

Under `index()` substring matching, safety rested on convention: `Wind-AE seed
out` tested before `Wind-AE seed`, and the `Solver` clause distinguished from
`Newton solver`/`Brent solver`/`Energy solver` purely by capitalization, so a
user writing `Newton Solver:` would false-match plain `Solver`. **Resolved:**
anchored `lbl_match` plus the `known_keys` list, which makes an unrecognized
non-blank line print a warning naming it — added for exactly this typo. The
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

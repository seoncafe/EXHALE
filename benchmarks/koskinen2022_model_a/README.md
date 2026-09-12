# Koskinen et al. (2022) Model A benchmark run

The EXHALE configuration of `docs/koskinen2022_model_a_comparison.tex`: hot
Uranus at 0.05 au, 1 microbar base (`base.inp`), molecular chemistry with H2
and H+ carried by the flow, no Lyman-Werner band, H/He only. Which directory
is which (the H2 comparison of 2026-09-12 was first drawn on the wrong one,
so this file says it explicitly):

* `matched_hnu_minus_I/` -- **the matched run.** Every condition of the
  memo's table in force: `Domain mode: Spherical`, `Outer radius [R_p]: 7.24`
  (Model A's 9.6 R_p top in base radii), `Rate/4 + Mdot`,
  `Ionization transport: True`, `Conduction: True`, `Atomic rate set:
  Koskinen2022`, the Ribas-band spectrum to 2 keV, photoelectron accounting h nu - I. Use this for any
  comparison with Model A. `matched_full_hnu/` and `matched_fraction_0p55/`
  are the same run continued with the two other photoelectron accountings
  (comparison options, not physics).
* `input.inp`, `base.inp`, `output/` at this level -- the **superseded**
  state of the morning of 2026-09-05: the default Roche domain (L1
  truncated, top 4.81 R_p) with the Gueymard-shaped spectrum, the final
  state after 2e5 marching steps (16 threads) continued from
  `output/*_IC.txt`, the 1e5-step transported state of
  `Update_EXHALE_stage1.md` section 169.1 (itself continued from the
  1e6-step `Rate/4` state of section 168); `run.log`,
  `EXHALE_resolved.out`. Made before the gravity and the spectrum were
  matched (memo section "What the matching found"); kept as the record of
  that step, not a Model A comparison. The matched run was started by
  interpolating this state onto the spherical grid.
* `before_section170/` -- that 1e5-step state, written before the
  heating/cooling assembly fix of section 170 (its ghost rows carry the
  section-169.3 writer inconsistency; physical cells are fine).
* `model_a_fig9_digitized.txt` -- their Figure 9 (heating and cooling),
  digitized 2026-09-05. `model_a_fig8_digitized.txt` -- their Figure 8 (H2,
  H, He densities), digitized 2026-09-12 every 0.1 R_p
  (`docs/p23_published_profiles.md` section 5.1 carries the 2026-09-05
  readings of every curve).
* `matched_fixture_states/pass12/` -- the state of the stationary-solver
  fixture re-pinned on THIS configuration
  (`backup/regression/carrier_model_a_newton`, the matched run's state
  mapped onto the current grid) after 12 bounded passes of the partitioned
  loop (memo section 5.4).
* `carrier_reload_states/` -- two states of the stationary-solver fixture
  `backup/regression/carrier_elem_newton` (`pass12/`: 12 bounded passes of
  the partitioned loop; `fixed_wind/`: the fixed point at a fixed wind),
  pinned for the memo's H2-extent section. That fixture is NOT the matched
  configuration (Roche domain to 4.71 r_base, power-law spectrum, 0.048 au,
  1.2 M_sun, `Rate/2 + Mdot/2`, He 2^3S on, no ionization transport, no
  conduction, default atomic rates); the memo says what it is there for.

Regenerate the memo's figures and table with
`python3 docs/k22_model_a_figures.py` (its default input is the matched run;
pass the three `matched_*/output` directories for the three accountings) and
`python3 docs/k22_h2_figure.py`. Every state here is a marching end point,
not a steady-solver root: `du_th` is set so that the run marches for
`EXHALE_MAXSTEPS`; the steady solver would put the composition back on its
local root (section 168).

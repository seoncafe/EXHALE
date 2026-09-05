# Koskinen et al. (2022) Model A benchmark run

The EXHALE configuration of `docs/koskinen2022_model_a_comparison.tex`: hot
Uranus at 0.048 au, 1 microbar base (`base.inp`), `Rate/4 + Mdot`, molecular
chemistry with H2 and H+ carried by the flow, no Lyman-Werner band, H/He only.

* `input.inp`, `base.inp` -- the run configuration
* `output/` -- the final state after 2e5 marching steps (16 threads) continued
  from `output/*_IC.txt`, which is the 1e5-step transported state of
  `Update_EXHALE.md` section 169.1 (itself continued from the 1e6-step
  `Rate/4` state of section 168); `run.log`, `EXHALE_resolved.out`
* `before_section170/` -- that 1e5-step state, written before the
  heating/cooling assembly fix of section 170 (its ghost rows carry the
  section-169.3 writer inconsistency; physical cells are fine)

Regenerate the memo's figures and table with
`python3 docs/k22_model_a_figures.py benchmarks/koskinen2022_model_a/output benchmarks/koskinen2022_model_a/before_section170`.
The state is a marching end point, not a steady-solver root: `du_th` is set
so that the run marches for `EXHALE_MAXSTEPS`; the steady solver would put
the composition back on its local root (section 168).

# p-winds re-retrieval of the LHS 1140 b He 10830 line with the broadening free

This directory re-runs the Cherubim et al. (2026) p-winds retrieval with one
change: the non-thermal broadening of the line, which their analysis holds at
the value the Lampon et al. (2020) turbulence formula returns for the fitted
temperature, is promoted to a free parameter.  The question it answers is
whether their retrieved composition and mass-loss rate survive when the model
is allowed to make the line as wide as the data ask for.

Everything else is theirs.  The forward model, the data, the likelihood, the
SED and its normalization, the system parameters and the priors are taken from
the script the authors released with the paper
(`../../../LHS1140b_zenodo/lhs1140b_2024B_pwinds.py`, Zenodo 15723779).  It is
a stricter reproduction than the earlier oracle in `../pwinds_oracle/`, which
predates that script; the two differ in the turbulence switch and the
five-phase averaging, both of which are used here.

The fit is done in two stages -- match the red-component equivalent width
first, then fit the line shape with the broadening free and the equivalent
width held -- with a joint likelihood grid as a control.  The write-up is
`../../docs/pwinds_broadening_free_refit.md`.

## Files

- `refit_lib.py` — the forward model, the data preparation and the likelihood.
  The one addition to the authors' model is `v_nt`, an isotropic non-thermal
  velocity dispersion (1-D Gaussian sigma) added in quadrature to the Doppler
  width of the line, i.e. to the opacity and not to the emergent spectrum.
  The module docstring records how it is applied and why the implementation is
  exact.
- `verify_reproduction.py` → `verify.log` — validation.  Reproduces the
  authors' own likelihood values at nodes of the grid they released
  (`lhs1140b_grid_likelihood_results_1132_v26e3.pickle`) to
  |delta ln L| < 0.03, and checks that the `v_nt = 0` model is bit-identical
  to calling p-winds with the authors' `turbulence_broadening=True`.
- `two_stage.py` → `two_stage_<grid>_stage{1,2}.txt`,
  `two_stage_<grid>.log` — the two-stage retrieval on two grids:
  `authors_prior` (inside the support of their own grid search) and
  `extended` (composition and temperature carried past it).
- `joint_grid.py` → `joint_grid.npz`, `joint_grid.log` — the control: the
  same likelihood on a joint grid without the equivalent-width constraint,
  with and without `v_nt`, and profile-likelihood intervals from it.
- `summarize.py` → `summary.log`, `fig1_spectrum_fits.pdf`,
  `fig2_broadening.pdf` — comparison tables and the two figures.
- `vnt_scan_published_vector.txt` — chi^2, equivalent width and line width
  against `v_nt` with the other parameters held at the published medians.

## Reproducing

```bash
cd .
OMP_NUM_THREADS=1 /opt/miniconda3/bin/python3 verify_reproduction.py
OMP_NUM_THREADS=1 /opt/miniconda3/bin/python3 two_stage.py authors_prior 64
OMP_NUM_THREADS=1 /opt/miniconda3/bin/python3 two_stage.py extended 64
OMP_NUM_THREADS=1 /opt/miniconda3/bin/python3 joint_grid.py 64
OMP_NUM_THREADS=1 /opt/miniconda3/bin/python3 summarize.py
```

p-winds 1.4.7 (`~/.local/lib/python3.11/site-packages/p_winds`), numpy 1.26.4,
scipy 1.17.1.  Measured on this machine: the two-stage grids take 51 and
59 min on 32 cores each, the joint grid 3.6 min (broadening fixed) plus
43 min (broadening free) on 64 cores.

Note on the equivalent-width solve: the equivalent width is not monotonic in
Mdot, so the measured value is generally reached on two mass-loss branches
that fit the profile very differently.  `solve_mdot` returns both and the
caller keeps the one with the better full-spectrum chi^2; taking only the
lower branch biases stage 2 and makes it composition dependent.

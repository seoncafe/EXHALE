# The base grid of the inputs in this tree

This tree is a preserved record, and its inputs are not edited. The 3
`input.inp` and `input_template.inp` files in it (outside `states/` and
`runs/`) state `Grid type: Mixed` and no `Base grid [dr,cells]:` line, so every
run made from them built its grid on the default of the code at the time:

- uniform base cell width `dr_base = 1.9999999494757503e-4` R_p, the
  single-precision neighbor of 2e-4 (binary64 bits `3F2A36E2E0000000`);
- `N_low_cells = 50` uniform base cells.

From 2026-09-19 the code's default width is the exact double `2.0e-4`
(`dr_base_default = 2.0d-4` in `src/modules/init/parameters.f90`). That builds
another grid: the cell centers move by 3.8e-9 to 6.6e-9 relative on the grids
measured, and `load_IC` refuses a state whose centers differ from the run's by
more than 1e-10. An input of this tree run unchanged on a present binary
therefore builds a grid on which its own stored states do not load.

To reproduce one of these runs, or to restart from one of its states, add this
line to a COPY of its input, next to `Grid type: Mixed`:

```
Base grid [dr,cells]: 1.9999999494757503e-4 50
```

gfortran's list-directed read returns that decimal to the bits above, and the
grid it builds equals the old default grid bit for bit (MEASURED,
`docs/lhs1140b_stationary_D1b_20260919.md`). The inputs of the live
directories were pinned with this line by `src/utils/pin_base_grid.py`.

# `xuvkzz_gm25` — three blocks of the LHS 1140 b memo re-solved on the current Penning coefficient

The He(2^3S)+H Penning coefficient was replaced by the Garcia Munoz (2025)
continuous form (`docs/Update_EXHALE.md` sec. 87) and the advection-corrected
post-process was fixed (sec. 88).  Three places in
`docs/lhs1140b_exhale_vs_pwinds.tex` were left marked as standing on the
retired coefficient.  This directory re-solves them.

| memo | block | arms here |
|---|---|---|
| Sect. `sec:xuvgrid`, Table `tab:xuvgrid` | the XUV grid, both models | `S*` (scalar base, He/H = 2.13), `C*` (flux-closed, He/H = 9.71) |
| Sect. `sec:kzzform`, Table `tab:kzzform` | the `K_zz` pressure form | `P11_*`, `P2_*`, with `L11p1` and `L2p09` as the "recorded, other seed" rows |
| Sect. `sec:outerinherit`, Table `tab:outerinherit` | the outer region a converged wind inherits | `L*`, including `L12p0ctl` and `L12p0jump` |

## What was and was not changed

`resolve_wind.sh` copies a stored arm's configuration and its converged
output, re-solves the wind from that output with four JFNK passes, then runs
one post-process pass and one transit synthesis at R = 68,000.  The
composition, the XUV scaling, `K_zz` and the lower atmosphere are the source
arm's and are untouched; the photochemical columns were not re-solved,
because the Penning coefficient is a wind-side rate.  The source directories
were read only.

`reseed_wind.sh` is the same solve seeded from a different solution, for
bounding the seed dependence of a grid point.

## Files

- `jobs.txt` — the arm-to-source map the batch was driven from.
- `<tag>.log` — the pass record of each arm: `info` and `||R||` per pass, and the steady-state `Mdot`.
- `summarize.py` -> `results.txt` — every number the three blocks need, each re-solved arm printed beside the stored arm it replaces.
- `tex_changes.py` — the list of edits the measurement implies for the memo, each quoting the memo's current text straight out of the file.  It does not edit anything.
- `exobase.py` — the exobase of each ladder arm, for the argument that discards the jump-seeded 12.01 rung.

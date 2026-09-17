# L22 step 2b, 2026-09-16: identities of the measurement

The frozen states reloaded are the three of step 2, unchanged; the md5 of the
files copied into each run directory as `output/Hydro_ioniz_IC.txt` and
`output/Ion_species_IC.txt` (MEASURED at the start of this item, and equal to
the ones `../L22_step2_20260916/IDENTITIES.md` records):

| case | file | md5 |
|---|---|---|
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `output/Hydro_ioniz.txt` | `c35719db84ecc19926f38a4db5718fa2` |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `output/Ion_species.txt` | `f91db8eecc98743b9c82cec59549b182` |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `output/Hydro_ioniz.txt` | `5a9de6ae680246c6badbec424d60b8dd` |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `output/Ion_species.txt` | `b0635219ad30ce5551a066e3b890c90e` |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `output/Hydro_ioniz.txt` | `1229b961323d7ac710146650da72454a` |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `output/Ion_species.txt` | `5c5bfcde1ea939adf39941bdd29002c6` |

| item | value |
|---|---|
| measuring binary | `EXHALE_L22b.x`, md5 `1152b33c11612518dc78a54b11afe5a1` (private build `build_L22b`, deleted at the end of the item) |
| tree binary, untouched | `EXHALE.x`, md5 `c2e9c9990b9f14f1be8cd77abca68945` |
| edited source | `src/modules/lower_atmosphere/diffusive_photochemistry.f90`, md5 `2ac1efac64119da315a4fdccd9c13ff9` before the edit |
| new suite | `src/tests/carrier_outer_boundary/` |
| threads | `OMP_NUM_THREADS=1` in every run |
| run controls | `EXHALE_PTC_DTAU0=1.0e8`, `EXHALE_OUTER_PASSES=3`, `EXHALE_CARRIER_DEBUG=1`, `EXHALE_CARRIER_ROW_TERMS=1` |

The interventions, each one environment key, all off by default:

| key | what it states |
|---|---|
| `EXHALE_L22B_OUTER_DIFF=copy` | the outer diffusive face is formed with a ghost that carries cell N's own particle mixing ratio (A1, zero gradient) |
| `EXHALE_L22B_OUTER_DIFF=extrap` | the same face with the ghost mixing ratio continued linearly from the last two cells, floored at zero (A1) |
| `EXHALE_L22B_OUTER_ADV=upwind` | the advective face at r_{N+1/2} carries the donor cell average and not its reconstruction (A3) |
| `EXHALE_L22B_JAC_RECON=1` | the advective Jacobian entries are a central difference of the face-flux divergence itself, restricted to the tridiagonal band (B) |
| `EXHALE_L22B_JAC_ACTION=<file>`, `EXHALE_L22B_JAC_CELLS=lo,hi` | the action of the assembled Jacobian against a central difference of the full residual, once per run (B) |
| `EXHALE_L22B_DISPLACEMENT=1` | the displacement of p, T, the mean mass per particle and the particle count between the entry and the returned composition of a relaxation pass (C) |

The row dumps in `row_terms/` are the certification of the frozen state under
each setting, `<case>_<setting>_frozen.carrier_row_terms.txt`. The control
dump reproduces the row measure step 2 recorded for the same cell
(`2.46194E-02` at cell 500 of `wm083`; `2.97153E-03` at cell 500 of `kz083`;
`2.47671E-02` at cell 306 of `wm055`), and every intervention was measured on
a state whose carrier densities agree with the control's to every digit
printed (largest relative difference over the 500 cells 0.0, MEASURED), so the
comparison is of operators and not of states.

## The regression comparison, and why it has no control

`backup/regression/mol_carrier` was run on scratch copies, single-threaded,
`EXHALE_MAXSTEPS=12000`, once with `EXHALE.x` and once with `EXHALE_L22b.x`
with none of the keys set. All four output files DIFFER (MEASURED). The
difference is the tree moving under other workers, not this item: a private
build compiles the LIVE tree, and `Hydro_ioniz_adv.txt` now differs in its
schema (`# derived_from:` against `# coupling:` and `# provenance:`), a header
written by `src/modules/files_IO/write_output.f90`, which this item never
opened. Largest relative movement over the 504 written rows: 1.98e-02 to
7.50e-02 on every column of `Hydro_ioniz.txt`, 1.86e-02 to 5.85e-02 on the
atomic stages of `Ion_species.txt`, 1.26e-02 to 2.87e-02 on the molecular
carriers. The two 12000-step marching logs are not kept here; the numbers
above are the record of them.

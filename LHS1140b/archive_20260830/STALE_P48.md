# Stale: transit outputs in this folder include the ghost rows

**stale: computed with ghost rows included (P48, Update_EXHALE_stage1 section 137);
depths move 0.2-3.8%; regeneration awaits instruction.**

Until 2026-09-03 every Python reader of an EXHALE profile returned the file's
**ghost rows** -- the two boundary rows at each end of every
`Hydro_ioniz*.txt` / `Ion_species*.txt`, a fixed base state below and a
zero-gradient / WENO3 extrapolation above -- as if the solver had converged
them. `EXHALE_transit.py` had the same trap independently, so every transit
curve in this folder was integrated over a domain two cells too tall.

## What is stale here

- the 4360 saved `tpm_*` curve and figure files under `LHS1140b/` (subdirectories
  included), and every number read off them;
- the stored outputs of these notebooks:
  - `LHS1140b/Cherubim_2026/LHS1140b_He10833_Zenodo_reproduction.ipynb`
  - `LHS1140b/LHS1140b_analysis.ipynb`

## Size of the shift

Measured on three converged states (section 137.2): line-center depths move by
-3.81 per cent (hot Uranus Ly-alpha, the largest single move) down to less than
0.2 per cent on the weak lines; 4 A band depths by -2.5 to +0.5 per cent. The
direction is almost always shallower, because the two upper ghost rows extended
the chord grid past the last physical cell.

## Not regenerated

Nothing here was regenerated and no value was edited. Re-running
`EXHALE_transit.py` and the notebooks is a separate, instructed step.

# Stale figures: computed with the profile files' ghost rows included

**stale: computed with ghost rows included (P48, Update_EXHALE_stage1 section 137);
depths move 0.2-3.8%; regeneration awaits instruction.**

**The `lhs1140b_*.pdf` set is no longer on this list: it was regenerated
2026-09-16 from `LHS1140b/models/`, the model tree of record on the binary
`c2e9c9990b9f`, by a `make_memo_figures.py` that reads every profile through
`exhale_io.loadtxt_cells`.** Everything else below is unchanged.

Until 2026-09-03 every Python reader of an EXHALE profile returned the file's
**ghost rows** -- the two boundary rows at each end of every
`Hydro_ioniz*.txt` / `Ion_species*.txt`, a fixed base state below and a
zero-gradient / WENO3 extrapolation above -- as if the solver had converged
them. The figures below were made under that behavior. Nothing here was
regenerated and no figure was edited or deleted.

## Transit spectra: the depths themselves move

Line-center depths shift by up to -3.81 per cent (hot Uranus Ly-alpha) and by
less than 0.2 per cent on the weak lines; 4 A band depths by -2.5 to +0.5 per
cent (measured tables: section 137.2).

- `tpm_He10830.pdf`, `tpm_Lya.pdf`, `tpm_Halpha.pdf`, `tpm_Hbeta.pdf` --
  captioned in `docs/EXHALE_user_manual.tex` (Fig. `tpmspec`), marked there.
- `tpm_MgII.pdf`, `tpm_CaII.pdf`, `tpm_NaI.pdf` -- captioned in
  `docs/transmission_spectrum.tex` (Fig. `metalspec`), marked there.
- `tpm_spectra.png`, `transmission_metals.pdf/.png` -- written by
  `examples/make_figures.py`, which drives `EXHALE_transit.py`.
- the `lhs1140b_*.pdf` set -- **regenerated 2026-09-14** and no longer on
  this list; `LHS1140b/make_memo_figures.py` reads through
  `exhale_io.loadtxt_cells` and draws the model tree of record,
  `LHS1140b/models/`.

## Radial profiles: the plotted range, not the physics

These plot the profile arrays directly, so the two ghost rows appear as an
extra point at each end of every curve -- at the base a fixed boundary state,
at the top an extrapolation past the last solution cell. Any quantity read off
the ends of these curves (a top-of-domain value, an axis limit, a base value)
is a boundary value, not a solution.

- `tutorial_overview.pdf/.png`, `metals_on_off.pdf/.png`,
  `cooling_breakdown.pdf/.png`, `spherical_vs_roche.pdf/.png`,
  `newton_convergence.pdf/.png` -- `examples/make_figures.py`.

## Not affected

`cool_*.pdf`, `rec_H.pdf`, `rec_He.pdf`, `he_rec_coupling_*.pdf` and the
`xsec_*.pdf` set are atomic-data and rate-fit figures that read no profile
file.

## Producers

`examples/make_figures.py` reads through `examples/exhale_io.py` and therefore
already carries the fix; `LHS1140b/make_memo_figures.py` was converted to
`exhale_io.loadtxt_cells` on 2026-09-03. Neither was run.

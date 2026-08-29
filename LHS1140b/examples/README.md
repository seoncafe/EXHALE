# LHS 1140 b: the two solutions of record

The two converged winds that carry the memo's headline result
(`docs/lhs1140b_exhale_vs_pwinds.pdf`): the He I 10830 equivalent width of
LHS 1140 b is reproduced at a very different helium abundance depending on
what is placed at the lower boundary, and each case here is solved **at its
own equivalent-width crossing**, so its line matches the measured
`1.108 +/- 0.030 %A` (red pair, vacuum 10832.60-10834.20 A, R = 68,000)
directly rather than by interpolation. Both use the GJ 1132 proxy spectrum
(`../sed/lhs1140_sed_gj1132_at_b.txt`), binary H/He element diffusion at
`K_zz = 1e9 cm2/s`, and the Garcia Munoz (2025) He(2^3S) rate
(`docs/Update_EXHALE.md` section 87).

| case | lower boundary | He/H | red-pair EW [%A] | log10 Mdot |
|---|---|---|---|---|
| `scalar_base_kzz1e9` | prescribed scalar base, metal-free | 1.6261 | 1.10806 | 7.80 |
| `photochem_profile_base` | solved photochemical column (Photochem), reservoir He/H = 9.048 | 9.058 at the matching level | 1.10800 | 7.39 |

The factor ~5.6 between the two abundances is the memo's central point: the
composition this line implies is conditional on the lower boundary, because
the C/N/O the photochemical column carries cool the wind and cut the line at
fixed helium. Neither case is a prediction on its own.

## Contents of each case

- `input.inp` — the runtime configuration, exactly as solved (only the
  spectrum-file path is rewritten relative to this directory).
- `output/` — the converged solution: `Hydro_ioniz.txt`, `Ion_species.txt`,
  the advection-corrected `*_adv.txt` profiles, heating/cooling breakdowns,
  and (diffusion on) `element_flux_profile.txt`.
- `tpm_*.txt` — transmission spectra synthesized from the `_adv` profiles
  (He I 10830 with its line metrics, H-alpha, H-beta, Ly-alpha).
- `EXHALE_setup.out`, `EXHALE_resolved.out` — the configuration echo and the
  resolved provenance record.
- `photochem_profile_base/` additionally carries
  `lower_atmosphere_profile.dat`, the Photochem column handed to EXHALE
  (its `solution_id` is echoed in `EXHALE_resolved.out`). Regenerating that
  file needs the Photochem build described in `../../README_photochem.md`;
  re-running the wind does not.

## Re-running

A run directory is the run: from inside a case,

```bash
../../../EXHALE.x        # reads ./input.inp, writes ./output/
```

re-solves the wind (the stored `output/` doubles as the initial condition
type `Load IC`), and

```bash
MPLBACKEND=Agg PYTHONPATH=../../.. python3 ../../../EXHALE_transit.py
```

re-synthesizes the spectra. For the equivalent width quoted above, set
`EXHALE_TRANSIT_RES_HETR=68000` (see `../winered_hires_y.sh`).

These two directories are the kept representatives of the full run set
(ladders, crossings, diagnosis runs), which is local-only
(`LHS1140b/exhale/`, gitignored) and reproduces from the scripts and
`results.txt` files it carries.

# The LHS 1140 b campaign of 2026-09-14, on the pre-L21 base boundary

These are the results of record of 2026-09-14, moved here unchanged on
2026-09-15 before the campaign was re-run. They are SUPERSEDED. Nothing
here is to be quoted as a current result; it is kept so the re-run can be
compared against it and so every statement written from it (the memo
`docs/lhs1140b_exhale_vs_pwinds.pdf`, `MODELS.md` sections 7 and 8 as they
stood, `docs/Update_EXHALE_stage2.md` section 11) can still be traced to the
files it was written from.

## Why they are superseded

Plan item L21 (`docs/lhs1140b_stationary_L21_20260915.md`) found two defects
in the base boundary condition that decide the state of the first cell on
this planet:

1. the entropy branch was keyed on the CELL-1 VELOCITY, which on LHS 1140 b
   is the odd-even artifact of plan item P44 and not the flow. It selected
   the interior isentrope and put it at the base, so every solution here
   carries a base at the interior entropy instead of the stated reservoir;
2. the wind window was 1e-6, above this planet's own base Mach number of
   4.6e-7, so the branch could not have been decided by the flow even had
   the velocity been clean.

With the branch keyed on the wind-window mass flux and the window at 1e-8,
the base receives the 226 K reservoir of `base.inp` as stated. On the
fiducial case that moves the base temperature from 418 K to 238 K and
log10 Mdot by -1.8 percent. The other regression bases are byte-identical.

## What is here

One directory per case, `<group>/HeH<value>/`, in the same relative layout
as `models/`. Each holds everything the 2026-09-14 run produced --
`output/`, the closure iterations `k00/ ... k05/`, `seed/`, `attempt_*/`,
the logs, the `tpm_*` line products and `REPRODUCE.md` -- together with a
copy of the input files the run was started from (`input.inp`, `base.inp`,
`lower_atmosphere_profile.dat`, `closure.json`, `input_template.inp`), so a
case here is self-contained. It is a record and not a runnable tree: the
runners of `models/` act on `models/`.

83 case directories: the 74 prescribed atomic cases that were solved and the
9 flux-closure rungs. Not here, and why:

- the three cases at 0.01 of the fiducial XUV never solved (item L17), so
  they produced nothing to preserve;
- the nine molecular cases were not part of the 2026-09-14 campaign and stay
  in `models/`. One of them,
  `molecular_scalar_gj1132_kzz1e9/HeH2.13`, was solved on 2026-09-15 by the
  L7f/L7g work on the CORRECTED boundary and on a private binary
  (md5 `48db9e876031`); it is not a pre-L21 result and was left where it is.

`MODELS.md` and `README_models.md` are the catalog and the tree-level
reproduction notes as they stood on 2026-09-14.

## The binaries

Each case's `REPRODUCE.md` names the binary it was solved with and its md5;
they are not all the same binary, the campaign having run across a series of
fixes. The re-run of 2026-09-15 is on one binary throughout, md5
`db87b88d1ce5`.

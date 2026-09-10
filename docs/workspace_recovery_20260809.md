# Workspace loss and recovery, 2026-08-09

- Scope: the whole `RT_Codes/ExoAtmosphere` workspace, not EXHALE alone
- Outcome: everything recovered except `aiolos_0/` and one day's uncommitted work
- Written because the recovery depended on facts about this repository that are
  not visible from the repository itself

## What happened

At 07:22:21 the contents of `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere` were
removed and a clone of `github.com/Nicholaswogan/photochem` was written in their
place. `EXHALE/`, `ATES/`, `aiolos/`, `aiolos_0/` and the workspace `CLAUDE.md`
went with it. The `git reflog` in the new tree reads

```
e1e8725 HEAD@{0}: clone: from https://github.com/Nicholaswogan/photochem
```

and every file carries the clone timestamp. The mtime of the parent `RT_Codes/`
was unchanged, so the directory itself was not recreated: its contents were
emptied and refilled. Only `VULCAN/` and `VULCAN_run_*/`, written 15 s later,
survived in place. A regression run in progress died when its working directory
disappeared.

## What the git remote did and did not protect

`git@github.com:seoncafe/EXHALE.git` was intact at `d4bd507` (2026-07-31 12:05,
"Put exp(-227/T) on the endothermic O + H+ charge-exchange row"). It tracks

```
benchmarks  cooling_data  docs  examples  inputdata  lart_runs
observational_data  paper  poster  python  src  Makefile  *.py  *.sh  README*
```

and nothing else. The following are deliberately untracked and were therefore
**not** recoverable from the remote:

| lost from the remote's point of view | why it matters |
|---|---|
| `WASP-121b/`, `HD209458b/`, `HD189733b/`, `WASP-52b/` | the run folders: `input.inp`, `metals.inp`, converged `output*/`, the analysis notebooks |
| `backup/` | the regression harness and every golden |
| `.gitignore` | never tracked, so the ignore rules themselves were lost |
| `benchmarks/*/output*`, `examples/*/output` | the shipped benchmark and example results |
| `TO_BE_DONE.md`, `LICENCE.md`, `git-update.sh` | ignored by name |

This is the lesson worth keeping: a clone restores the code and the documents,
and none of the run data or the test references.

## How each piece was recovered

Three sources, merged newest-wins.

1. **`git@github.com:seoncafe/EXHALE.git` at `d4bd507`**, all tracked content.
2. **`EXHALE_bkg/`**, a full tree snapshot from 2026-07-17 restored from the
   user's own backup: the original `.gitignore`, the whole `backup/` tree
   (36 case directories, `run_check.sh`, `run_fcheck.sh`, `test_roundtrip.sh`),
   the benchmark and example outputs, `TO_BE_DONE.md`, `LICENCE.md`.
3. **Session scratchpads under `/tmp`**, which sit on local disk and so survived
   an NFS loss: `figs_run/` held the 2026-07-24 converged profiles and
   `input.inp` for the four paper planets, `paper_runs/` and `gate_stage/` the
   staged runs.

`ATES/` and `aiolos/` were restored by the user from separate backups.
`ATES/EXHALE_v1.0/` (2026-07-02) also carries planet-folder notebooks and
`analyze_*.py` scripts and was used where the newer sources had none.

## Verification

Two checks, in order of strength.

**The regression cases reproduce their pre-loss step counts exactly.** Rebuilt
from the restored source at `d4bd507`, with the case inputs taken from
`EXHALE_bkg`, and run single-threaded:

| case | steps | du |
|---|---|---|
| `wasp_full` | 13487 | 9.9645e-04 |
| `wasp_he23off` | 13480 | 9.9501e-04 |

Both match what the pre-loss goldens recorded. Since these runs are long and
path-sensitive, an exact match on both is strong evidence that source, inputs and
build all came back identical.

**The 2026-07-31 work is present in the code, not just in the log.** Checked at
the source rather than by reading commit messages: `charge_exchange.f90` carries
`exp(-227/T)` on the endothermic `O + H+` row (A13, not A14);
`docs/HUANG2023_TABLE4_OXYGEN_ERRATUM.md` is present with its "Status in EXHALE"
section; `paper/ms.tex` holds the corrected paragraph (227.7 K, the
`Stancil1999` citation, the 12.752 helium check) and `paper/ms.bib` both new
entries.

## What is different afterwards

- **`aiolos_0/`**, the pristine AIOLOS reference, was not recovered and has been
  given up. To recover the upstream version of an AIOLOS file, diff against
  upstream directly.
- **One day of uncommitted work** (the He I two-photon correction of
  2026-08-09) was lost and re-applied from the session transcript. It is
  recorded in `Update_EXHALE_stage1.md` §37 and
  `atomic_data_EXHALE_vs_MoCHII.md`. The user confirmed there was no other
  EXHALE work between 2026-07-31 12:05 and the loss.
- **The regression goldens were re-snapshotted** at `d4bd507` before the
  two-photon change was re-applied, so the baseline is the same physics the old
  goldens held, not a new one.
- One trap worth naming: the regression case `input.inp` is **not** the planet
  folder's `input.inp`. It carries `Deexc heat: True` and omits `Solver`, `CFL`
  and `Transonic IC`. A case rebuilt from `WASP-121b/input.inp` is a different
  run, and the goldens taken from it would have been silently wrong. The
  original inputs came back from `EXHALE_bkg`.

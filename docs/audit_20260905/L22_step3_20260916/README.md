# Record of item L22 step 3, increments I3 to I5 (2026-09-17 KST)

The design is `docs/lhs1140b_stationary_L22_step3_design_20260916.md`, approved
by the user with all nine decisions. The memo these files belong to is the
section "Step 3, increments I3 to I5" of
`docs/lhs1140b_stationary_L22_20260916.md`.

| file | what it holds |
|---|---|
| `identities.txt` | the md5 of every source file this increment changed, before and after; the md5 of the private binary and of the control built from the same tree snapshot with those files at their entry text; the md5 of every state the solves were reloaded from |
| `restart_option_change_control.log` | the RED run of `src/tests/grid_and_gates/restart_option_change.sh` against the control binary |
| `restart_option_change_increment.log` | the GREEN run of the same suite against the private binary |
| `suites.txt` | the suites of section 7 of the design, run with the private objects |
| `grid_and_gates.log` | the whole `grid_and_gates` suite, which carries `coupled_carrier_h2` and `restart_option_change` |
| `fixtures.txt` | the three Newton fixtures run with the control binary and with the private one, with the new key value absent |
| `i3_block.log`, `i3_alt.log` | the T2 solves of `molecular_scalar_gj1132_kzz1e9/HeH2.13` reloaded from its certified state, the block against the alternation |
| `i4_*.log` | the T3 solves of the three refusing cases under `Coupled carrier solve: On stall` |

The solves were run on `lart3` under
`LHS1140b/models/.L22/`, which is not a catalog directory: nothing under
`LHS1140b/models/<case>/` was written.

## The hosts

The I3 pair (`i3_block`, `i3_alt`) ran on `lart3`, both started together and
kept side by side: the cost comparison of one route against the other is only
readable if the two carry the same load, and that host carried other work at a
load average of 109 to 128 on 72 cores throughout.

The I4 solves were restarted on `lart4` at 04:30 KST on 2026-09-17, on the
user's instruction that a new long solve goes to the idle host (72 cores, load
average 2.4 at the restart). Their `lart3` logs up to that moment are kept
beside them as `run_lart3.log` in each run directory, and the passes those logs
record are the ones quoted for the `lart3` leg in the memo.

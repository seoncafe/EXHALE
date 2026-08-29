# Path-independence test at He/H = 3.0

The wind alone was re-solved for the converged He/H = 3.0 arm
(`heh3/k03`), keeping that arm's own converged lower-atmosphere profile
and changing only the seed: `heh11p1/k01` instead of the 2.09 chain the
ladder used.  The seed differs by a factor 3.7 in composition and comes
from *above*, a direction the ladder never takes.

Measured on the `_adv` profiles kept here (`r_drop` is the radius where
`T` first falls below half of `T(12 R_p)`):

| | T(12 R_p) [K] | r_drop [R_p] | N(2^3S) | N(2^3S, r > 10 R_p) |
|---|---|---|---|---|
| ladder arm `heh3/k03` | 1034 | 27.63 | 5.6129e+10 | 6.4940e+08 |
| this run (11.1 seed)  | 1043 | 27.63 | 5.7526e+10 | 6.8642e+08 |
| difference            | +0.9 % | identical | +2.5 % | +5.7 % |

`info = 0`, `||R|| = 1.997E-04`.  The outer structure is reproduced, so
this rung does not inherit its outer state from the seed -- unlike the
discarded He/H = 12.01 arm (memo section "The outer region a converged
wind inherits", `docs/Update_EXHALE.md` section 83).

Only the two `_adv` profiles and the input are kept; the run was made in
a scratch directory and the rest was not preserved.  The same test was
not completed at He/H = 5.0, 8.0 or 2.09.

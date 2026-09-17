# The discrete kinetic-energy form of the low-Mach damping stress

What this suite measures, and why it exists: the header of
`src/modules/flux/low_mach_dissipation.f90` argues that the gated
fourth-difference momentum stress and its energy pair remove kinetic energy
with a definite sign, by an integration by parts in which the boundary terms
vanish. The boundary terms do vanish, exactly, because the stress is set to
zero at the faces `j = 0` and `j = N`. What the argument drops is the variable
coefficient: on the actual nonuniform grid the face weight
`w = r_edg^2 eps4 g(M) rho_f lambda_f` runs over decades, and the gate `g`
steps to zero wherever the face Mach number crosses `M_th`.

Writing `d_j = v_{j+1} - v_j`, the discrete kinetic-energy rate of the stress
is exactly `E' = d^T W L d` with `W = diag(w)` and `L` the second difference.
For constant `w` that is `-w sum (d_{j+1} - d_j)^2 <= 0`. For varying `w` the
symmetric part of `W L` is not negative semidefinite, and the suite measures
where it is not.

## What the driver does

`low_mach_stress_energy_form.f90` is standalone: it re-forms the stress from
the equations of the module header and links no production object, so it can
be run while the tree is being edited. It assembles the matrix `A` of
`E' = -v^T A v` over the N physical cells, symmetrizes it and calls LAPACK
`dsyev`. Nonnegative dissipation is `lambda_min >= 0`; the reference on every
row is `-100 eps lambda_max`, the accuracy of the symmetric eigenvalue problem
itself.

Four configurations:

1. `uniform` (N = 64 and N = 500), constant coefficient. The present stress and
   the corrected form must both be dissipative; this is the case in which the
   header's argument is exact.
2. `smooth_coefficient`, a coefficient falling four decades over 500 cells on a
   uniform grid. The present stress must be INDEFINITE, weakly.
3. `gate_edge`, the gate closed over the outer half of a uniform grid, which is
   the sharpest coefficient variation the term's own gate can produce. The
   present stress must be INDEFINITE, and the negative mode must sit on the
   cells astride the gate edge.
4. `state`, the grid and the frozen state of an EXHALE `Hydro_ioniz.txt` given
   as the first argument (default: the LHS 1140 b 0.02-XUV transient of item
   L25). The present stress is reported both on the whole matrix, with a
   zero-gradient ghost closure, and on the block `3..N-2` that no ghost value
   reaches; the corrected form is a verdict.

The corrected form the suite measures alongside is the dissipative operator
`dv/dt = -M^-1 B^T K B v` written as a conservative flux,
`A_j D_j = k_j q_j - k_{j+1} q_{j+1}` with `q_j = 2 v_j - v_{j+1} - v_{j-1}`,
`k >= 0` cell-centered and `k_1 = k_N = 0`. For constant `k` it is
algebraically identical to the present stress; its energy rate is
`-sum_j k_j q_j^2 <= 0` on every grid, with no boundary term and no ghost
value entering. It is measured here and is NOT in the production code.

## Running it

```
src/tests/low_mach_stress_energy/run.sh [Hydro_ioniz.txt]
```

`EXHALE_TEST_OUTDIR` redirects the build directory. The measurements behind
the verdicts are in `docs/lhs1140b_stationary_L25_20260916.md` section 1.3.

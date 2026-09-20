# `h2_level_ladder`: a reduced statistical equilibrium of the H2 ground-state ladder

`run.sh` builds and runs one driver,
`h2_ladder_statistical_equilibrium.f90`, on the model in
`h2_ladder_model.f90`. Both are test material; neither is compiled into
`EXHALE.x`.

## What it is for

`h2_vibrational_heat_fraction` of
`src/modules/lower_atmosphere/h2_vibrational_relaxation.f90` gives the
share of an H2 vibrational excitation that becomes heat rather than
leaving the gas in the infrared quadrupole lines, as

```
f = C1 / (C1 + A_max)
```

with `C1` the v = 1 collisional coefficient of each collider and `A_max`
the largest total spontaneous decay rate over the whole bound ladder.
That is a first-event branching model. It is not a bound over the
ladder: a first-event fraction is not the energy fraction of a multistep
cascade, and the v = 1 collisional coefficient does not bound every
level's collisional loss from below. What it approximates is the net
collisional heat of a solved ladder,

```
Q = sum over level pairs (u, l) and colliders M of
    [ n_u C_ul(M) - n_l C_lu(M) ] n_M (E_u - E_l) ,
```

and this suite solves one and evaluates that expression on it.

## The model, and what is uncertain in it

- **Resolved levels**: the 54 rovibrational levels of X (v <= 3) that
  Lique (2015) distributes state-to-state rate coefficients for.
- **Colliders**: atomic hydrogen (Lique 2015), helium (Jozwiak et al.
  2024, both directions on a 20 to 8000 K grid), and H2 (Le Bourlot,
  Pineau des Forets & Flower 1999). The H2 state-to-state table is not in
  this tree; the helium matrix supplies its shape and it is normalized to
  the published thermal v = 1 -> v' = 0 coefficient the production module
  carries. That is an assumption of the model, and H2 carries 6e-05 of
  the collider sum at the cells it is run on.
- **Radiative rates**: Roueff et al. (2019) table 2, the same file the
  production ladder and infrared line list are reduced from.
- **Upward rates**: from the downward ones by detailed balance, with
  g = g_I (2J+1), the weights of the code's own ladder.
- **One high-level group**: every bound level above the resolved set,
  which is where the nascent molecules of the three-body association R15
  (v = 10 to 14) and of the H3+ dissociative recombination R6 (v = 5 to
  6) are born. No published set reaches those levels. The group carries
  its own energy and empties into the top resolved level at a
  collisional and a radiative rate that are **explicitly uncertain**, and
  the spread of the answer over their stated range is the model's
  uncertainty, together with the chemical destruction rate of an H2
  molecule.

## What it asserts

- The all-level maximum decay rate reduced here is the production
  module's, up to the weak lines the production line list trims.
- **Detailed balance and the thermal limit**: in a collision-only bath
  with no pumping and no radiative escape, the populations are Boltzmann
  and the net collisional heat vanishes. This is the one test that
  isolates the upward rates, which are built and never read.
- **The energy of an injected formation partitions and sums**: gas heat,
  escaping radiation, chemical destruction and stored internal energy
  add up to the injection, in the steady state and over one
  backward-Euler step from an empty ladder.
- **The bracket on f** at three cells of the certified molecular
  fiducial, and the direction of the comparison: the scalar form's
  radiated share is at least the largest the ladder gives, so the scalar
  keeps more energy out of the gas than any member of the bracket and
  replacing it could only add heat.

The numbers are tabulated in
`docs/lhs1140b_stationary_L31_energy_cycles_20260917.md`.

## Running it

```bash
src/tests/h2_level_ladder/run.sh
EXHALE_TEST_OUT=/somewhere/private src/tests/h2_level_ladder/run.sh
```

It needs no build tree: the driver compiles the production sources it
uses through `src/tests/physics_probe/source_closure.py`, and the
statistical equilibrium is one dense LAPACK solve.

The collision tables of Lique (2015) and Jozwiak et al. (2024) are
third-party material distributed in the workspace `references/` tree
beside the repository, not inside it. Their absence is reported as a
SKIP with the path that was looked for.

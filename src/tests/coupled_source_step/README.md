# Coupled temperature-composition source step

Tests of step B3c of `docs/PLAN_20260906_rev2.md`, which implements T1.4 to
T1.9, T2.1 and T2.2 of `docs/b1_target_system_20260906.md`.

Run: `src/tests/coupled_source_step/run.sh`

`EXHALE_OBJDIR` selects another build, `EXHALE_TEST_OUT` where the driver is
written, `EXHALE_COUPLED_EXE` the binary the whole-binary rows use (without it
those rows are skipped).

## What each row asserts

### The ledger identity (AT-1a), `closed_reacting_cell`

One cell, no radiation field, no transport, `heat = cool = 0`, and half the
hydrogen recombining over the step. The reservoir gives that energy back, so
the cell must warm by exactly `du_form/(n_k c_v)` and `u_th + u_form` must not
move. Both are asserted, and the same cell is put through the update the code
had before B3c: there the energy row was anchored on the pressure rebuilt at
the NEW particle count, so with no source at all the step changed the cell's
energy by `X = (n_k(c_new) - n_k(c_old)) u_pp(T_old)`. `X` is measured and
printed beside the temperature each of the two returns, and the row
`the_projection_created_energy_from_nothing` keeps that as the permanent RED
reference.

### The photoevent ledger (AT-1b), `atomic_photoevents` and `h2_photoevents`

One absorbed photon assigns its energy to each recipient exactly once. For H I,
He I, He II, He 2^3S and a metal stage the reservoir share is the difference of
the species formation-energy table and equals the code's own photoionization
threshold, so the two cannot disagree about what an event leaves behind. For
the four H2 channels the four recipients of `h2_channel_energy_recipients` sum
to the photon, and each channel's reservoir share is checked against the same
table.

### The oxygen and associative ledger (AT-1d), `oxygen_and_associative_ledger`

The collisional oxygen channels and the associative He(2^3S) branch are
differences of the one table: a forward and a reverse channel sum to zero, the
eliminated O(1D) carries its excitation to the products of its single sink O6,
and the O4 + O6 pair deposits what the direct H2O photodissociation deposits at
the same photon energy. The associative branch releases 8.066 eV; in an atomic
gas, where the code's own closure sends the HeH+ straight back to He + H, the
cycle returns the metastable's whole 19.82 eV instead, and that is asserted too.

### The reservoir density, `reservoir_density`

`sum_s n_s eps_s` against the sum written out term by term, and the statement
that an excited population is counted on top of its parent column and never
instead of it: doubling the He I column changes nothing, and H(n=2) adds its
excitation to a cell whose H I column is unchanged.

### The whole-binary rows

`the_density_argument_is_intent_in` asserts the declaration in
`ionization_equilibrium.f90`. The property T2.1 asks for is that the density
handed to the chemistry is the density it returns, and that is enforced by the
language rather than by a comparison: no path through the routine can write an
`intent(in)` argument.

The remaining rows run 60 steps of `mol_base_handoff`, the configuration in
which the mass sum has the most ways to disagree with the density (four
molecular carriers, the oxygen columns and the metal block all weigh into
`calc_rho`): every step of the trajectory reaches the fixed point of the
coupled pair, and the species mass sum agrees with the supplied density.

Step 0 is exempt from the fixed-point row and the exemption is named in the
script: at step 0 of a base handoff the composition jumps from the imposed base
to the equilibrium of the state the handoff sets up, which is not a step of the
trajectory, and it is the largest composition change the run ever takes.

### The geometric decay of the pass sequence

The loop's error decays geometrically with a ratio stable from its second
pass, and the code reads that ratio twice: once for the stopping test, which
is taken on the ESTIMATED REMAINING ERROR of the pair, and once for the
extrapolation, which can place the next iterate at the sequence's limit and
is default off.

`accepted_state_within_its_tolerances` is the row the stopping test rests
on, and it is a measurement rather than a model: with
`EXHALE_CSM_ERR_PROBE` the loop holds the state its test accepted, walks the
same sequence on to an increment of 1e-12 (four decades below the
tolerance, so that state is the fixed point at this precision), records the
max-norm distance between the two in each of the two norms of the test, and
then retakes the accepted pass. The row asserts the worst distance over the
run against the two tolerances. `the_probe_does_not_move_the_trajectory`
asserts what makes that measurement belong to the run: the probe writes the
same files as a run without it, which is the statement that the pair
(T, f_sp) it restores is the whole iterate of the loop.

`the_increment_test_is_the_fallback` asserts that the test the loop falls
back to where the sequence gives no estimate is a single path: refusing
every estimate (`EXHALE_CSM_GEOM_REFUSE`) and putting the first pass that
could carry one out of reach (`EXHALE_CSM_GEOM_PASS`) give the same run,
line for line and pass for pass, with every exit taken on the increment.
`the_error_test_lowers_the_pass_count` measures the count against that same
fallback run, and `the_error_estimate_is_one_symbol` reads the source: one
routine forms the ratio, each of the two Rayleigh quotients is assigned in
exactly one place, and both consumers read what it returns.

`refused_extrapolation_costs_nothing` asserts that a candidate refused at
the guards is free: the guards run before any state is moved and before any
pass is spent, so a run in which every candidate is refused has to be the
run with the extrapolation off. `extrapolated_state_is_admissible` asserts
that every constraint on a candidate holds on every candidate taken, which
the ratio test that builds it makes true by construction rather than by a
check after the fact. `extrapolation_does_not_cost_passes` states what the
extrapolation is worth now that the stopping test is taken on the error: the
passes it used to save are the ones that test no longer spends, so the count
with it on must not be above the count without it.

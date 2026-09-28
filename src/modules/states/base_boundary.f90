   module base_boundary
   ! The lower boundary: a characteristic condition imposed AT THE FACE
   ! r_edg(0), with the number of conditions set by the sign and magnitude of
   ! the Mach number there.
   !
   ! WHAT THIS REPLACES.  The previous closure wrote three primitive
   ! components into a ghost CELL with expressions evaluated at ghost-cell
   ! quantities, and Rec_BC then reused the same three numbers as the state AT
   ! THE FACE.  The two locations differ by half a cell of hydrostatic
   ! stratification, dr/(2H), so the face state was systematically too heavy
   ! and too high in pressure.  Dividing that by dr to make a momentum flux
   ! difference leaves an imbalance that is a FIXED FRACTION of the local
   ! weight rho g at every resolution: measured 0.2514, 0.2507, 0.2503, 0.2502
   ! of rho g at cell 1 over a 8x refinement of the hot-Uranus base grid, with
   ! the interior truncation error at 4.9e-6.  The velocity artifact it drives
   ! is first order in dr (-247.1, -123.2, -61.3, -30.5 cm/s over the same
   ! ladder), which is the base sawtooth
   ! reproduced by Euler plus gravity plus that boundary alone, with no
   ! chemistry, no radiation and no wind.  Separately, copying the interior
   ! velocity into the ghost (the one-way valve) made the boundary reflect
   ! 95 percent of an acoustic pulse leaving through it; holding the same
   ! numbers with a stated velocity instead left 3 percent.
   !
   ! THE CONDITION.  At r_edg(0) the eigenvalues are v-c, v and v+c, and a
   ! wave enters the domain when its eigenvalue is positive (the radial
   ! coordinate increases outward, so gas entering at the base moves outward
   ! and inflow is v > 0).  So
   !
   !     0 < v < c   subsonic inflow      2 reservoir conditions, 1 interior
   !    -c < v < 0   flow reversal        1 reservoir condition,  2 interior
   !         v > c   supersonic inflow    3 reservoir conditions, 0 interior
   !         v < -c  supersonic outflow   0 reservoir conditions, 3 interior
   !
   ! THE ORDER OF THE THREE STEPS.  The acoustic matching comes first: the
   ! interior's outgoing relation (C-) and the reservoir's pressure fix the
   ! face pressure and the face velocity.  The direction of the contact is
   ! then read off the matched state, and only then is the contact upwinded
   ! -- the entropy and composition of the face are the reservoir's where
   ! gas enters the domain and the interior's where it leaves.  The
   ! direction is a RESULT of the matching and never an input to it, and it
   ! is never the velocity of the first interior CELL, which at a base
   ! carrying a collocated mode is not the velocity of the face.
   !
   ! The interior relation is the linearized v-c compatibility condition
   !
   !     p_b - rho_i c_i v_b  =  p_i - rho_i c_i v_i ,                    (C-)
   !
   ! the LODI form (Thompson 1987, J. Comput. Phys. 68, 1; Poinsot & Lele
   ! 1992, J. Comput. Phys. 101, 104; as stated in Carlson 2011,
   ! NASA/TM-2011-217181, section 2).  The closed-form invariant
   ! v - 2c/(gamma-1) is NOT used: it assumes one constant gamma, and the
   ! caloric EOS of this code gives every cell its own gamma_eff(x_H2, T), so
   ! that invariant is wrong wherever the H2 fraction or the temperature
   ! varies across the boundary stencil -- which is the whole molecular base.
   ! c_i is the INTERIOR sound speed at the interior gamma_eff, because the
   ! outgoing wave is an interior wave.  An earlier experiment that took the
   ! velocity from v - 2c/(gamma-1) while still pinning rho and p at T0 drove
   ! a spurious base inflow whenever cell 1 was hot; that is what a
   ! compatibility relation closed with an over-stated reservoir pair does.
   !
   ! THE RESERVOIR IS (p, s).  The two reservoir conditions are the pressure
   ! and the specific entropy of the lower atmosphere at the base level, not
   ! the pressure and the density.  The difference matters when the base
   ! COMPOSITION changes: with (p, s) the temperature follows the composition,
   ! with (p, T) the entropy does.  The entropy is carried implicitly, as the
   ! isentrope through (p_base, T_base) at the base composition, which is what
   ! base.inp and the "Lower atmosphere profile:" reader state.
   !
   ! WHO OWNS THE LEVEL WHEN NOTHING FLOWS THROUGH IT.  The reservoir does.
   ! At a pressure-balanced contact at rest the two sides may carry two
   ! entropies and characteristic theory alone does not choose between them;
   ! this model chooses the reservoir, and the physical assumption behind
   ! the choice is stated rather than carried implicitly:
   !
   !     the lower atmosphere below the base level is dense and radiatively
   !     controlled, so on the time scales of a stationary solution it is a
   !     HEAT BATH, and a column standing at rest above it takes the level's
   !     own temperature by radiation and conduction.
   !
   ! That is the assumption that makes T_base a temperature OF THE LEVEL and
   ! not only of the gas that happens to be moving upward through it, which
   ! is what base.inp and the lower-atmosphere profile state it as.  Its
   ! validity is the subsonic base of a level lying inside the radiatively
   ! controlled lower atmosphere (1 microbar for the hot Uranus of Koskinen
   ! et al. 2022, the photochemical column's matching level for LHS 1140 b);
   ! it says nothing about a supersonic base, and nothing about a level
   ! placed above the region where radiation controls the temperature.
   ! During a reverse flow the assumption does not apply to the ADVECTED
   ! trace: gas leaving the domain carries its own entropy and composition
   ! out, and the reservoir then states only the one condition its single
   ! entering characteristic allows, the pressure.  The conduction operator
   ! keeps the bath's temperature at the level in both directions, because
   ! the bath does not stop existing when the gas above it drains; that
   ! condition is written at the operator (viscous_conduction.f90) and not
   ! taken from the advected trace.
   !
   ! THE BASE LEVEL IS THE FACE.  "Base BC: pressure" states the pressure at
   ! r = 1 (the planet radius the handoff quotes: 1 microbar in Koskinen et
   ! al. 2022 section 3.1), r_base_level, and define_grid places the first
   ! face of every grid on that level, r_edg(0) = r_base_level = 1 exactly.
   ! So the reservoir's (p, s) IS the state of the face side it owns, with no
   ! transport along any isentrope, and the contact between the reservoir
   ! and the domain neither moves with the cell width nor encloses a slab of
   ! reservoir gas.  The two lower ghost cells lie below the level and hold
   ! the reservoir continued downward (base_ghost_averages).
   !
   ! Until 2026-09-27 the level was the center of a ghost cell and the face
   ! half a cell above it (Mixed grid); the reservoir was carried to the face
   ! along its own hydrostatic isentrope, which put a slab of 226 K reservoir
   ! gas of thickness dr/2 under the domain and made the base pressure of the
   ! domain p_0 [1 - dr/(2 H_res)]: a first-order error in dr, MEASURED as a
   ! -2.2 % / -1.15 % density offset of the whole LHS 1140 b molecular lower
   ! layer on 500 / 1000 cells (md/ew_grid_order_20260926.md, section 3).
   ! base_boundary_states refuses a grid whose first face is not the level.

   use global_parameters
   use grav_func, only: phi
   use caloric_eos, only: internal_energy_per_particle,                  &
                          heat_capacity_per_particle,                    &
                          adiabatic_index_at_T

   implicit none

   ! ---- the PRESCRIBED reservoir, set once by set_base_reservoir ----
   !
   ! These four numbers are the lower atmosphere's statement about the base
   ! level, and they are the only boundary input the physical column cannot
   ! reconstruct.  They are fixed for the whole run: no solve and no sweep
   ! writes them back (see base_ghost_particle_count below for the count the
   ! sweep measures, which is a different quantity under a similar name).
   real*8 :: base_reservoir_p = 1.0d0   ! p at r_base_level          [p0]
   real*8 :: base_reservoir_T = 1.0d0   ! T at r_base_level          [T0]
   real*8 :: base_reservoir_nhat = 1.0d0! particles per unit mass at the base
   real*8 :: r_base_level = 1.0d0       ! radius the reservoir is stated at

   ! WHICH PRESCRIPTION THE FOUR NUMBERS ABOVE ARE.  A state file carries
   ! them, and a reader has to know what they meant when they were written:
   ! version 1 is (p, T, particles per unit mass, level radius) in the code's
   ! own units, with p = n_hat rho T at the level and the entropy carried
   ! implicitly by the isentrope through (p, T) at the base composition.  A
   ! later version that changes any of those meanings raises the number.
   integer, parameter :: base_reservoir_prescription_version = 1

   ! WHICH BOUNDARY MODEL THE FACE STATE BELOW IS BUILT BY.  One string, so
   ! that a state written by this model and a state written by another are
   ! told apart by their provenance and not by a difference in the numbers.
   ! The model itself is the one this module implements and is described at
   ! the head of the file: the (p, s) reservoir stated at the face, which
   ! is the reservoir level itself (v4; until v3 the face sat dr/2 above
   ! the level and the reservoir was carried to it along its own
   ! hydrostatic isentrope), the linearized C- relation of the
   ! first interior cell, and the contact upwinded on the direction the
   ! matching returns, with the reservoir owning the level at rest.
   !
   ! v2 replaced v1's cubic smoothstep handover of the entropy source, whose
   ! value at zero was the average of the two isentropes.  A state written
   ! under one model and read under the other is the same physical column,
   ! but its residual and its certificate are the other model's; load_IC
   ! says so at the read and loads the state anyway.
   !
   ! v3 adds the ghost's composition to the model.  Under v2 the boundary
   ! took whatever composition the sweep returned for the two lower ghost
   ! cells, and that composition was a function of the one the sweep was
   ! ENTERED at: on one frozen physical column, reservoir, spectrum and
   ! executable, two entry compositions returned ghosts 4.7e-3 to 5.5e-2
   ! apart in their trace ions, both accepted by every test the sweep
   ! applies, and the base continuity row of the first physical cell read
   ! 8.335e-08 in one and 2.923e-09 in the other.  v3 states two things
   ! the state had no way to carry before: WHICH composition the ghost solve
   ! starts from (base_ghost_composition_seed_id below), and that the
   ! composition it returns is a fixed point of that same solve to a stated
   ! accuracy (ionization_equilibrium,
   ! ghost_composition_fixed_point_move and ghost_count_fixed_point_move),
   ! so that the ghost no longer remembers the seed.
   !
   ! v4 puts the face on the level (see "THE BASE LEVEL IS THE FACE" above):
   ! the reservoir is the face state it owns, where v3 carried it from the
   ! level to a face half a cell above along its own isentrope. The grid
   ! moves with it, so a v3 state reaches a v4 run only through
   ! map_state_to_grid.py, as a seed.
   character(len=*), parameter :: base_boundary_model_id =                &
        'characteristic_face_ps_reservoir_C_minus_contact_upwind'//       &
        '_ghost_fixed_point_seed_reservoir_row_face_at_level_v4'

   ! WHICH COMPOSITION THE GHOST SOLVE STARTS FROM.  Part of the model, not
   ! of a reader's handling of rows it drops: the seed selects which of two
   ! states 5 per cent apart in the trace ions the ghost solve returns
   ! whenever the returned ghost is not held to a fixed point, so it is
   ! never neutral and it is stated here.
   !
   ! THE SEED OF RECORD is the composition of the gas the reservoir holds at
   ! the base level, whose partition within each element is the FIRST
   ! PHYSICAL CELL's, at the first sweep of a run, and the ghost the
   ! previous sweep of the run returned at every sweep after it.  The first
   ! is the only seed that is a function of the physical column and the
   ! declared inputs alone, which is why load_IC installs it in place of the
   ! ghost rows a restart pair carries; the second is the state's own ghost,
   ! which is where it exists.  The two are far apart as seeds -- the
   ! reservoir row carries 31 times the ghost's own H II and 26 times less
   ! H3+ -- and the fixed point is what makes the choice cost nothing but
   ! passes.
   !
   ! EXHALE_GHOST_COMPOSITION_SEED replaces it for a measurement, and a run
   ! with that key set is entered at a composition the model does not state.
   character(len=*), parameter :: base_ghost_composition_seed_id =        &
        'the reservoir row at the first sweep, the ghost the previous'//  &
        ' sweep returned after it'

   ! ---- the SOLVED ghost counts: diagnostics of the state, never inputs ----
   !
   ! What the composition sweep measures IN the lower ghost, as opposed to
   ! what the reservoir above prescribes FOR it.  The two are counts of
   ! different quantities AT DIFFERENT RADII: the prescribed count is a
   ! count of nuclei corrected for H2 binding, resolved at startup from the
   ! base mixing ratio AT THE LEVEL r_base_level, and the solved count is
   ! the ghost cell's own heavy-particle count at the composition the sweep
   ! returned, at the ghost's own radius one cell lower.  MEASURED on the
   ! hot-Uranus molecular state of LHS 1140 b: 8.402582811444896e-01 at the
   ! level against 8.794250686178137e-01 in the ghost, 4.7 per cent apart.
   ! That distance is the density ratio across one base cell of that grid,
   ! whose cells are 0.048 pressure scale heights: MEASURED, the same ratio
   ! is 4.87e-02 on the LHS 1140 b He/H 9.7 state and 5.50e-03 on the
   ! hot-Uranus grid of carrier_model_a_newton, and it follows the grid and
   ! not the chemistry.
   !
   ! THE QUANTITY THE TWO DO SHARE is the particles per unit mass, which is
   ! what the ghost construction below continues the level with, and there
   ! the reservoir's prescribed value and the ghost's own stand at 2.7e-08
   ! (hot Uranus) and 1.3e-07 (LHS 1140 b He/H 9.7) (MEASURED, the same
   ! memo).
   !
   ! THEY ARE NOT TWO ESTIMATES OF ONE QUANTITY and the solved one is never
   ! fed back to set_base_reservoir: a count that belongs to the ghost's
   ! radius, restated as the pressure of the level above it, moves the base
   ! pressure by 4.45 per cent and takes the cell-1 continuity row of that
   ! state from 1.2e-08 to 1.0 (the same memo, section 6.2).  They are
   ! reported side by side so that each is read at the radius it belongs
   ! to, and neither is silently substituted for the other.
   real*8  :: base_ghost_particle_count  = 0.0d0  ! heavy particles  [n0]
   real*8  :: base_ghost_electron_count  = 0.0d0  ! electrons        [n0]
   logical :: base_ghost_counts_measured = .false.

   ! THE CLOSURE OF THE GHOST'S OWN H2 PARTITION, as the sweep reports it:
   ! the largest residual of x_H2 - q_H2,base (1 - x_ion) left in a lower
   ! ghost cell, the passes the closure took, and whether it was ever run.
   ! Written by ionization_equilibrium through set_base_ghost_closure.
   real*8  :: base_ghost_closure_residual = 0.0d0
   integer :: base_ghost_closure_passes   = 0
   logical :: base_ghost_closure_measured = .false.

   ! WHAT THE GHOST COMPOSITION FIXED POINT REACHED in the sweep that solved
   ! the ghost: the move of the ghost's species densities over the last
   ! application of its own solve, in the composition's own units (the
   ! species' number density per unit mass of the cell), the move of its
   ! heavy-particle plus electron count over the same application against its
   ! own value (the thermodynamic leg), the applications the sweep took, and
   ! the two accuracies they were held to.
   ! Written by ionization_equilibrium through set_base_ghost_fixed_point.
   ! A measurement of the returned state; no boundary expression reads it.
   real*8  :: base_ghost_fixed_point_move    = 0.0d0
   real*8  :: base_ghost_thermal_move        = 0.0d0
   real*8  :: base_ghost_fixed_point_tol     = 0.0d0
   real*8  :: base_ghost_count_tol           = 0.0d0
   integer :: base_ghost_fixed_point_passes  = 0
   logical :: base_ghost_fixed_point_measured = .false.

   ! ---- the one designed-in width, and what it is a width OF ----
   !
   ! Width, in the interior's Mach number at the face, of the transition
   ! into the SUPERSONIC OUTFLOW branch at M_i = -1, where the count of
   ! characteristics that leave the domain changes from two to three and the
   ! reservoir loses its last condition.  That change of count is a kink in
   ! the residual, and the steady Newton solve needs a differentiable one;
   ! the width is a NUMERICAL REGULARIZATION of that kink and carries no
   ! physics of its own.
   !
   ! It is NOT a width of the entropy branch.  The contact's entropy source
   ! is upwinded on the direction of the mass flux the level carries, which
   ! is a discrete fact about a state and not a quantity that is half true:
   ! a contact with different entropies on its two sides has
   ! direction-dependent traces, and an intermediate face density is an
   ! interpolation between them rather than a law of thermal contact.
   !
   ! THE VALUE STANDS FAR FROM EVERY OPERATING POINT, which is the whole
   ! requirement on it: a base at |M_i| within 1e-8 of unity is a base the
   ! condition has no model for either way, and no state of this code has
   ! been measured inside the transition.  EXHALE_BASE_MACH_BLEND=<value>
   ! overrides it for a control experiment; nothing else reads that variable.
   real*8, save :: base_face_mach_blend = 1.0d-8

   ! WHICH QUANTITY THE DIRECTION OF THE CONTACT IS READ FROM.
   !
   ! The contact moves with the gas, so the side that owns the face is the
   ! side the MASS FLUX at the level comes from.  Three quantities are
   ! candidates for its sign and they are not the same number:
   !
   !   the first interior cell's own product rho_1 v_1 r_1^2.  AT A BASE
   !   THAT CARRIES A COLLOCATED ODD-EVEN VELOCITY MODE THIS IS NOT A FLUX:
   !   MEASURED on the certified LHS 1140 b atomic state, the Riemann face
   !   flux is +1.00000 F_wind at every face of the grid including this one
   !   while the cell-centred product reads -2.00 F_wind at cell 1 and -2.27
   !   at cell 2 (the hot-Uranus base reads -196 F_wind against a face
   !   flux of +0.98).  It is not used.
   !
   !   the matched face velocity v_b of (C-).  This is the face's own
   !   velocity and it is what characteristic theory names, but at the
   !   operating point of this code it does not RESOLVE the direction: v_b
   !   carries (p_res - p_i)/(rho_i c_i), and a converged state leaves that
   !   pressure difference at about 1e-4 of p, which at a base Mach number
   !   of 3e-7 is a velocity two orders above the flow's own.  MEASURED on
   !   the three certified LHS 1140 b states (atomic HeH2.13, molecular
   !   HeH2.13, molecular HeH9.7): v_b = -7.25, -2.10 and -3.41 cm/s while
   !   the level carries the wind's own outward flux and the base face
   !   budget is flat to 1.6e-5, 3.6e-4 and 4.7e-4 of it.  So v_b answers
   !   only where nothing better does.
   !
   !   the mass flux the state carries, read WHERE the cell-centred product
   !   IS a flux: the wind window r >= r_flux that the convergence gate is
   !   already defined on.  In a state whose mass flux is conserved that is
   !   the flux through this face as well, and it is the only one of the
   !   three that resolves the direction on a converged wind.
   !
   ! So the direction is the window's where the window has standing (below),
   ! and the matched face velocity's where it has none -- a cold start, a
   ! column at rest, a state that is not a wind.  The face flux itself is
   ! not available here: it is produced by the Riemann solve this boundary
   ! feeds, and reading it would be the circular dependency this module
   ! exists to break.
   !
   ! LIMITATION, stated rather than hidden: during a transient the base can
   ! carry a flux the wind does not yet, and then this reads the wind's sign
   ! and not the base's.  The branch it selects is still an admissible
   ! boundary condition, and no converged state has the two disagreeing.

   ! WHEN THE WIND WINDOW HAS STANDING TO SPEAK, AND ON WHAT SCALE.
   !
   ! Two conditions, both of them properties of the window's own flux
   ! profile: it must carry a flux at all, and it must carry ONE flux -- a
   ! window that carries a wind carries the same rho v r^2 at every altitude
   ! in it.  The second is a property of the SHAPE of the profile and has no
   ! scale of its own, so it is read as the window's relative spread
   !
   !     d_window = sqrt( <F^2> - <F>^2 ) / |<F>| ,   F = rho v r^2 ,
   !
   ! the relative standard deviation of the flux over r >= r_flux, averaged
   ! over the physical cells of the window.
   !
   ! WHY THE STANDARD DEVIATION AND NOT (max - min)/|mean|.  The range is not
   ! a differentiable function of the state: its two extrema are carried by
   ! particular cells, and where the cell that carries one of them changes
   ! the range keeps its value but changes its slope.  MEASURED on the
   ! certified LHS 1140 b state of `.L14/x003_HeH2.13` with one window cell
   ! swept through the argmax: the two one-sided derivatives
   ! of the face density converge at three step sizes to -242.2 and +0.0128,
   ! a slope jump of 1.9e4 at a point where the value itself is continuous.
   ! The second moment has no argmax: it is a polynomial in the state divided
   ! by the mean, so it is smooth wherever the mean is nonzero.
   !
   ! THE THRESHOLD IS THE MEASURED SEPARATION OF THE TWO POPULATIONS, chosen
   ! the way flux_spread_th_default was (parameters.f90): the geometric middle
   ! of the gap, to one digit, fitted to no single state.  MEASURED 2026-09-17
   ! on d_window, over this window, on the physical cells of every written
   ! state of the two LHS 1140 b catalogs and of the regression fixtures:
   !
   !   certified stationary states, 302 of them
   !                               4.06e-05 to 1.82e-04, median 6.09e-05
   !   the stalled 0.02 transient  1.90e-04
   !   ---------------------------- the gap, a factor 73 wide ---------------
   !   a converged carrier reload  1.39e-02 (carrier_elem_newton)
   !   a partly relaxed reload     1.51e-02 (carrier_model_a_newton's IC)
   !   states that are not winds   4.43e-02 (wasp_he23off), 4.45e-02
   !                               (wasp_full), 1.63e-01 (hp_front),
   !                               1.83e-01 (mol_base_handoff), 2.24e-01
   !                               (carrier_elem_newton's IC), 7.87e-01
   !                               (lower_profile), 8.22e-01
   !                               (atomic_elem_newton), 2.40e+00
   !                               (hydrostatic_column)
   !
   ! sqrt(1.90e-04 * 1.39e-02) = 1.6e-03, and 2e-03 to one digit: 10 times
   ! above the loosest state that is a wind and 7 times below the tightest
   ! state that is not one, so the two populations sit decades away from the
   ! threshold on either side and no state is near it.
   real*8, save :: base_wind_window_spread  = 2.0d-3

   ! THE AMPLITUDE BELOW WHICH THE WINDOW CARRIES NO FLUX AT ALL.
   !
   ! d_window is scale free: it is unchanged when the whole window flux is
   ! multiplied by a constant, so on a window at rest its shape would still
   ! decide the branch while its flux says nothing.  The window therefore
   ! has standing only above a face-mapped Mach number of this size, and
   ! below it the matched face velocity answers.
   !
   ! THE VALUE IS SET SO THAT IT ACTS ON NO STATE THE CODE MEETS.  MEASURED
   ! face-mapped window Mach numbers: 1.13e-08 on the certified 0.03
   ! LHS 1140 b state, 7.49e-09 on the stalled 0.02 transient, 3.75e-06 on a
   ! cold start of the same case, and 3.7e-09 at the
   ! cold start of the wasp_full regression case, which is the smallest any
   ! state has been measured at.  1e-12 stands 3700 times below that and 1e4
   ! below the certified operating point, and the rounding floor of M_wind on
   ! a window of a few hundred cells is sixteen decades below the flux
   ! itself.
   !
   ! EXHALE_BASE_WIND_AMPLITUDE=<value> overrides it and
   ! EXHALE_BASE_WIND_SPREAD=<value> the threshold above, which is how the
   ! sensitivity of the answer to the two scales is measured.
   real*8, save :: base_wind_window_amplitude = 1.0d-12
   logical, save :: base_branch_options_read = .false.


   ! WHICH COMBINATION of the reservoir's (p, s) is held instantaneously.
   !
   !   0 -> the reservoir states its PRESSURE,            p_b = p_res.
   !   1 -> the reservoir states its INCOMING ACOUSTIC INVARIANT at rest,
   !                                                      p_b + rho c v_b = p_res.
   !
   ! THE VALUE IS 0 AND THE ALTERNATIVE IS RECORDED BECAUSE IT WAS MEASURED
   ! AND REJECTED, not because it is a tuning knob.  Both forms are admissible
   ! pairs with the entropy and both are exactly well balanced at rest, where
   ! v_i = 0 makes them identical.  They differ on a MOVING solution, and the
   ! difference is not small:
   !
   !                                weight 0    weight 1
   !   |R|, acoustic reflection      0.956       0.0024
   !   R_mom(1)/rho g, at rest       3.53e-4     1.75e-4
   !   R_mom(1)/rho g, steady wind   1.97e-4     1.32e-1
   !
   ! The last row is the one that decides it.  Driving an EXACT steady
   ! isentropic inflow through the boundary (the manufactured solution of the
   ! design note's test C, base Mach 1e-3), weight 0 leaves cell 1 at 0.68 of
   ! the interior truncation error -- the boundary is as accurate as the
   ! scheme it bounds -- while weight 1 leaves it 460 times ABOVE it.  The
   ! reason is that p_b + rho c v_b = p_res makes the static pressure at the
   ! face fall below the reservoir's by rho c v_b, i.e. by a fraction M of the
   ! pressure, where a real reservoir's Bernoulli drop is rho v^2/2, a
   ! fraction M^2.  Weight 1 over-states it by 1/M and turns the wind's own
   ! mass flux into a base pressure error, which the momentum residual then
   ! multiplies by H/dr.
   !
   ! So the reflection weight 0 leaves standing is NOT a defect of this pair;
   ! it is what any reservoir that holds a STATIC quantity does to sound
   ! (density instead of pressure gives the same |R|), and the ghost closure
   ! this module replaced measured 0.953.  Making the boundary transparent to
   ! sound without corrupting the steady state needs a RELAXATION form, in
   ! which the incoming amplitude is driven toward the reservoir at a finite
   ! rate; that form is not a function of the instantaneous state, so it does
   ! not have a fixed point the steady residual can express.  It is an open
   ! item, not a setting.
   real*8, parameter :: base_incoming_invariant_weight = 0.0d0

   ! A compatibility relation is a linear relation; far from a solution it can
   ! ask for any velocity at all.  The face state is capped at this Mach
   ! number and the events are counted, so a cold start cannot put a
   ! hypersonic ghost into the Riemann solver and no cap can act unreported.
   real*8, parameter :: base_face_mach_max = 5.0d0

   ! ---- diagnostics, reported at the end of a run ----
   ! Evaluations at which the first interior cell was not an admissible gas
   ! state, so the boundary was the reservoir alone (see base_boundary_states).
   integer :: n_base_interior_inadmissible = 0
   integer :: n_base_face_mach_limited   = 0
   integer :: n_base_reversal_evals      = 0
   integer :: n_base_supersonic_outflow  = 0
   real*8  :: base_face_mach_last        = 0.0d0
   real*8  :: base_face_blend_last       = 0.0d0
   ! The quantities the branch is decided by, kept from the last evaluation
   ! so that a run can be asked WHICH branch its base is on and WHY
   ! (report_base_face_state).
   real*8  :: base_face_Mi_last          = 0.0d0
   real*8  :: base_face_Mwind_last       = 0.0d0
   real*8  :: base_face_swind_last       = 0.0d0
   real*8  :: base_face_d_window_last    = 0.0d0
   real*8  :: base_face_vi_last          = 0.0d0
   real*8  :: base_face_vb_last          = 0.0d0
   real*8  :: base_face_rho_res_last     = 0.0d0
   real*8  :: base_face_rho_rev_last     = 0.0d0
   real*8  :: base_face_T_res_last       = 0.0d0
   real*8  :: base_face_T_i_last         = 0.0d0

   ! ---- THE COMPOSITION THE LOWER GHOST CELLS ENTER A SWEEP WITH ----
   !
   ! A MEASUREMENT KEY, DEFAULT OFF, read by nothing but the composition
   ! sweep's own entry state. load_IC replaces the lower ghost rows of a
   ! restart with the composition of the gas the reservoir holds at the base
   ! level, whose partition within each element is the FIRST PHYSICAL
   ! CELL's, and the sweep then solves the ghost's own ionization from that
   ! entry state.
   ! Whether the ghost the sweep RETURNS depends on where it started is a
   ! property of the ghost system and not of the state, and the only way to
   ! measure it is to start the same solve from a stated composition.
   !
   ! EXHALE_GHOST_COMPOSITION_SEED names the source:
   !   unset (the default)         the reservoir row, the rule above
   !   the_state_file_ghost_rows   the lower ghost rows the restart pair
   !                               carries, which load_IC otherwise drops
   !   the_previous_sweep_ghost    the ghost composition the previous sweep
   !                               of this run returned
   !   <path>                      a file of rows in the format
   !                               EXHALE_GHOST_SEED_WRITE writes
   !
   ! VALIDITY: a diagnostic of the ghost solve alone. A run with it set is
   ! entered at a composition the boundary model does not state, so its
   ! result is a measurement and never a certificate. With the key unset
   ! not one number moves: the seed is the entry state's own row.
   character(len=256), save :: ghost_seed_source_env  = ''
   logical, save :: ghost_seed_source_read  = .false.
   real*8,  save :: state_file_ghost_row(1-Ng:0,n_species)     = 0.0d0
   logical, save :: state_file_ghost_row_have    = .false.
   real*8,  save :: previous_sweep_ghost_row(1-Ng:0,n_species) = 0.0d0
   logical, save :: previous_sweep_ghost_row_have = .false.
   real*8,  save :: supplied_ghost_row(1-Ng:0,n_species)       = 0.0d0
   logical, save :: supplied_ghost_row_have      = .false.
   logical, save :: supplied_ghost_row_tried     = .false.

   ! ---- THE GHOST STATE A SWEEP RETURNED, cell by cell ----
   !
   ! Filled by the composition sweep for the two lower ghost cells and read
   ! by the ghost record (boundary_state_trace.f90), which is off unless
   ! EXHALE_GHOST_RECORD is set. Number densities are in cm^-3; the x are
   ! fractions of the element's own nuclei, in the form the ghost's imposed
   ! partition is written in (two hydrogen nuclei per H2 and per H2+, three
   ! per H3+, one per HeH+).
   real*8,  save :: ghost_state_ntot(1-Ng:0)      = 0.0d0
   real*8,  save :: ghost_state_ne(1-Ng:0)        = 0.0d0
   real*8,  save :: ghost_state_x_h2(1-Ng:0)      = 0.0d0
   real*8,  save :: ghost_state_x_h2p(1-Ng:0)     = 0.0d0
   real*8,  save :: ghost_state_x_h3p(1-Ng:0)     = 0.0d0
   real*8,  save :: ghost_state_x_hehp(1-Ng:0)    = 0.0d0
   real*8,  save :: ghost_state_x_hii(1-Ng:0)     = 0.0d0
   real*8,  save :: ghost_state_x_heii(1-Ng:0)    = 0.0d0
   real*8,  save :: ghost_state_x_heiii(1-Ng:0)   = 0.0d0
   ! Charge budget: (sum_i q_i n_i - n_e)/n_e at the returned composition,
   ! zero to rounding wherever the electron count is the composition's own.
   real*8,  save :: ghost_state_charge_gap(1-Ng:0) = 0.0d0
   ! Elemental budget: the largest relative departure, over the tracked
   ! elements, of the ghost's nucleus ratio n_El/n_H from the first physical
   ! cell's. The reservoir row states He/H and each El/H, so this is the
   ! distance between the gas below the base and the gas above it.
   real*8,  save :: ghost_state_element_gap(1-Ng:0) = 0.0d0
   ! The reaction residual the cell was ACCEPTED at, the residual of the
   ! imposed molecular partition at the returned state, and the passes the
   ! closure took.
   real*8,  save :: ghost_state_reaction_res(1-Ng:0)  = 0.0d0
   real*8,  save :: ghost_state_partition_res(1-Ng:0) = 0.0d0
   integer, save :: ghost_state_closure_passes(1-Ng:0) = 0
   logical, save :: ghost_state_recorded = .false.

   contains

   !------------------------------------------!

   subroutine report_base_face_state(tag)
   ! What the base face condition did at its last evaluation: which branch,
   ! with what weight, and the two densities the weight mixes.  Written
   ! where a diagnostic asks for it; the run itself does not print it.
   character(len=*), intent(in) :: tag
   write(*,'(A)') ' [base face] '//trim(tag)
   write(*,'(A,ES13.6,A,ES13.6)') '   interior at the face: v_i =',        &
        base_face_vi_last*v0, ' cm/s,  M_i =', base_face_Mi_last
   write(*,'(A,ES13.6,A,ES13.6)') '   the wind window: M_wind =',         &
        base_face_Mwind_last, ',  its own relative flux spread d =',       &
        base_face_d_window_last
   write(*,'(A,F4.1,A,ES10.3,A,ES10.3,A)')                                &
        '     does the window have standing in the direction? ',          &
        base_face_swind_last, '  (1 where |M_wind| >=',                    &
        base_wind_window_amplitude, ' and d <=',                          &
        base_wind_window_spread, '; else the face velocity decides)'
   write(*,'(A,ES13.6,A,ES13.6)') '   face state:           v_b =',        &
        base_face_vb_last*v0, ' cm/s,  M_b =', base_face_mach_last
   write(*,'(A,F4.1,A)') '   contact upwind w_rev =',                      &
        base_face_blend_last, '  (0 = the reservoir states the entropy'//  &
        ' of the level,'
   write(*,'(A)') '     1 = the interior trace leaves through the face)'
   write(*,'(A,ES13.6,A,ES13.6,A,F8.5)') '   rho_res =',                   &
        base_face_rho_res_last*n0, '  rho_rev =',                          &
        base_face_rho_rev_last*n0, '  [mH/cm3], ratio ',                   &
        base_face_rho_rev_last/max(base_face_rho_res_last,1.0d-99)
   write(*,'(A,F10.3,A,F10.3,A)') '   T_res =',                            &
        base_face_T_res_last*T0, ' K   T_i =', base_face_T_i_last*T0,      &
        ' K  (the reservoir''s and the interior''s, at the face)'
   end subroutine report_base_face_state

   !------------------------------------------!

   subroutine set_base_reservoir(p_level, T_level, nhat_level, r_level)
   ! Record the lower atmosphere's (p, T, composition) at the base level.
   ! The entropy is not stored as a number: it is carried by the isentrope
   ! through this state, which continue_hydrostatic_isentrope follows.
   !
   ! CALLED ONCE, at startup, with the prescription the input states. The
   ! counts a composition sweep measures in the ghost are NOT this: they go
   ! to set_base_ghost_counts below and stay there.
   real*8, intent(in) :: p_level, T_level, nhat_level, r_level
   base_reservoir_p    = p_level
   base_reservoir_T    = T_level
   base_reservoir_nhat = nhat_level
   r_base_level        = r_level
   end subroutine set_base_reservoir

   !------------------------------------------!

   subroutine set_base_ghost_counts(n_heavy, n_electron)
   ! The heavy-particle and electron counts of the lower ghost at the
   ! composition a sweep returned, in units of n0. A measurement of the
   ! state, kept for the report and for the consistency the run can be asked
   ! about; no boundary expression reads it.
   real*8, intent(in) :: n_heavy, n_electron
   base_ghost_particle_count  = n_heavy
   base_ghost_electron_count  = n_electron
   base_ghost_counts_measured = .true.
   end subroutine set_base_ghost_counts

   !------------------------------------------!

   subroutine set_base_ghost_closure(residual, passes)
   ! What the ghost's H2 partition closed to, from the sweep that solved it.
   real*8,  intent(in) :: residual
   integer, intent(in) :: passes
   base_ghost_closure_residual = residual
   base_ghost_closure_passes   = passes
   base_ghost_closure_measured = .true.
   end subroutine set_base_ghost_closure

   !------------------------------------------!

   subroutine set_base_ghost_fixed_point(move, thermal_move, passes,      &
                                         move_tol, count_tol)
   ! What the ghost composition fixed point reached, from the sweep that
   ! solved the ghost. move_tol is in the composition's own units and
   ! count_tol against the count's own value.
   real*8,  intent(in) :: move, thermal_move, move_tol, count_tol
   integer, intent(in) :: passes
   base_ghost_fixed_point_move     = move
   base_ghost_thermal_move         = thermal_move
   base_ghost_fixed_point_passes   = passes
   base_ghost_fixed_point_tol      = move_tol
   base_ghost_count_tol            = count_tol
   base_ghost_fixed_point_measured = .true.
   end subroutine set_base_ghost_fixed_point

   !------------------------------------------!

   function ghost_composition_seed_source() result(src)
   ! EXHALE_GHOST_COMPOSITION_SEED, read once for the run and announced only
   ! when it is set, so that a run without the key is the run without the
   ! code. An empty string is the default rule (the entry state's own row).
   character(len=256) :: src
   if (.not. ghost_seed_source_read) then
      call get_environment_variable('EXHALE_GHOST_COMPOSITION_SEED',      &
                                    ghost_seed_source_env)
      ghost_seed_source_read = .true.
      if (len_trim(ghost_seed_source_env) .gt. 0)                         &
         write(*,'(A)') ' (base_boundary) EXHALE_GHOST_COMPOSITION_'//    &
              'SEED: the lower ghost cells enter the composition sweep'// &
              ' at '//trim(ghost_seed_source_env)
   endif
   src = ghost_seed_source_env
   end function ghost_composition_seed_source

   !------------------------------------------!

   logical function ghost_composition_seed_armed() result(on)
   character(len=256) :: src
   src = ghost_composition_seed_source()
   on  = (len_trim(src) .gt. 0)
   end function ghost_composition_seed_armed

   !------------------------------------------!

   subroutine set_state_file_ghost_rows(f_rows)
   ! The lower ghost rows a restart pair carries, kept before load_IC
   ! replaces them with the reservoir row. Nothing reads them unless the
   ! seed key names them.
   real*8, intent(in) :: f_rows(1-Ng:0,n_species)
   state_file_ghost_row      = f_rows
   state_file_ghost_row_have = .true.
   end subroutine set_state_file_ghost_rows

   !------------------------------------------!

   subroutine set_previous_sweep_ghost_rows(f_rows)
   ! The lower ghost composition the last sweep of this run returned.
   real*8, intent(in) :: f_rows(1-Ng:0,n_species)
   previous_sweep_ghost_row      = f_rows
   previous_sweep_ghost_row_have = .true.
   end subroutine set_previous_sweep_ghost_rows

   !------------------------------------------!

   subroutine lower_ghost_seed_rows(f_entry, f_seed, source_used)
   ! THE COMPOSITION THE LOWER GHOST CELLS ARE HANDED, per unit mass, in the
   ! f_sp layout. With the key unset this returns f_entry unchanged, which is
   ! the state the caller already holds, so the call is a no-op. A named
   ! source that has nothing stored yet also returns f_entry, and says so.
   real*8, intent(in)  :: f_entry(1-Ng:0,n_species)
   real*8, intent(out) :: f_seed(1-Ng:0,n_species)
   character(len=*), intent(out) :: source_used
   character(len=256) :: src
   logical :: ok
   f_seed     = f_entry
   source_used = 'the entry state'
   src = ghost_composition_seed_source()
   if (len_trim(src) .eq. 0) return
   if (trim(src) .eq. 'the_state_file_ghost_rows') then
      if (state_file_ghost_row_have) then
         f_seed      = state_file_ghost_row
         source_used = 'the state file ghost rows'
      else
         source_used = 'the entry state (no state file ghost rows kept)'
      endif
   else if (trim(src) .eq. 'the_previous_sweep_ghost') then
      if (previous_sweep_ghost_row_have) then
         f_seed      = previous_sweep_ghost_row
         source_used = 'the previous sweep ghost'
      else
         source_used = 'the entry state (no previous sweep yet)'
      endif
   else
      if (.not. supplied_ghost_row_tried) then
         call read_supplied_ghost_rows(trim(src), ok)
         supplied_ghost_row_tried = .true.
         supplied_ghost_row_have  = ok
      endif
      if (supplied_ghost_row_have) then
         f_seed      = supplied_ghost_row
         source_used = 'the file '//trim(src)
      else
         source_used = 'the entry state (the file '//trim(src)//          &
                       ' was not read)'
      endif
   endif
   end subroutine lower_ghost_seed_rows

   !------------------------------------------!

   subroutine read_supplied_ghost_rows(path, ok)
   ! A file of lower ghost rows, in the format write_ghost_seed_rows writes:
   ! one line per ghost cell, the token 'ghost', the cell index, the species
   ! count and that many mass fractions. Rows for other indices are ignored.
   character(len=*), intent(in)  :: path
   logical,          intent(out) :: ok
   integer :: u, ios, j, ns, k, n_got
   character(len=16)    :: tok
   character(len=32768) :: line
   real*8 :: row(n_species)
   ok    = .false.
   n_got = 0
   open(newunit = u, file = path, status = 'old', action = 'read',        &
        iostat = ios)
   if (ios .ne. 0) then
      write(*,'(A)') ' (base_boundary) the ghost seed file '//trim(path)//&
           ' could not be opened; the entry state is used'
      return
   endif
   do
      read(u,'(A)',iostat=ios) line
      if (ios .ne. 0) exit
      if (len_trim(line) .eq. 0) cycle
      line = adjustl(line)
      if (line(1:1) .eq. '#') cycle
      read(line,*,iostat=ios) tok, j, ns, (row(k), k = 1, n_species)
      if (ios .ne. 0) cycle
      if (trim(tok) .ne. 'ghost') cycle
      if (ns .ne. n_species) cycle
      if (j .lt. 1-Ng .or. j .gt. 0) cycle
      supplied_ghost_row(j,:) = row
      n_got = n_got + 1
   enddo
   close(u)
   ok = (n_got .eq. Ng)
   if (.not. ok) write(*,'(A,I0,A,I0,A)') ' (base_boundary) the ghost'//  &
        ' seed file carried ', n_got, ' of ', Ng, ' rows; the entry'//    &
        ' state is used'
   end subroutine read_supplied_ghost_rows

   !------------------------------------------!

   subroutine write_ghost_seed_rows(f_rows)
   ! EXHALE_GHOST_SEED_WRITE=<path>, default off: the lower ghost rows as
   ! they were installed, in the format the seed key reads back. It exists so
   ! that a perturbed seed is made from a measured one rather than assembled
   ! by hand.
   real*8, intent(in) :: f_rows(1-Ng:0,n_species)
   character(len=256) :: path
   integer :: u, ios, j, k
   call get_environment_variable('EXHALE_GHOST_SEED_WRITE', path)
   if (len_trim(path) .eq. 0) return
   open(newunit = u, file = trim(path), status = 'replace',               &
        action = 'write', iostat = ios)
   if (ios .ne. 0) return
   write(u,'(A)') '# the lower ghost composition, mass fractions in the'//&
        ' f_sp layout'
   write(u,'(A)') '# ghost <cell> <n_species> <f(1)> ... <f(n_species)>'
   do j = 1-Ng, 0
      write(u,'(A,2I6)', advance = 'no') 'ghost ', j, n_species
      do k = 1, n_species
         write(u,'(ES26.16E3)', advance = 'no') f_rows(j,k)
      enddo
      write(u,'(A)') ''
   enddo
   close(u)
   write(*,'(A)') ' (base_boundary) the installed lower ghost'//          &
        ' composition was written to '//trim(path)
   end subroutine write_ghost_seed_rows

   !------------------------------------------!

   subroutine set_base_ghost_state_record(j, ntot, ne_cell, x_h2, x_h2p,  &
                                          x_h3p, x_hehp, x_hii, x_heii,  &
                                          x_heiii, charge_gap,           &
                                          element_gap, reaction_res,     &
                                          partition_res, passes)
   ! The ghost a sweep returned, for the ghost record. A measurement: no
   ! boundary expression reads any of it.
   integer, intent(in) :: j, passes
   real*8,  intent(in) :: ntot, ne_cell, x_h2, x_h2p, x_h3p, x_hehp
   real*8,  intent(in) :: x_hii, x_heii, x_heiii, charge_gap, element_gap
   real*8,  intent(in) :: reaction_res, partition_res
   if (j .lt. 1-Ng .or. j .gt. 0) return
   ghost_state_ntot(j)           = ntot
   ghost_state_ne(j)             = ne_cell
   ghost_state_x_h2(j)           = x_h2
   ghost_state_x_h2p(j)          = x_h2p
   ghost_state_x_h3p(j)          = x_h3p
   ghost_state_x_hehp(j)         = x_hehp
   ghost_state_x_hii(j)          = x_hii
   ghost_state_x_heii(j)         = x_heii
   ghost_state_x_heiii(j)        = x_heiii
   ghost_state_charge_gap(j)     = charge_gap
   ghost_state_element_gap(j)    = element_gap
   ghost_state_reaction_res(j)   = reaction_res
   ghost_state_partition_res(j)  = partition_res
   ghost_state_closure_passes(j) = passes
   ghost_state_recorded          = .true.
   end subroutine set_base_ghost_state_record

   !------------------------------------------!

   subroutine report_base_boundary_model(tag)
   ! WHICH LOWER BOUNDARY THIS RUN IS SOLVING, in one block: the model, the
   ! prescribed reservoir with its version, and, where a sweep has measured
   ! them, the ghost's own counts beside the prescribed ones and the closure
   ! of the ghost's molecular partition.
   !
   ! The reservoir's particles per unit mass and the ghost's own are printed
   ! together because the model does not reconcile them: the ghost pressure
   ! base_ghost_averages builds is base_reservoir_nhat rho T while the
   ! caloric map that turns that ghost into an energy and a sound speed uses
   ! the composition's own count. The distance between them is a statement
   ! about the model and is reported rather than split between two
   ! expressions.
   character(len=*), intent(in) :: tag
   write(*,'(A)') ' [base boundary model] '//trim(tag)
   write(*,'(A)') '   model  '//base_boundary_model_id
   write(*,'(A)') '   ghost composition seed: '//                        &
        base_ghost_composition_seed_id
   if (ghost_composition_seed_armed()) write(*,'(A)') '   NOTE: EXHALE'//&
        '_GHOST_COMPOSITION_SEED replaces that seed for this run, which'//&
        ' is a measurement and not a certificate'
   write(*,'(A,I0,A)') '   prescribed reservoir (version ',               &
        base_reservoir_prescription_version, '):'
   write(*,'(A,ES23.16,A,ES23.16)') '     p [p0] ', base_reservoir_p,     &
        '   T [T0] ', base_reservoir_T
   write(*,'(A,ES23.16,A,ES23.16)') '     particles per unit mass ',      &
        base_reservoir_nhat, '   level radius [Rp] ', r_base_level
   if (base_ghost_counts_measured) then
      write(*,'(A,ES23.16,A,ES23.16)') '   solved ghost counts [n0]:'//   &
           ' heavy ', base_ghost_particle_count, '   electrons ',         &
           base_ghost_electron_count
   endif
   if (base_ghost_closure_measured) then
      write(*,'(A,ES12.5,A,I0,A)') '   ghost H2 partition closed to ',    &
           base_ghost_closure_residual, ' in at most ',                   &
           base_ghost_closure_passes, ' passes'
   endif
   if (base_ghost_fixed_point_measured) then
      write(*,'(A,ES12.5,A,ES12.5)') '   ghost composition is a fixed'//  &
           ' point of its own solve to ', base_ghost_fixed_point_move,    &
           ' against ', base_ghost_fixed_point_tol
      write(*,'(A,ES12.5,A,ES12.5,A,I0,A)') '     its particle plus'//    &
           ' electron count to ', base_ghost_thermal_move, ' against ',   &
           base_ghost_count_tol, ', reached in ',                         &
           base_ghost_fixed_point_passes, ' applications'
   endif
   end subroutine report_base_boundary_model

   !------------------------------------------!

   subroutine continue_hydrostatic_isentrope(jcomp, nhat, r_a, rho_a, T_a,  &
                                             r_b, rho_b, T_b)
   ! Continue a state along the hydrostatic isentrope that passes through it:
   ! at rest and at constant entropy the Bernoulli function is uniform,
   !
   !     h + phi = const ,      h = nhat ( u(T) + T )    [per unit mass]
   !
   ! with u the internal energy per particle of the caloric EOS, so
   !
   !     u(T_b) + T_b = u(T_a) + T_a - (phi(r_b) - phi(r_a))/nhat .
   !
   ! Solved by Newton on T with dh/dT = c_v + 1 (the heat capacity at constant
   ! pressure per particle), which is exact for the mixture and closed form
   ! for a monatomic gas.  The density then follows from the isentrope itself,
   ! d ln rho = c_v(T) d ln T:
   !
   !     rho_b = rho_a (T_b/T_a)^c_v(T_mid) .
   !
   ! VALIDITY.  The midpoint c_v makes the density exact to third order in
   ! ln(T_b/T_a).  This routine is only ever used over half a cell and over
   ! the two ghost cells (cell 1 to the face, the face to the ghosts below
   ! it), i.e. |ln(T_b/T_a)| < (gamma-1)/gamma * Ng*dr/H,
   ! which on the hot-Uranus base grid is 6e-3 and on any grid a run of this
   ! code uses is below 0.1; the density error is then below 1e-7 relative.
   ! It is NOT valid as a general isentrope integrator over a scale height.
   integer, intent(in)  :: jcomp
   real*8,  intent(in)  :: nhat, r_a, rho_a, T_a, r_b
   real*8,  intent(out) :: rho_b, T_b
   real*8 :: h_target, f, cp, T, cv_mid
   integer :: it

   h_target = internal_energy_per_particle(jcomp, T_a) + T_a               &
              - (phi(r_b) - phi(r_a))/nhat

   ! Newton from T_a.  The Bernoulli function is strictly increasing in T
   ! (c_v > 0), so the root is unique; three iterations reach round-off over
   ! the increments this routine is used for.
   T = T_a
   do it = 1, 20
      cp = heat_capacity_per_particle(jcomp, T) + 1.0d0
      f  = internal_energy_per_particle(jcomp, T) + T - h_target
      if (abs(f) .le. 1.0d-15*abs(h_target)) exit
      T  = T - f/cp
      ! A potential drop deep enough to extrapolate through zero temperature
      ! is outside the range stated above; fall back on the isothermal
      ! continuation, which is the c_v -> infinity limit of the same relation
      ! and keeps the state admissible.
      if (.not. (T .gt. 0.0d0)) then
         T = T_a
         exit
      endif
   enddo
   T_b = T

   cv_mid = heat_capacity_per_particle(jcomp, 0.5d0*(T_a + T_b))
   if (T_b .eq. T_a) then
      ! Isothermal fallback of the guard above: hydrostatic at fixed T.
      rho_b = rho_a*exp(-(phi(r_b) - phi(r_a))/(nhat*T_a))
   else
      rho_b = rho_a*exp(cv_mid*log(T_b/T_a))
   endif

   end subroutine continue_hydrostatic_isentrope

   !------------------------------------------!

   real*8 function isentropic_density_at_pressure(jcomp, rho_a, p_a, T_a,  &
                                                  p_b)
   ! Density on the isentrope through (rho_a, p_a) at pressure p_b:
   ! d ln rho = d ln p / gamma_eff, evaluated with gamma_eff frozen at T_a.
   !
   ! VALIDITY.  Freezing gamma_eff costs a relative density error of order
   ! (d ln gamma/d ln T) (ln p_b/p_a)^2.  Both pressures here are pressures AT
   ! THE SAME FACE -- the reservoir's and the interior's -- so their ratio
   ! departs from unity only by the boundary's own error, of order dr/H
   ! (4e-3 on the hot-Uranus base grid), and the frozen-index error is below
   ! 1e-5 of that.
   integer, intent(in) :: jcomp
   real*8,  intent(in) :: rho_a, p_a, T_a, p_b
   real*8 :: gam
   gam = adiabatic_index_at_T(jcomp, T_a)
   isentropic_density_at_pressure = rho_a*exp(log(p_b/p_a)/gam)
   end function isentropic_density_at_pressure

   !------------------------------------------!

   real*8 function characteristic_branch_weight(x)
   ! C1 handover from 0 at x >= 1 to 1 at x <= -1, the cubic smoothstep.
   ! x is the face Mach number scaled by the blend width, so this is the
   ! weight with which a characteristic that is LEAVING the domain takes over
   ! the datum a characteristic that was ENTERING used to carry.
   ! Its one use is the SUPERSONIC-OUTFLOW transition at M_i = -1, where the
   ! count of outgoing characteristics changes from two to three and the
   ! whole face state passes from the reservoir to the interior.  That
   ! transition is a change of the characteristic count and its width is a
   ! declared numerical regularization (base_face_mach_blend).  The entropy
   ! source of the contact is NOT handed over by this function: it is
   ! upwinded on the direction of the level's mass flux, which is a discrete
   ! fact and takes no intermediate value.
   real*8, intent(in) :: x
   if (x .ge.  1.0d0) then
      characteristic_branch_weight = 0.0d0
   else if (x .le. -1.0d0) then
      characteristic_branch_weight = 1.0d0
   else
      characteristic_branch_weight = 0.5d0 - 0.75d0*x + 0.25d0*x*x*x
   endif
   end function characteristic_branch_weight

   !------------------------------------------!

   subroutine read_base_branch_options()
   ! The three scales of the condition, overridable for a control
   ! experiment and read once.  Nothing else reads these variables and a
   ! production run sets none of them.
   character(len=32) :: env
   integer :: ios
   real*8  :: v
   if (base_branch_options_read) return
   base_branch_options_read = .true.
   env = ' '
   call get_environment_variable('EXHALE_BASE_MACH_BLEND', env)
   if (len_trim(env) .gt. 0) then
      read(env,*,iostat=ios) v
      if (ios .eq. 0 .and. v .gt. 0.0d0) base_face_mach_blend = v
   endif
   env = ' '
   call get_environment_variable('EXHALE_BASE_WIND_SPREAD', env)
   if (len_trim(env) .gt. 0) then
      read(env,*,iostat=ios) v
      if (ios .eq. 0 .and. v .gt. 0.0d0) base_wind_window_spread = v
   endif
   env = ' '
   call get_environment_variable('EXHALE_BASE_WIND_AMPLITUDE', env)
   if (len_trim(env) .gt. 0) then
      read(env,*,iostat=ios) v
      if (ios .eq. 0 .and. v .gt. 0.0d0) base_wind_window_amplitude = v
   endif
   end subroutine read_base_branch_options

   !------------------------------------------!

   subroutine wind_window_mass_flux(W, F_wind, d_window, have_F)
   ! The mass flux the state carries AND HOW COHERENT IT IS, read where the
   ! cell-centred product IS that flux: over the wind window r >= r_flux, the
   ! window the convergence gate is defined on (flux_spread_of_state), the
   ! mean of rho v r^2 and its relative standard deviation
   !
   !     d_window = sqrt( <F^2> - <F>^2 ) / |<F>| ,   j = j_flux .. N .
   !
   ! The base layer is deliberately excluded: that is where the collocated
   ! mode lives and where the product is not a flux (see the note at
   ! base_branch_on_wind_flux).
   !
   ! THE SIGN IS KEPT: a window whose flux reverses inside itself is not one
   ! wind, and the second moment about the SIGNED mean reads that reversal as
   ! a spread of order one, where the moments of the magnitudes would not.
   !
   ! WHY THE SECOND MOMENT AND NOT THE RANGE: the range is not differentiable
   ! where the cell carrying an extremum changes, and this quantity is an
   ! argument of a residual.  The design note above carries the measured slope
   ! jump and the recalibrated threshold.
   !
   ! The variance is accumulated about the mean in a second pass rather than
   ! as <F^2> - <F>^2 in one: on a coherent window the two terms agree to
   ! nine digits and the one-pass form loses the difference to cancellation,
   ! which is the whole of d_window on exactly the states the weight has to
   ! resolve.
   !
   ! have_F is false where the window is empty or the mean is not a usable
   ! number.  An exactly zero mean is NOT such a case: it is the limit the
   ! amplitude weight of the caller handles, and F_wind = 0 carries it.
   real*8, intent(in)  :: W(3,1-Ng:N+Ng)
   real*8, intent(out) :: F_wind, d_window
   logical, intent(out):: have_F
   integer :: j, n_in
   real*8  :: f_j, var
   ! The value a window that says nothing is given: above twice the
   ! threshold, so the shape weight is zero, and FINITE, so that forming the
   ! weight's argument raises no overflow.
   real*8, parameter :: d_window_silent = 3.0d0
   F_wind = 0.0d0;  n_in = 0
   d_window = d_window_silent*base_wind_window_spread
   do j = max(j_flux,1), N
      f_j    = W(1,j)*W(2,j)*r(j)*r(j)
      F_wind = F_wind + f_j
      n_in   = n_in + 1
   enddo
   have_F = (n_in .gt. 1)
   if (have_F) then
      F_wind = F_wind/dble(n_in)
      have_F = (F_wind .eq. F_wind) .and. (abs(F_wind) .lt. huge(1.0d0))
   endif
   if (have_F) then
      var = 0.0d0
      do j = max(j_flux,1), N
         f_j = W(1,j)*W(2,j)*r(j)*r(j) - F_wind
         var = var + f_j*f_j
      enddo
      var = var/dble(n_in)
      if (F_wind .ne. 0.0d0) then
         ! Capped far above every measured value (the largest is 2.40, the
         ! hydrostatic_column fixture) so that the number stays a diagnostic
         ! and the cap only keeps the weight's argument finite where the
         ! mean is at the rounding floor of a window carrying no flux.
         d_window = min(sqrt(max(var,0.0d0))/abs(F_wind), 1.0d3)
         if (.not. (d_window .eq. d_window)) have_F = .false.
      else
         ! The shape of a window of zero mean says nothing the amplitude
         ! weight does not already silence; leave it at the value that closes
         ! the gate.
         d_window = d_window_silent*base_wind_window_spread
      endif
   endif
   end subroutine wind_window_mass_flux

   !------------------------------------------!

   subroutine characteristic_base_face_state(W1, nhat1, T1, F_wind,       &
                                             d_window, have_F, Wface, Wi)
   ! The primitive state at the base face r_edg(0), and the interior state
   ! continued to the same face (returned so the ghost construction and the
   ! diagnostics read the same numbers the solve used).
   !
   ! W1, nhat1, T1 are the CELL AVERAGE of the first interior cell, its
   ! particle count per unit mass and its temperature.  The interior state at
   ! the face is built from cell 1 alone, along cell 1's own hydrostatic
   ! isentrope, and NOT from the reconstruction stencil: the stencil reaches
   ! into the ghost, and a boundary condition that reads its own output is the
   ! circular dependency this module exists to break.  Continuing along the
   ! hydrostatic isentrope rather than extrapolating linearly is what makes
   ! the condition WELL BALANCED: in a hydrostatic atmosphere at rest the
   ! interior face pressure equals the reservoir's and (C-) returns v_b = 0
   ! exactly, at any resolution.
   real*8, intent(in)  :: W1(3), nhat1, T1
   ! The mass flux the state carries in the wind window, the relative
   ! standard deviation of that flux inside the window, and whether they are
   ! usable numbers: the quantities the direction of the contact is read
   ! from where the window has standing.
   real*8, intent(in)  :: F_wind, d_window
   logical, intent(in) :: have_F
   real*8, intent(out) :: Wface(3), Wi(3)
   real*8 :: rho_i, T_i, p_i, v_i, c_i, M_i, M_wind, s_wind
   real*8 :: rho_res, T_res, p_res
   real*8 :: rho_b, v_b, p_b, c_b, rho_rev
   real*8 :: w_rev, w_out, rb, r1, v_floor
   logical :: reverse_flow

   rb = r_edg(0)
   r1 = r(1)

   ! ---- the interior state at the face ----
   call continue_hydrostatic_isentrope(1, nhat1, r1, W1(1), T1, rb,       &
                                       rho_i, T_i)
   p_i = nhat1*rho_i*T_i
   ! Velocity along the same continuation: a steady flow carries rho v r^2,
   ! so this is the interior's velocity at the face and not its value half a
   ! cell above it.
   v_i = W1(1)*W1(2)*r1*r1/(rho_i*rb*rb)
   c_i = sqrt(adiabatic_index_at_T(1, T_i)*p_i/rho_i)
   M_i = v_i/c_i

   Wi(1) = rho_i;  Wi(2) = v_i;  Wi(3) = p_i

   ! ---- the reservoir state at the face ----
   ! The face is the level (base_boundary_states checks it), so this is the
   ! stated reservoir itself, with no transport along any isentrope.
   call reservoir_state_at_level(rho_res, T_res, p_res)

   ! ---- the two reservoir conditions, and the one interior relation ----
   !
   ! Pressure is a reservoir condition in every subsonic branch.  Entropy is
   ! one only while gas enters: on reversal the entropy at the face is the
   ! interior's, advected out, and the reservoir has nothing to say about it.
   ! The upwind choice below is that decision.  (The Riemann solver downstream
   ! upwinds the contact wave as well, so this choice acts on the acoustic
   ! part of the flux; imposing it here is what makes the STATE, and hence
   ! the ghost cells that every other module reads, carry the right entropy.)
   ! (C-) alone leaves one relation between p_b and v_b; the second reservoir
   ! condition closes it.  Written as one expression so the two choices differ
   ! by a weight and not by a branch: at weight 0 this is p_b = p_res, and at
   ! weight 1 it is the pair of (C-) with p_b + rho c v_b = p_res solved
   ! together.
   p_b = p_res + 0.5d0*base_incoming_invariant_weight                     &
                 *(p_i - rho_i*c_i*v_i - p_res)

   ! ---- the acoustic matching, before the contact is upwinded ----
   !
   ! (C-): the outgoing acoustic invariant of the first interior state,
   ! closed with the reservoir's pressure.  The face velocity is a RESULT of
   ! the matching, so it is available before the entropy source is chosen
   ! and the choice is not circular: v_b does not depend on rho_b.
   v_b = v_i + (p_b - p_i)/(rho_i*c_i)

   ! ---- the direction of the contact, and the trace it selects ----
   !
   ! The upwind side of the contact is the side the level's mass flux comes
   ! from.  Its sign is read from the wind window where the window has
   ! standing, and from the matched face velocity where it has none; the
   ! comment at base_wind_window_spread carries the two conditions and the
   ! measured populations they separate, and the block "WHICH QUANTITY THE
   ! DIRECTION OF THE CONTACT IS READ FROM" says why the interior CELL
   ! velocity is not one of the candidates.
   !
   ! THE CHOICE IS DISCRETE AND HAS NO WIDTH.  A contact with different
   ! entropies on its two sides has direction-dependent traces; a face
   ! density between them is an interpolation and not a law of thermal
   ! contact, and the average of the two isentropes that a symmetric
   ! handover returns at zero has no physics behind it.  Where the window
   ! has standing the states measured above sit decades from its two
   ! thresholds, but not every state does: the atomic x0.01, He/H = 9.7
   ! state g0002 of LHS 1140 b sits one decade away (d = 2.04e-4 against
   ! 2e-3), and a finite-difference probe of about two standard arcs
   ! crosses the switch (MEASURED 2026-09-23), so the residual a Newton
   ! differentiates can carry this step.  Where the window does not have
   ! standing the matched face velocity of a state that
   ! moves is of order 1 to 100 cm/s away from zero.  The one state that IS
   ! at the switch is the one it was built for, a column in exact
   ! hydrostatic balance with the level, and the floor below decides it for
   ! the reservoir.
   !
   ! AT REST AND AT INFLOW THE RESERVOIR OWNS THE LEVEL.  That is the
   ! assumption stated at the head of the file: the lower atmosphere is a
   ! heat bath, so a column at rest above it carries the level's own
   ! temperature.  Only a REVERSE flow hands the trace to the interior, and
   ! then the reservoir keeps exactly the one condition its single entering
   ! characteristic allows, the pressure p_b above.
   M_wind = 0.0d0
   s_wind = 0.0d0
   if (have_F) then
      M_wind = F_wind/(rho_i*rb*rb)/c_i
      if (abs(M_wind) .ge. base_wind_window_amplitude .and.               &
          d_window .le. base_wind_window_spread) s_wind = 1.0d0
   endif
   if (s_wind .gt. 0.0d0) then
      reverse_flow = (M_wind .lt. 0.0d0)
   else
      ! ZERO TO THE MATCHING'S OWN PRECISION IS ZERO, AND THE RESERVOIR
      ! OWNS IT.  v_b is formed as v_i + (p_b - p_i)/(rho_i c_i), and at a
      ! pressure-balanced contact that difference is a cancellation: what is
      ! left is the rounding of the two pressures, so the SIGN of v_b there
      ! is the sign of an arithmetic remainder and not of a flow.  MEASURED
      ! on a discrete hydrostatic column built on the reservoir's own
      ! isentrope, v_b comes out at -1.2e-16 of the code velocity unit,
      ! which a bare sign test would read as a reverse flow and hand the
      ! level to the interior -- the one state choice A exists to answer.
      ! The floor below is that remainder and not a tuned width: it is a few
      ! rounding units of the terms the subtraction is made of, it scales
      ! with the state, and it vanishes with the arithmetic precision.
      v_floor = 4.0d0*epsilon(1.0d0)                                      &
                *((abs(p_i) + abs(p_b))/(rho_i*c_i) + abs(v_i))
      reverse_flow = (v_b .lt. -v_floor)
   endif
   w_rev = 0.0d0
   if (reverse_flow) w_rev = 1.0d0
   if (w_rev .gt. 0.0d0) n_base_reversal_evals = n_base_reversal_evals + 1
   rho_rev = isentropic_density_at_pressure(1, rho_i, p_i, T_i, p_b)
   rho_b   = (1.0d0 - w_rev)*rho_res + w_rev*rho_rev

   ! ---- supersonic branches ----
   !
   ! Supersonic OUTFLOW takes every characteristic out of the domain, so
   ! nothing may be prescribed and the face state is the interior's.
   w_out = characteristic_branch_weight((M_i + 1.0d0)/base_face_mach_blend)
   if (w_out .gt. 0.0d0) then
      n_base_supersonic_outflow = n_base_supersonic_outflow + 1
      rho_b = (1.0d0 - w_out)*rho_b + w_out*rho_i
      v_b   = (1.0d0 - w_out)*v_b   + w_out*v_i
      p_b   = (1.0d0 - w_out)*p_b   + w_out*p_i
   endif
   !
   ! Supersonic INFLOW would take all three characteristics into the domain
   ! and the reservoir would owe three conditions.  A (p, s) reservoir states
   ! two; this model has no third, because it has no lower-atmosphere
   ! velocity.  The state below is therefore NOT a solution of a well-posed
   ! problem, and the run says so through check_base_inflow_is_subsonic in
   ! BC_Apply rather than inventing a velocity here.

   ! Sound speed of the face state, at the BASE composition.  On the
   ! reversal branch the density is the interior's and cell 1's composition
   ! would be the consistent one, so this reads a sound speed of the
   ! reservoir's mixture at the interior's density there.  It is a stated
   ! approximation and its only consumers are the Mach cap below, which is a
   ! guard against a linear relation asking for a hypersonic face far from a
   ! solution, and the diagnostic; no flux and no residual reads it.
   c_b = sqrt(adiabatic_index_at_T(0, p_b/(base_reservoir_nhat*rho_b))    &
              *p_b/rho_b)
   if (abs(v_b) .gt. base_face_mach_max*c_b) then
      v_b = sign(base_face_mach_max*c_b, v_b)
      n_base_face_mach_limited = n_base_face_mach_limited + 1
   endif

   base_face_mach_last  = v_b/c_b
   base_face_blend_last = w_rev
   base_face_Mi_last      = M_i
   base_face_Mwind_last     = M_wind
   base_face_swind_last     = s_wind
   base_face_d_window_last  = d_window
   base_face_vi_last      = v_i
   base_face_vb_last      = v_b
   base_face_rho_res_last = rho_res
   base_face_rho_rev_last = rho_rev
   base_face_T_res_last   = T_res
   base_face_T_i_last     = T_i

   Wface(1) = rho_b;  Wface(2) = v_b;  Wface(3) = p_b

   end subroutine characteristic_base_face_state

   !------------------------------------------!

   subroutine reservoir_state_at_level(rho_res, T_res, p_res)
   ! The lower atmosphere's own state at the base level, which is the face:
   ! the prescribed (p, T) and the density p = n_hat rho T gives them.
   real*8, intent(out) :: rho_res, T_res, p_res
   T_res   = base_reservoir_T
   p_res   = base_reservoir_p
   rho_res = base_reservoir_p/(base_reservoir_nhat*base_reservoir_T)
   end subroutine reservoir_state_at_level

   !------------------------------------------!

   subroutine reservoir_face_state(Wface)
   ! The reservoir at the face, at rest. This is the whole boundary when the
   ! interior has no admissible state to supply a characteristic from, and
   ! it is the state every other branch starts from.
   real*8, intent(out) :: Wface(3)
   real*8 :: rho_res, T_res, p_res
   call reservoir_state_at_level(rho_res, T_res, p_res)
   Wface(1) = rho_res
   Wface(2) = 0.0d0
   Wface(3) = p_res
   end subroutine reservoir_face_state

   !------------------------------------------!

   subroutine require_base_face_at_level()
   ! The face the condition is imposed at must be the level the reservoir
   ! is stated at: every grid define_grid builds has r_edg(0) = 1 exactly,
   ! and set_base_reservoir is called with r_level = 1. A grid built some
   ! other way (a test that writes r_edg itself) states its reservoir at its
   ! own face. Any other pairing would put a slab between the reservoir and
   ! the domain, which is the error the face-at-level model removes, so it
   ! is refused and not bridged.
   if (r_edg(0) .ne. r_base_level) then
      write(*,'(A)') ' (base_boundary) ERROR: the lower face of the grid'// &
           ' is not the level the reservoir is stated at.'
      write(*,'(A,ES23.16,A,ES23.16)') '   r_edg(0) = ', r_edg(0),         &
           '   r_base_level = ', r_base_level
      error stop 1
   endif
   end subroutine require_base_face_at_level

   !------------------------------------------!

   real*8 function base_reservoir_temperature_at(rq) result(T_res)
   ! The temperature of the lower atmosphere at radius rq, on the
   ! reservoir's own hydrostatic isentrope through (p_base, T_base) at the
   ! base composition.
   !
   ! THIS IS THE TEMPERATURE OF THE LEVEL AND NOT OF THE ADVECTED TRACE.  It
   ! is what an operator that exchanges energy with the lower atmosphere
   ! reads at the base face -- thermal conduction does, viscous_conduction.f90
   ! -- so that the bath's temperature there is a statement of the boundary
   ! and does not follow the direction the contact happens to be upwinded
   ! on.  At rest and at inflow it is the temperature the advective ghost
   ! carries anyway, because the ghost is then the reservoir's own
   ! continuation; during a reverse flow the two differ and this is the one
   ! the bath states.
   !
   ! VALIDITY is continue_hydrostatic_isentrope's: the midpoint heat
   ! capacity makes the density exact to third order in ln(T/T_base), which
   ! holds over the two ghost cells below the level this is asked for and
   ! not over a scale height.
   real*8, intent(in) :: rq
   real*8 :: rho_q
   call continue_hydrostatic_isentrope(0, base_reservoir_nhat,            &
                                       r_base_level,                      &
                                       base_reservoir_p                   &
                                       /(base_reservoir_nhat              &
                                         *base_reservoir_T),              &
                                       base_reservoir_T, rq,              &
                                       rho_q, T_res)
   end function base_reservoir_temperature_at

   !------------------------------------------!

   subroutine base_ghost_averages(Wface, Wghost, Wface_lower)
   ! Ghost CELL AVERAGES that are consistent with the face state, and the
   ! point state at the next face down.
   !
   ! The ghost is not a copy of the face: it is the volume average of the
   ! hydrostatic isentrope that passes through the face state, continued
   ! BELOW it, with the velocity carrying a constant rho v r^2 so the
   ! reconstruction sees no spurious velocity gradient at the boundary.
   ! Writing the face value into the ghost instead is the half-cell error
   ! this module removes; the two differ by dr/(2H), which is first order in
   ! dr in the state and ZEROTH order in the momentum residual it produces.
   !
   ! The quadrature is 8-point Gauss-Legendre, volume weighted, which is far
   ! beyond the reconstruction's order and costs two cells per boundary call.
   real*8, intent(in)  :: Wface(3)
   real*8, intent(out) :: Wghost(3,1-Ng:0)
   real*8, intent(out) :: Wface_lower(3)
   real*8 :: T_face, flux_r2, a_, b_, xx, num_r, num_p, den, rho_q, T_q
   real*8 :: rho_l, T_l
   integer :: j, k
   real*8, parameter :: xg(8) = (/                                        &
      -0.9602898564975363d0, -0.7966664774136267d0,                       &
      -0.5255324099163290d0, -0.1834346424956498d0,                       &
       0.1834346424956498d0,  0.5255324099163290d0,                       &
       0.7966664774136267d0,  0.9602898564975363d0 /)
   real*8, parameter :: wg(8) = (/                                        &
       0.1012285362903763d0,  0.2223810344533745d0,                       &
       0.3137066458778873d0,  0.3626837833783620d0,                       &
       0.3626837833783620d0,  0.3137066458778873d0,                       &
       0.2223810344533745d0,  0.1012285362903763d0 /)

   T_face  = Wface(3)/(base_reservoir_nhat*Wface(1))
   flux_r2 = Wface(1)*Wface(2)*r_edg(0)*r_edg(0)

   do j = 1-Ng, 0
      if (j .eq. 1-Ng) then
         a_ = r(j) - 0.5d0*dr_j(j);  b_ = r_edg(j)
      else
         a_ = r_edg(j-1);  b_ = r_edg(j)
      endif
      num_r = 0.0d0;  num_p = 0.0d0;  den = 0.0d0
      do k = 1, 8
         xx = 0.5d0*(a_ + b_) + 0.5d0*(b_ - a_)*xg(k)
         call continue_hydrostatic_isentrope(0, base_reservoir_nhat,      &
                                             r_edg(0), Wface(1), T_face,  &
                                             xx, rho_q, T_q)
         num_r = num_r + wg(k)*rho_q*xx*xx
         num_p = num_p + wg(k)*base_reservoir_nhat*rho_q*T_q*xx*xx
         den   = den   + wg(k)*xx*xx
      enddo
      Wghost(1,j) = num_r/den
      Wghost(3,j) = num_p/den
      ! Constant rho v r^2 through the ghost, evaluated at the cell centre
      ! with the cell's own average density.
      Wghost(2,j) = flux_r2/(Wghost(1,j)*r(j)*r(j))
   enddo

   call continue_hydrostatic_isentrope(0, base_reservoir_nhat, r_edg(0),  &
                                       Wface(1), T_face, r_edg(1-Ng),     &
                                       rho_l, T_l)
   Wface_lower(1) = rho_l
   Wface_lower(3) = base_reservoir_nhat*rho_l*T_l
   Wface_lower(2) = flux_r2/(rho_l*r_edg(1-Ng)*r_edg(1-Ng))

   end subroutine base_ghost_averages

   !------------------------------------------!

   subroutine base_boundary_states(W, Wface, Wghost, Wface_lower)
   ! The single entry point.  Everything the lower boundary states is derived
   ! here, from the first interior CELL AVERAGE, so that the two places that
   ! need it -- the ghost cell averages of Apply_BC_W and the face state of
   ! Rec_BC -- are answers from one function and cannot drift apart.  That
   ! drift is what the previous closure had: the valve read the cell-1 average
   ! in one path and a reconstructed face velocity a full cell higher in the
   ! other.
   real*8, intent(in)  :: W(3,1-Ng:N+Ng)
   real*8, intent(out) :: Wface(3), Wghost(3,1-Ng:0), Wface_lower(3)
   real*8 :: nhat1, T1, Wi(3), F_wind_layer, d_window_layer
   logical :: have_F_layer

   call require_base_face_at_level()

   ! An inadmissible first cell cannot state an outgoing characteristic, and
   ! dividing by it would make the boundary itself the source of the NaN. The
   ! reservoir alone is then the whole boundary -- which is what a boundary
   ! condition is when the interior has nothing to say -- and the run's own
   ! positivity detector reports the cell.
   if (.not. (W(1,1) .gt. 0.0d0 .and. W(3,1) .gt. 0.0d0 .and.             &
              n_part_cell1 .gt. 0.0d0)) then
      call reservoir_face_state(Wface)
      call base_ghost_averages(Wface, Wghost, Wface_lower)
      n_base_interior_inadmissible = n_base_interior_inadmissible + 1
      return
   endif

   ! Particle count per unit mass of cell 1, and its temperature.  n_part_cell1
   ! is n_tot(1) + n_e(1) in units of n0, refreshed by the composition solve
   ! (get_species_densities, the single policy point), and p = n_part T in
   ! these units, so both follow without duplicating any bookkeeping.
   nhat1 = n_part_cell1/W(1,1)
   T1    = W(3,1)/n_part_cell1

   call read_base_branch_options()
   call wind_window_mass_flux(W, F_wind_layer, d_window_layer,            &
                              have_F_layer)
   call characteristic_base_face_state(W(:,1), nhat1, T1, F_wind_layer,   &
                                       d_window_layer, have_F_layer,      &
                                       Wface, Wi)
   call base_ghost_averages(Wface, Wghost, Wface_lower)

   end subroutine base_boundary_states

   ! End of module
   end module base_boundary

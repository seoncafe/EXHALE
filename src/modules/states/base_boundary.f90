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
   ! ladder), which is the base sawtooth of docs/p44_base_sawtooth.md section 7
   ! reproduced by Euler plus gravity plus that boundary alone, with no
   ! chemistry, no radiation and no wind.  Separately, copying the interior
   ! velocity into the ghost (the one-way valve) made the boundary reflect
   ! 95 percent of an acoustic pulse leaving through it; holding the same
   ! numbers with a stated velocity instead left 3 percent.
   ! Design note and every measurement: docs/phaseC_characteristic_base_bc.md.
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
   ! a spurious base inflow whenever cell 1 was hot
   ! (docs/heitr_metals_and_convergence_notes.md section 2.5); that is what a
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
   ! THE BASE LEVEL IS r_base_level, NOT THE FACE.  "Base BC: pressure" states
   ! the pressure at r = 1 (the planet radius the handoff quotes: 1 microbar in
   ! Koskinen et al. 2022 section 3.1), and the reservoir state is carried from
   ! there to r_edg(0) along its own hydrostatic isentrope.  So the level the
   ! user states does not move when the grid changes, and the half-cell error
   ! this module exists to remove is not reintroduced at the level itself.

   use global_parameters
   use grav_func, only: phi
   use caloric_eos, only: internal_energy_per_particle,                  &
                          heat_capacity_per_particle,                    &
                          adiabatic_index_at_T

   implicit none

   ! ---- the reservoir, set once by set_base_reservoir ----
   real*8 :: base_reservoir_p = 1.0d0   ! p at r_base_level          [p0]
   real*8 :: base_reservoir_T = 1.0d0   ! T at r_base_level          [T0]
   real*8 :: base_reservoir_nhat = 1.0d0! particles per unit mass at the base
   real*8 :: r_base_level = 1.0d0       ! radius the reservoir is stated at

   ! ---- the one designed-in constant ----
   !
   ! Width, in face Mach number, of the window over which the ENTROPY source
   ! of the face state is handed over from the reservoir (inflow) to the
   ! interior (reversal).  The branch is a genuine change in the number of
   ! conditions, so it is a kink in the residual, and the steady Newton solver
   ! needs a differentiable one.
   !
   ! The window must be far BELOW the physical operating point, or it would
   ! blend the reservoir entropy away where the boundary is supposed to state
   ! it.  On the hot-Uranus gate the base carries the wind's own mass flux,
   ! rho_b v_b r^2 = F_wind, which at the base density is v_b ~ 1.6 cm/s
   ! against c ~ 2.6e5 cm/s: a face Mach number of 6e-6.  A hot Jupiter's is
   ! larger.
   !
   ! 1e-6 WAS NOT FAR BELOW EVERY OPERATING POINT, and the planet that showed
   ! it is LHS 1140 b (item L21, 2026-09-15,
   ! docs/lhs1140b_stationary_L21_20260915.md).  On the certified atomic state
   ! of `atomic_scalar_gj1132_kzz1e9/HeH2.13` the base carries the wind's own
   ! flux at v ~ 0.064 cm/s against c = 1.38e5 cm/s, a face Mach number of
   ! 4.6e-7: BELOW the window, not far above it.  A boundary whose operating
   ! point sits inside its own handover blends away a part of the datum it
   ! exists to state, and on that state the two branches are not close --
   ! rho_rev/rho_res = 0.531 and T_i/T_res = 1.88 -- so the smoothstep was
   ! smearing a factor of two and not the difference between two isentropes
   ! that coincide.
   !
   ! 1e-8 stands 46 times below the LHS 1140 b operating point and 600 times
   ! below the hot-Uranus one, which is what the paragraph above asks for.
   ! The handover is C1 either way; a narrower window makes the derivative
   ! inside it larger and no converged state sits inside it, which is the
   ! whole design requirement.
   !
   ! THIS WIDTH CARRIES ONE STATEMENT, and did not always.  It is the width
   ! over which the reversal handover of a FACE MACH NUMBER is smoothed, and
   ! nothing else.  It was for a time also read as the scale below which the
   ! WIND WINDOW was deemed to say nothing, and that cost was measured: at the
   ! cold start of the wasp_full regression case the window carried
   ! M_wind = +3.70e-09 -- no wind at all -- while the first cell was at
   ! M_i = +1.15e-02, so |M_wind| was 0.37 of this width, the window was given
   ! a 9 per cent say, and the state moved 1.8e-04 in ONE step and 3.954e-04
   ! by step 300.  The window's standing is now read off the window's own flux
   ! spread instead (base_wind_window_spread below), which has no scale to
   ! borrow, and with it that case is bit for bit what it was before item L21
   ! (docs/lhs1140b_stationary_L21_20260915.md sections 11.4 and 12).
   ! Its influence is measured, not assumed: see the design note's section 7
   ! and the item above.  EXHALE_BASE_MACH_BLEND=<value> overrides it for a
   ! control experiment; nothing else reads that variable.
   real*8, save :: base_face_mach_blend = 1.0d-8

   ! WHICH QUANTITY THE HANDOVER IS KEYED ON.
   !
   ! The branch asks whether gas ENTERS the domain through the base face, and
   ! the honest answer is the sign of the mass flux the face carries.  This
   ! routine used the cell-centred product of the first cell, rho_1 v_1 r_1^2,
   ! continued to the face.  AT A BASE THAT CARRIES A COLLOCATED ODD-EVEN
   ! VELOCITY MODE THAT PRODUCT IS NOT A FLUX AND ITS SIGN IS NOT THE FLOW'S:
   ! measured on the certified LHS 1140 b atomic state, the Riemann face flux
   ! is +1.00000 F_wind at EVERY face of the grid including this one, while
   ! the cell-centred product reads -2.00 F_wind at cell 1 and -2.27 at cell 2
   ! (item L21; the hot-Uranus base of docs/p44_base_sawtooth.md section 3
   ! reads -196 F_wind against a face flux of +0.98).  Keyed on it, this
   ! boundary ran an INFLOW face 97 per cent in the reversal branch and kept
   ! the reservoir's entropy out of the state it exists to state.
   !
   ! flux_spread_of_state made the same correction for the convergence gate
   ! and says why there ("WHY THE FACE AND NOT rho v r^2 AT THE CELL CENTRE").
   ! The face flux itself is not available here -- it is produced by the
   ! Riemann solve this boundary feeds, and reading it would be the circular
   ! dependency this module exists to break, or would put history into a
   ! residual that has to stay a function of its argument.  What IS available
   ! and carries no collocated mode is the same conserved flux read WHERE the
   ! cell-centred product is the flux: the wind window r >= r_flux that the
   ! flux gate is already defined on.  In a state whose mass flux is
   ! conserved that is the flux through this face as well.
   !
   ! LIMITATION, stated rather than hidden: during a transient the base can
   ! carry a flux the wind does not yet, and then this reads the wind's sign
   ! and not the base's.  The branch it selects is still an admissible
   ! boundary condition, and no converged state has the two disagreeing.
   ! EXHALE_BASE_BRANCH_ON_CELL1=1 restores the old discriminant.
   logical, save :: base_branch_on_wind_flux = .true.

   ! WHEN THE WIND WINDOW HAS STANDING TO SPEAK, AND ON WHAT SCALE.
   !
   ! The question the discriminant asks of the window is not "is its flux
   ! large" but "is its flux ONE flux" -- a window that carries a wind carries
   ! the same rho v r^2 at every altitude in it.  That is a property of the
   ! shape of the flux profile and it has no scale of its own, so the weight
   ! of the window in the branch is read off the window's OWN relative spread
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
   ! swept through the argmax (item L26, section 8 of
   ! docs/lhs1140b_stationary_L26_20260916.md): the two one-sided derivatives
   ! of the face density converge at three step sizes to -242.2 and +0.0128,
   ! a slope jump of 1.9e4 at a point where the value itself is continuous.
   ! A residual that carries that kink is not differentiable, and the Newton
   ! solve differences through it.  The second moment has no argmax: it is a
   ! polynomial in the state divided by the mean, so it is smooth wherever
   ! the mean is nonzero, and its own square root is reached only at
   ! d_window = 0, where the weight below is saturated and flat.
   !
   ! THE FUNCTIONAL IS NO LONGER THE CONVERGENCE GATE'S.  flux_spread_of_state
   ! keeps the range, because a gate is evaluated once and is allowed a kink;
   ! this boundary sits inside a residual and is not.  The threshold below is
   ! therefore recalibrated on this functional and not carried over.
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
   ! above the loosest state that is a wind, and the tightest state that is
   ! not a wind sits 3.5 times above TWICE it, where the window falls silent.
   ! Both populations are in the saturated parts of the weight and no state
   ! sits in the handover.  The window
   ! speaks in full at or below the threshold and is silent at or above twice
   ! it, on the cubic smoothstep the branch itself uses, so the handover is
   ! C1 in d_window as well.
   !
   ! WHY THE FULL WEIGHT IS REACHED AT d_window = threshold AND NOT AT ZERO.
   ! The cell-centred product is a reconstruction of the state, so its spread
   ! over the window has a floor set by the reconstruction and not by the
   ! wind: the certified states sit at 4.1e-05 to 1.8e-04 and no state
   ! reaches zero.  A weight that reached one only at d_window = 0 would
   ! therefore be below one on EVERY state, and a converged wind would carry
   ! a few parts in ten thousand of the branch it does not belong to.
   real*8, save :: base_wind_window_spread  = 2.0d-3

   ! THE AMPLITUDE THE WINDOW'S SAY IS REGULARIZED ON, and why a shape
   ! criterion alone is not admissible inside a residual.
   !
   ! d_window is scale free: it is unchanged when the whole window flux is
   ! multiplied by a constant.  A weight built on it alone therefore has no
   ! limit at the zero window: along a uniform flux eps F the weight is 1 at
   ! every eps and the window decides, while along a flux of the same
   ! amplitude and a spread above the gate the weight is 0 at every eps and
   ! the first interior cell decides.  Both paths reach the same state, so
   ! the face density has two limits there.  MEASURED before this constant
   ! existed (item L26, sections 4 and 5): the two limits differ by 37.4 per
   ! cent of the face density on the two LHS 1140 b wind states and by 0.41
   ! per cent on a cold start, and the gap equals |w_i - 1/2| times the
   ! distance between the two isentropes the branch mixes, to the last
   ! printed digit.
   !
   ! So the window's say carries an AMPLITUDE factor as well as a shape
   ! factor,
   !
   !     s_wind = A(M_wind) C(d_window) ,
   !
   ! with A a cubic smoothstep that is 0 at M_wind = 0 and 1 at
   ! |M_wind| >= this constant.  A vanishes QUADRATICALLY in M_wind (the
   ! smoothstep's derivative vanishes at both ends), so s_wind and its
   ! derivative both go to zero as the window flux goes to zero ALONG EVERY
   ! DIRECTION, whatever the shape factor is doing, and the first interior
   ! cell takes over smoothly.  The zero window is then an ordinary point of
   ! the closure and not a special case of it.
   !
   ! THE VALUE IS A REGULARIZATION AND NOT A GATE, and it is set so that it
   ! acts on no state the code meets.  A is saturated at 1 for
   ! |M_wind| >= 1e-12.  MEASURED face-mapped window Mach numbers: 1.13e-08
   ! on the certified 0.03 LHS 1140 b state, 7.49e-09 on the stalled 0.02
   ! transient, 3.75e-06 on a cold start of the same case, and READ from
   ! item L21, 3.7e-09 at the cold start of the wasp_full regression case,
   ! which is the smallest any state has been measured at.  1e-12 stands
   ! 3700 times below that and 1e4 below the certified operating point, and
   ! the rounding floor of M_wind on a window of a few hundred cells is
   ! sixteen decades below the flux itself, so A is exactly 1 on every state
   ! and the closure's behaviour on them is C alone.
   !
   ! EXHALE_BASE_WIND_AMPLITUDE=<value> overrides it, which is how the
   ! sensitivity of the answer to this scale is measured;
   ! EXHALE_BASE_WIND_SPREAD=<value> overrides the threshold above.
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
   ! item, not a setting: see the decision section of
   ! docs/phaseC_characteristic_base_bc.md.
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
   ! (report_base_face_state).  Item L21.
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
   write(*,'(A,ES13.6,A,ES10.3,A,ES10.3,A)')                              &
        '     its weight in the branch s =',                              &
        base_face_swind_last, '  (s = A C; C = 1 at d <=',                &
        base_wind_window_spread, ' and 0 at twice it, A = 1 at |M_wind|'//&
        ' >=', base_wind_window_amplitude, ' and 0 at zero)'
   write(*,'(A,ES13.6,A,ES13.6)') '   face state:           v_b =',        &
        base_face_vb_last*v0, ' cm/s,  M_b =', base_face_mach_last
   write(*,'(A,ES10.3,A,ES10.3,A)') '   branch weight w_rev =',            &
        base_face_blend_last, '  (handover width ', base_face_mach_blend,  &
        ' in face Mach; 0 = the reservoir states the entropy,'
   write(*,'(A)') '     1 = the interior does)'
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
   real*8, intent(in) :: p_level, T_level, nhat_level, r_level
   base_reservoir_p    = p_level
   base_reservoir_T    = T_level
   base_reservoir_nhat = nhat_level
   r_base_level        = r_level
   end subroutine set_base_reservoir

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
   ! the two ghost cells, i.e. |ln(T_b/T_a)| < (gamma-1)/gamma * Ng*dr/H,
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
   ! Used to hand the entropy source from the reservoir to the interior
   ! across the reversal branch, and the whole face state from the reservoir
   ! to the interior across the supersonic-outflow branch.
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
   ! The control experiments of items L21 and L26, read once.  Nothing else
   ! reads these variables and a production run sets none of them.
   character(len=32) :: env
   integer :: ios
   real*8  :: v
   if (base_branch_options_read) return
   base_branch_options_read = .true.
   env = ' '
   call get_environment_variable('EXHALE_BASE_BRANCH_ON_CELL1', env)
   if (trim(env) .eq. '1') base_branch_on_wind_flux = .false.
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
   ! usable numbers: the quantities the entropy handover is keyed on.
   real*8, intent(in)  :: F_wind, d_window
   logical, intent(in) :: have_F
   real*8, intent(out) :: Wface(3), Wi(3)
   real*8 :: rho_i, T_i, p_i, v_i, c_i, M_i, M_wind, s_wind, w_i
   real*8 :: rho_res, T_res, p_res
   real*8 :: rho_b, v_b, p_b, c_b, rho_rev
   real*8 :: w_rev, w_out, rb, r1

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
   call continue_hydrostatic_isentrope(0, base_reservoir_nhat,            &
                                       r_base_level,                      &
                                       base_reservoir_p                   &
                                       /(base_reservoir_nhat              &
                                         *base_reservoir_T),              &
                                       base_reservoir_T, rb,              &
                                       rho_res, T_res)
   p_res = base_reservoir_nhat*rho_res*T_res

   ! ---- the two reservoir conditions, and the one interior relation ----
   !
   ! Pressure is a reservoir condition in every subsonic branch.  Entropy is
   ! one only while gas enters: on reversal the entropy at the face is the
   ! interior's, advected out, and the reservoir has nothing to say about it.
   ! The handover is the weight below.  (The Riemann solver downstream
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

   ! THE DIRECTION THE HANDOVER IS DECIDED BY (see base_branch_on_wind_flux).
   ! The conserved flux mapped onto this face by the interior density, which
   ! is the same continuation v_i is built from -- only the flux it carries
   ! is read where the cell-centred product is a flux.
   ! THE WIND DECIDES ONLY WHEN IT HAS SOMETHING TO SAY, AND IT TAKES OVER
   ! SMOOTHLY.  A cold start's wind window carries no wind: the isothermal
   ! initial condition is at rest there, so the mean flux is whatever the
   ! first steps put into it and its sign is not a statement about the base.
   ! The wind is therefore consulted with a weight that vanishes where the
   ! face Mach number it implies vanishes, and the first interior cell
   ! answers there, which is what this boundary did before.
   !
   ! THE HANDOVER IS C1 AND NOT A THRESHOLD, and that is the correction of
   ! the Codex review of 2026-09-15.  A test `if |M_wind| > blend then
   ! M_branch = M_wind else M_i` jumps between two numbers OF OPPOSITE SIGN
   ! as |M_wind| crosses the threshold -- on this very base M_i is -8.1e-07
   ! and M_wind is +4.0e-07 -- so the entropy source of the face, and with it
   ! the ghost density, would step discontinuously and the residual would
   ! carry a jump the Newton solver cannot differentiate.  The weight below
   ! is the same cubic smoothstep the branch itself uses, read as a function
   ! of |M_wind|/blend: 0 at 0, 1/2 at the blend width, 1 at twice it.  Its
   ! derivative vanishes at both ends, so the |.| does not put a kink at
   ! M_wind = 0 either: d(s M_wind)/dM_wind -> 0 from both sides.
   !
   ! ONE BOUNDARY FOR BOTH PATHS.  This was for a time keyed on whether the
   ! evaluation was a stationary one, so that the marching goldens would not
   ! move; that was wrong and the measurement says so (item R1 of the review
   ! of 2026-09-15, docs/lhs1140b_stationary_L21_20260915.md section 11).  A
   ! state is a steady state of ONE operator: with two boundaries,
   ! R_stat(U*) = 0 says nothing about R_time(U*), and on the certified
   ! LHS 1140 b state the two differ by NINE DECADES -- the same state reads
   ! mass 2.17e-09 and energy 2.13e-08 under the flux-keyed face and 1.31 and
   ! 1.05 under the local one, because the local face turns the base into an
   ! OUTFLOW carrying -1.19 F_wind while the wind above it carries +1.000.
   ! So the rule below does not ask which route is evaluating.  It asks only
   ! whether the wind window carries a flux that says anything (have_F, and
   ! the weight s, which vanishes with |M_wind|); where it does not, the
   ! first interior cell answers, which is the cold start and the early
   ! transient.
   !
   ! WHAT REMAINS, stated rather than hidden: the residual of the base cells
   ! depends on cells at r >= r_flux, which the banded preconditioner of the
   ! Newton solve does not contain.  The Krylov products see it through the
   ! finite-difference action, the preconditioner does not, and that is one
   ! more term in the model's error that item L17 already measures at a
   ! factor 27.  A small boundary block in the preconditioner is the obvious
   ! answer and is recorded as a proposal, not made here.
   ! WHAT IS BLENDED IS THE WEIGHT AND NOT THE MACH NUMBER.  Blending the
   ! two Mach numbers first would hand the weight a quantity that carries
   ! (1 - s) M_i, and |M_i| can be hundreds of blend widths at a base inside
   ! a collocated mode, so the blended Mach sweeps hundreds of widths while s
   ! moves by a per cent: continuous, but with a derivative of order
   ! |M_i|/blend, which for a Newton solve is nearly as bad as the jump it
   ! was meant to remove (MEASURED: with M_i = -547 blend, the blended Mach
   ! crossed zero between two samples 3 per cent apart in the wind speed).
   ! The weight of each discriminant is bounded in [0,1] by construction, so
   ! blending THEM is bounded too, and it is still C1: s vanishes with a
   ! vanishing derivative at |M_wind| = 0, so the |.| introduces no kink.
   w_i    = characteristic_branch_weight(M_i/base_face_mach_blend)
   M_wind = 0.0d0
   s_wind = 0.0d0
   w_rev  = w_i
   if (base_branch_on_wind_flux .and. have_F) then
      M_wind = F_wind/(rho_i*rb*rb)/c_i
      ! THE WINDOW'S SAY IS AN AMPLITUDE TIMES A SHAPE, s = A C.  C asks
      ! whether the window carries ONE flux and is scale free; A asks whether
      ! it carries a flux at all and vanishes quadratically with it, so the
      ! product and its derivative vanish as the window flux goes to zero
      ! along EVERY direction and the zero window is an ordinary point.  The
      ! design note at base_wind_window_amplitude carries the measurement
      ! that A is saturated at 1 on every state the code meets, so on those
      ! states this is C alone.
      s_wind = characteristic_branch_weight(1.0d0 - 2.0d0*abs(M_wind)      &
                                            /base_wind_window_amplitude)   &
             * characteristic_branch_weight(                               &
                  (2.0d0*d_window - 3.0d0*base_wind_window_spread)         &
                  /base_wind_window_spread)
      w_rev  = s_wind*characteristic_branch_weight(M_wind                  &
                                            /base_face_mach_blend)         &
             + (1.0d0 - s_wind)*w_i
   endif
   if (w_rev .gt. 0.0d0) n_base_reversal_evals = n_base_reversal_evals + 1
   rho_rev = isentropic_density_at_pressure(1, rho_i, p_i, T_i, p_b)
   rho_b   = (1.0d0 - w_rev)*rho_res + w_rev*rho_rev

   ! (C-): the outgoing acoustic invariant of the first interior state.
   v_b = v_i + (p_b - p_i)/(rho_i*c_i)

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

   ! Sound speed of the face state, at the BASE composition. On the reversal
   ! branch the density came from the interior isentrope and cell 1's
   ! composition would be the consistent one; the difference is confined to
   ! |M_i| < base_face_mach_blend, where the two isentropes agree to the
   ! boundary's own error, and c_b is read only by the Mach cap below and by
   ! the diagnostic.
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

   subroutine reservoir_face_state(Wface)
   ! The reservoir carried to the face, at rest. This is the whole boundary
   ! when the interior has no admissible state to supply a characteristic
   ! from, and it is the state every other branch starts from.
   real*8, intent(out) :: Wface(3)
   real*8 :: rho_res, T_res
   call continue_hydrostatic_isentrope(0, base_reservoir_nhat,            &
                                       r_base_level,                      &
                                       base_reservoir_p                   &
                                       /(base_reservoir_nhat              &
                                         *base_reservoir_T),              &
                                       base_reservoir_T, r_edg(0),        &
                                       rho_res, T_res)
   Wface(1) = rho_res
   Wface(2) = 0.0d0
   Wface(3) = base_reservoir_nhat*rho_res*T_res
   end subroutine reservoir_face_state

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

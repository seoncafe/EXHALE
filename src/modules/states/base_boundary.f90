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
   ! larger.  1e-6 therefore leaves every intended operating point entirely
   ! inside the inflow branch, and the residual jump it smooths is the
   ! difference between two isentropes that coincide at a converged state.
   !
   ! Its influence is measured, not assumed: see the design note's section 7.
   real*8, parameter :: base_face_mach_blend = 1.0d-6

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

   contains

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

   subroutine characteristic_base_face_state(W1, nhat1, T1, Wface, Wi)
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
   real*8, intent(out) :: Wface(3), Wi(3)
   real*8 :: rho_i, T_i, p_i, v_i, c_i, M_i
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

   w_rev = characteristic_branch_weight(M_i/base_face_mach_blend)
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
   ! |M_i| < base_face_mach_blend = 1e-6, where the two isentropes agree to the
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
   real*8 :: nhat1, T1, Wi(3)

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

   call characteristic_base_face_state(W(:,1), nhat1, T1, Wface, Wi)
   call base_ghost_averages(Wface, Wghost, Wface_lower)

   end subroutine base_boundary_states

   ! End of module
   end module base_boundary

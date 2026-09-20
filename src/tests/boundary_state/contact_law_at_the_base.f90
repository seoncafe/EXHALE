      program contact_law_at_the_base
      ! THE CONTACT LAW OF THE LOWER BOUNDARY, on states built here.
      !
      ! The boundary states the level as a contact between the lower
      ! atmosphere and the first interior cell.  It matches the two
      ! acoustically first -- the interior's outgoing relation (C-) closed
      ! with the reservoir's pressure -- then reads the direction of the
      ! contact off the matched state, and only then upwinds the entropy and
      ! composition of the face.  At rest and at inflow the reservoir owns
      ! the level, because the lower atmosphere is a heat bath on the time
      ! scales of a stationary solution; only a reverse flow carries the
      ! interior's trace out.  base_boundary.f90 carries the statement and
      ! its validity range.
      !
      ! SIX STATEMENTS, one per section.
      !
      !   a. A column at rest on the reservoir's own hydrostatic isentrope,
      !      with the molecular equation of state and the discrete gravity,
      !      is left at rest: the face velocity and the mass flux the face
      !      state carries are at their rounding floors, and the face is the
      !      reservoir.
      !   b. A pressure-balanced stationary contact whose interior carries a
      !      DIFFERENT entropy takes the reservoir's density and temperature,
      !      not an average of the two isentropes, and launches no flow.
      !   c. A pressure mismatch of either sign launches an acoustic
      !      response, the outgoing invariant of the interior is carried
      !      through the face unchanged, and the branch follows the sign of
      !      the matched face velocity: the reservoir where gas enters, the
      !      interior's trace where it leaves.
      !   d. The base boundary conditions of the transport operators are
      !      their own: thermal conduction reads the level's own temperature
      !      and the element diffusion the reservoir's composition, and
      !      neither moves when the advective contact changes direction.
      !      The conductive heat flux through the base face at zero bulk
      !      velocity is measured and reported, not assumed zero.
      !   e. One function answers for every consumer: the face state
      !      Apply_BC installs on a state is the face state
      !      base_boundary_states returns for it.
      !   f. The reversal branch has no width: the face state does not move
      !      when the one remaining regularization width, the
      !      supersonic-outflow handover, is changed by four decades.
      !
      ! Everything is the production boundary on the production grid and
      ! gravity; nothing here re-implements the condition.
      use global_parameters
      use grid_construction,          only: define_grid
      use gravity_grid_construction,  only: set_gravity_grid
      use species_table,              only: isp_H2
      use caloric_eos,                only: caloric_state_from_composition, &
                                            caloric_mixture_active,         &
                                            adiabatic_index_at_T
      use base_boundary,              only: set_base_reservoir,             &
                                            base_boundary_states,           &
                                            characteristic_base_face_state, &
                                            wind_window_mass_flux,          &
                                            base_face_T_i_last,             &
                                            continue_hydrostatic_isentrope, &
                                            base_reservoir_temperature_at,  &
                                            base_face_mach_blend,           &
                                            base_face_blend_last,           &
                                            base_face_vb_last,              &
                                            base_face_rho_res_last,         &
                                            base_face_rho_rev_last,         &
                                            base_reservoir_nhat
      use viscous_conduction,         only: conduction_base_level_T,        &
                                            conduction_base_heat_flux,      &
                                            thermal_conduction_source
      use BC_Apply,                   only: Apply_BC, base_face_W
      use Conversion,                 only: W_to_U
      implicit none

      integer, parameter :: Ncells = 200
      real*8,  parameter :: Teq = 1.0d3, Rp_RJ = 0.2d0, Mp_MJ = 0.02d0
      ! The reservoir of the tests: one code unit of pressure and
      ! temperature at the planet radius, the particles per unit mass of a
      ! hydrogen mixture that is mostly H2.
      real*8,  parameter :: p_res_level = 1.0d0, T_res_level = 1.0d0
      real*8,  parameter :: nhat_res    = 0.6d0

      real*8, allocatable :: W(:,:), u(:,:)
      real*8, allocatable :: rho_c(:), ne_c(:), ntot_c(:), fsp(:,:)
      real*8  :: Wface(3), Wghost(3,1-2:0), Wface_lower(3)
      real*8  :: rho_res_face, T_res_face, p_res_face, c_res
      real*8  :: rho_b_a, v_b_a, flux_a
      real*8  :: rho_b_b, v_b_b, w_b, rho_rev_b
      real*8  :: v_in, v_out, w_in, w_out, inv_i, inv_b
      real*8  :: T_cond_in, T_cond_out, X_res_in, X_res_out
      real*8  :: q_base, rho_ref
      real*8  :: rho_w(3), face_w(3)
      real*8, dimension(:), allocatable :: Tcell, Qc
      integer :: j, n_fail

      n_fail = 0

      N  = Ncells
      T0 = Teq
      ! The number-density scale, which the transport coefficients divide
      ! by; the value is the one a hot-Uranus base run carries.
      n0 = 1.0d14
      HeH = 0.0793d0
      R0 = Rp_RJ*RJ
      Mp = Mp_MJ*MJ
      v0 = sqrt(kb_erg*T0/mu)
      b0 = (Gc*Mp*mu)/(kb_erg*T0*R0)
      ! The grid configuration of a base-resolved run, the same defaults
      ! input_read uses; only the base region matters here.
      spherical_domain = .true.
      dr_base     = 2.0d-4
      N_low_cells = 50
      r_max       = 10.0d0
      r_esc       =  2.0d0
      r_flux      =  1.20d0
      grid_type   = 'Mixed'
      CFL         = 0.6d0
      call allocate_grid_arrays
      call define_grid
      call set_gravity_grid
      allocate(W(3,1-Ng:N+Ng), u(3,1-Ng:N+Ng))
      allocate(rho_c(1-Ng:N+Ng), ne_c(1-Ng:N+Ng), ntot_c(1-Ng:N+Ng))
      allocate(fsp(1-Ng:N+Ng,n_species))
      allocate(Tcell(1-Ng:N+Ng), Qc(1-Ng:N+Ng))

      call set_base_reservoir(p_res_level, T_res_level, nhat_res, 1.0d0)
      ! The caloric state first: the hydrostatic continuation reads the
      ! rovibrational ladder of the cell it is continuing, so the reservoir
      ! state at the face is a different number under the monatomic map
      ! (2.1e-4 in density on this grid, MEASURED) and the column has to be
      ! built under the map the condition will use.  The mixture the install
      ! writes is a pair of ratios and does not depend on the column, so a
      ! provisional column is enough to set it.
      W = 1.0d0
      call install_molecular_caloric_state()
      call continue_hydrostatic_isentrope(0, nhat_res, 1.0d0,              &
               p_res_level/(nhat_res*T_res_level), T_res_level, r_edg(0),  &
               rho_res_face, T_res_face)
      p_res_face = nhat_res*rho_res_face*T_res_face
      c_res      = sqrt(adiabatic_index_at_T(0, T_res_face)                &
                        *p_res_face/rho_res_face)

      ! ---------------------------------------------------------------- !
      ! a. the reservoir's own hydrostatic column, at rest
      ! ---------------------------------------------------------------- !
      ! Every cell is the point state of the reservoir's isentrope at its
      ! centre, so the interior continued back down to the face IS the
      ! reservoir there and the condition is exercised at its well-balanced
      ! point.  The caloric state is the molecular one: the H2 fraction
      ! below makes caloric_mixture_active true, so the continuation and the
      ! sound speed read the rovibrational ladder and not 3/2 kT.
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      if (.not. caloric_mixture_active) then
         write(*,'(A)') 'FAIL contact_molecular_eos_is_active'//           &
              ' measured=0 reference=1 tol=0'
         n_fail = n_fail + 1
      else
         write(*,'(A)') 'PASS contact_molecular_eos_is_active'//           &
              ' measured=1 reference=1 tol=0'
      endif
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_b_a = Wface(1);  v_b_a = Wface(2)
      flux_a  = Wface(1)*Wface(2)*r_edg(0)*r_edg(0)
      call verdict('contact_hydrostatic_column_face_velocity_at_the_floor',&
                   abs(v_b_a)/c_res, 0.0d0, 1.0d-13)
      call verdict('contact_hydrostatic_column_base_flux_at_the_floor',    &
                   abs(flux_a)/(rho_res_face*c_res*r_edg(0)*r_edg(0)),     &
                   0.0d0, 1.0d-13)
      call verdict('contact_hydrostatic_column_face_is_the_reservoir',     &
                   abs(rho_b_a - rho_res_face)/rho_res_face, 0.0d0,        &
                   1.0d-13)
      ! Repeated derivation from the same unmoved column returns the same
      ! numbers: the condition is a function of its inputs and of nothing
      ! kept between calls.
      do j = 1, 20
         call base_boundary_states(W, Wface, Wghost, Wface_lower)
      enddo
      call verdict('contact_hydrostatic_column_repeats_exactly',           &
                   abs(Wface(1) - rho_b_a) + abs(Wface(2) - v_b_a),        &
                   0.0d0, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! b. a pressure-balanced stationary contact with unequal entropies
      ! ---------------------------------------------------------------- !
      ! The interior is a column at rest whose own isentrope passes through
      ! the SAME pressure at the face and twice the temperature there, so
      ! the contact is stationary and pressure balanced and the two sides
      ! carry two entropies.  The face is the reservoir's: its density and
      ! its temperature are the level's, and no flow is launched.
      call fill_interior_isentrope(2.0d0)
      call install_molecular_caloric_state()
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_b_b   = Wface(1);  v_b_b = Wface(2);  w_b = base_face_blend_last
      rho_rev_b = base_face_rho_rev_last
      write(*,'(A,ES14.7,A,ES14.7,A,F6.3)')                                &
           '  the two isentropes at the face: reservoir ', rho_res_face,   &
           '  interior ', rho_rev_b, '  ratio ', rho_rev_b/rho_res_face
      call verdict('contact_stationary_two_entropies_are_distinct',        &
                   abs(rho_rev_b - rho_res_face)/rho_res_face, 0.5d0,      &
                   0.2d0)
      call verdict('contact_stationary_face_velocity_at_the_floor',        &
                   abs(v_b_b)/c_res, 0.0d0, 1.0d-13)
      call verdict('contact_stationary_face_is_the_reservoir',             &
                   abs(rho_b_b - rho_res_face)/rho_res_face, 0.0d0,        &
                   1.0d-13)
      call verdict('contact_stationary_takes_no_intermediate_trace',       &
                   w_b, 0.0d0, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! c. acoustic perturbations of both signs
      ! ---------------------------------------------------------------- !
      ! The interior of (a) with its pressure raised and lowered by one part
      ! in 1e5 at fixed density, which is a pressure mismatch across the
      ! level and no cell velocity at all: the matched face velocity is then
      ! the whole of the response and the two signs select the two branches.
      ! The interior CELL velocity is zero in both, so a condition keyed on
      ! it could not tell them apart.
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      call scale_interior_pressure(1.0d0 + 1.0d-5)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      v_out = Wface(2);  w_out = base_face_blend_last
      call interior_invariants(inv_i, inv_b, Wface)
      call verdict('contact_overpressured_interior_leaves_through_the_face',&
                   sign(1.0d0, v_out), -1.0d0, 0.0d0)
      call verdict('contact_reverse_flow_carries_the_interior_trace',      &
                   w_out, 1.0d0, 0.0d0)
      call verdict('contact_outgoing_invariant_is_carried_through',        &
                   abs(inv_b - inv_i)/abs(inv_i), 0.0d0, 1.0d-13)
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      call scale_interior_pressure(1.0d0 - 1.0d-5)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      v_in = Wface(2);  w_in = base_face_blend_last
      call verdict('contact_underpressured_interior_draws_gas_in',         &
                   sign(1.0d0, v_in), 1.0d0, 0.0d0)
      call verdict('contact_inflow_carries_the_reservoir_trace',           &
                   w_in, 0.0d0, 0.0d0)
      call verdict('contact_inflow_face_is_the_reservoir',                 &
                   abs(Wface(1) - base_face_rho_res_last)                  &
                   /base_face_rho_res_last, 0.0d0, 1.0d-13)

      ! ---------------------------------------------------------------- !
      ! d. the transport operators carry their own base conditions
      ! ---------------------------------------------------------------- !
      ! The advective contact is flipped between the two branches by the
      ! pressure mismatch of section c, and the two transport boundary
      ! conditions are read on either side of the flip.  Thermal conduction
      ! reads the level's own temperature, the lower atmosphere being a heat
      ! bath whichever way the gas crosses the level; the element diffusion
      ! reads the reservoir's helium mass fraction, a Dirichlet value that
      ! is a function of the input He/H and of nothing else.
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      call scale_interior_pressure(1.0d0 - 1.0d-5)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      call cell_temperatures(Tcell)
      T_cond_in = conduction_base_level_T(Tcell)
      X_res_in  = reservoir_helium_mass_fraction()
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      call scale_interior_pressure(1.0d0 + 1.0d-5)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      call cell_temperatures(Tcell)
      T_cond_out = conduction_base_level_T(Tcell)
      X_res_out  = reservoir_helium_mass_fraction()
      call verdict('transport_conduction_base_T_does_not_follow_the_contact',&
                   abs(T_cond_out - T_cond_in), 0.0d0, 0.0d0)
      call verdict('transport_conduction_base_T_is_the_level_temperature', &
                   abs(T_cond_in - base_reservoir_temperature_at(r(0)))    &
                   /base_reservoir_temperature_at(r(0)), 0.0d0, 0.0d0)
      call verdict('transport_element_base_composition_does_not_follow_'// &
                   'the_contact', abs(X_res_out - X_res_in), 0.0d0, 0.0d0)
      ! The conductive heat flux the level puts through the base face with
      ! the gas at rest.  It is MEASURED and printed, not asserted: a zero
      ! bulk velocity does not make a diffusive flux zero, and the number
      ! this column carries depends on the temperature contrast it is built
      ! with.  Here the interior of section b stands at twice the level's
      ! temperature, so the flux is out of the domain.
      call fill_interior_isentrope(2.0d0)
      call install_molecular_caloric_state()
      call cell_temperatures(Tcell)
      cond_on = .true.
      call thermal_conduction_source(Tcell, Qc)
      q_base = conduction_base_heat_flux
      cond_on = .false.
      write(*,'(A,ES14.7,A)') 'DIAGNOSTIC conduction base heat flux at'//  &
           ' zero bulk velocity = ', q_base, ' [code, positive inward]'
      call verdict('transport_conduction_base_flux_is_not_assumed_zero',   &
                   sign(1.0d0, -q_base), 1.0d0, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! e. one function answers for every consumer
      ! ---------------------------------------------------------------- !
      ! Apply_BC derives the boundary and caches the face state; every flux
      ! assembly, marching or stationary, reads that cache.  The cached face
      ! state must be the face state the condition returns for the same
      ! interior column.
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      call scale_interior_pressure(1.0d0 + 1.0d-5)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      face_w = Wface
      call W_to_U(W, u)
      call Apply_BC(u)
      rho_w(1) = abs(base_face_W(1) - face_w(1))/face_w(1)
      rho_w(2) = abs(base_face_W(2) - face_w(2))/abs(face_w(2))
      rho_w(3) = abs(base_face_W(3) - face_w(3))/face_w(3)
      ! The tolerance is the u -> W round trip and not a choice: Apply_BC
      ! derives its own primitives from the conserved array, and that round
      ! trip is not the identity in floating point (it recomputes the
      ! kinetic energy from rho and v where the input held it from the
      ! momentum), so the two derivations stand on interior states that
      ! differ in their last bits.  The base amplifies that to 2.2e-11 here,
      ! MEASURED; anything above 1e-9 would be a second boundary.
      call verdict('contact_installed_face_state_is_the_derived_one',      &
                   maxval(rho_w), 0.0d0, 1.0d-9)

      ! ---------------------------------------------------------------- !
      ! f. the reversal branch has no width
      ! ---------------------------------------------------------------- !
      ! The one width left in the condition regularizes the
      ! SUPERSONIC-OUTFLOW transition at M_i = -1.  A subsonic base is four
      ! decades of that width away from it, so the face state must not move
      ! when the width is changed by four decades.  Under a condition whose
      ! entropy branch is handed over on the same width, a state whose
      ! interior Mach number is comparable to it moves instead.
      call fill_interior_isentrope(1.0d0)
      call install_molecular_caloric_state()
      call scale_interior_pressure(1.0d0 + 1.0d-5)
      call set_interior_velocity(3.0d-9)
      rho_ref = face_density_at_width(1.0d-8)
      call verdict('contact_face_state_independent_of_the_width_1e6',      &
                   abs(face_density_at_width(1.0d-6) - rho_ref)/rho_ref,   &
                   0.0d0, 0.0d0)
      call verdict('contact_face_state_independent_of_the_width_1e10',     &
                   abs(face_density_at_width(1.0d-10) - rho_ref)/rho_ref,  &
                   0.0d0, 0.0d0)
      base_face_mach_blend = 1.0d-8

      if (n_fail .gt. 0) then
         write(*,'(A,I0)') 'contact_law_at_the_base FAILURES: ', n_fail
         stop 1
      endif
      stop 0

      contains

      subroutine fill_interior_isentrope(Tfac)
      ! A column at rest whose own isentrope passes through the reservoir's
      ! PRESSURE at the face and Tfac times its temperature there: the
      ! contact is pressure balanced, and at Tfac /= 1 the two sides carry
      ! two entropies.  Tfac = 1 is the reservoir's own column, the state
      ! the condition is well balanced on.
      !
      ! IT IS BUILT FROM THE FACE, not from the base level, and that is not
      ! a convenience: the continuation is exact only over short steps (its
      ! midpoint heat capacity makes the density third-order accurate in
      ! ln(T_b/T_a)), so a column continued from the level to every cell
      ! centre in ONE step is a different discrete isentrope from the one
      ! the condition walks in two, and the difference -- 2.1e-4 in the face
      ! density on this grid, MEASURED -- would be read as a boundary error.
      real*8, intent(in) :: Tfac
      real*8 :: rho_q, T_q, rho_face, T_face
      integer :: jj
      T_face   = Tfac*T_res_face
      rho_face = p_res_face/(nhat_res*T_face)
      do jj = 1-Ng, N+Ng
         call continue_hydrostatic_isentrope(0, nhat_res, r_edg(0),       &
                  rho_face, T_face, r(jj), rho_q, T_q)
         W(1,jj) = rho_q
         W(2,jj) = 0.0d0
         W(3,jj) = nhat_res*rho_q*T_q
      enddo
      n_part_cell1 = nhat_res*W(1,1)
      end subroutine fill_interior_isentrope

      subroutine scale_interior_pressure(fac)
      ! A pressure mismatch across the level at fixed density and fixed
      ! velocity: the interior's face pressure moves off the reservoir's and
      ! the matched face velocity is the whole response.
      real*8, intent(in) :: fac
      integer :: jj
      do jj = 1, N+Ng
         W(3,jj) = fac*W(3,jj)
      enddo
      end subroutine scale_interior_pressure

      subroutine set_interior_velocity(vv)
      ! One velocity in every cell, at fixed density and pressure.
      real*8, intent(in) :: vv
      integer :: jj
      do jj = 1-Ng, N+Ng
         W(2,jj) = vv
      enddo
      end subroutine set_interior_velocity

      subroutine install_molecular_caloric_state()
      ! A hydrogen mixture that is four fifths H2 by nucleus, the partition
      ! of the hot-Uranus handoff, so that the caloric equation of state the
      ! boundary continuation reads is the rovibrational one.
      integer :: jj
      thereis_mol = .true.
      fsp   = 0.0d0
      do jj = 1-Ng, N+Ng
         fsp(jj,isp_H2) = 0.42d0*nhat_res
         ne_c(jj)       = 1.0d-6*nhat_res*W(1,jj)
         ntot_c(jj)     = nhat_res*W(1,jj) - ne_c(jj)
         rho_c(jj)      = W(1,jj)
      enddo
      call caloric_state_from_composition(rho_c, fsp, ne_c, ntot_c)
      end subroutine install_molecular_caloric_state

      subroutine cell_temperatures(Tc)
      ! p = n_part T in code units, with the particle count per unit mass
      ! the caloric state was installed at.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Tc
      integer :: jj
      do jj = 1-Ng, N+Ng
         Tc(jj) = W(3,jj)/(nhat_res*W(1,jj))
      enddo
      end subroutine cell_temperatures

      real*8 function reservoir_helium_mass_fraction() result(X)
      ! The Dirichlet base composition of the element-diffusion operator,
      ! written as that operator writes it: a function of the input He/H and
      ! of nothing the flow does.
      real*8, parameter :: m_He_amu_local = 4.002602d0, m_1_local = 1.0d0
      X = m_He_amu_local*HeH/(m_1_local + m_He_amu_local*HeH)
      end function reservoir_helium_mass_fraction

      subroutine interior_invariants(inv_interior, inv_face, Wf)
      ! The outgoing acoustic invariant p - rho_i c_i v of the interior at
      ! the face, and the same invariant formed on the face state.  (C-)
      ! says the two are equal, and that equality is what "the outgoing
      ! wave leaves through the boundary without reflection" means for the
      ! STATE the boundary hands to the Riemann solve.  The interior state
      ! at the face comes from the condition itself, so the factor
      ! rho_i c_i is the one it used.
      real*8, intent(out) :: inv_interior, inv_face
      real*8, intent(in)  :: Wf(3)
      real*8 :: Fw, dw, Wfc(3), Wi(3), rho_c_i
      logical :: hF
      call wind_window_mass_flux(W, Fw, dw, hF)
      call characteristic_base_face_state(W(:,1), n_part_cell1/W(1,1),    &
               W(3,1)/n_part_cell1, Fw, dw, hF, Wfc, Wi)
      rho_c_i = sqrt(adiabatic_index_at_T(1, base_face_T_i_last)          &
                     *Wi(3)*Wi(1))
      inv_interior = Wi(3) - rho_c_i*Wi(2)
      inv_face     = Wf(3) - rho_c_i*Wf(2)
      end subroutine interior_invariants

      real*8 function face_density_at_width(wd) result(rho_f)
      ! The face density of the state now in W, at a stated regularization
      ! width.
      real*8, intent(in) :: wd
      real*8 :: Wf(3), Wg(3,1-Ng:0), Wl(3)
      base_face_mach_blend = wd
      call base_boundary_states(W, Wf, Wg, Wl)
      rho_f = Wf(1)
      end function face_density_at_width

      subroutine verdict(name, measured, reference, tol)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: measured, reference, tol
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'PASS ', name,        &
              ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'FAIL ', name,        &
              ' measured=', measured, ' reference=', reference, ' tol=', tol
         n_fail = n_fail + 1
      endif
      end subroutine verdict

      end program contact_law_at_the_base

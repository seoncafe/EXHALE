      module steady_residual_mod
      ! Finite-volume STEADY residual  R = du/dt  (zero at a true steady
      ! state), shared by the EXHALE_RESIDUAL diagnostic, the in-loop
      ! convergence monitor, and the steady-state Newton/PTC solver.
      !
      !   R(:,1) = dF - S                          (mass)
      !   R(:,2) = dF - S - F_mu                   (momentum)
      !   R(:,3) = dF_E - S_E - (heat-cool)
      !                       - (w F_mu + q_mu + conduction)   (energy)
      !
      ! dF, S come from the existing Reconstruct + RK_rhs. heat/cool are
      ! supplied by the caller, which decides how they were obtained: the
      ! in-loop monitor reuses the current step's values; the diagnostic
      ! and the Newton residual first call ioniz_eq (local
      ! ionization-equilibrium elimination) to get heat/cool consistent
      ! with u.
      !
      ! The molecular-transport terms are the operator-split stage that
      ! viscous_conduction_step relaxes in the marching loop, evaluated by
      ! the SAME routine (viscous_conduction_sources) from the same
      ! tridiagonal operator, so the Newton solver solves exactly the system
      ! the marching relaxes. They vanish identically unless "Viscosity:" /
      ! "Conduction:" are set. The caller passes n_part = n_tot + n_e (the
      ! adimensional particle count) rather than T itself, so that T = p/n_part
      ! stays a function of the unknowns and the temperature dependence of
      ! the conduction operator is picked up by the residual's linearization.
      !
      ! The reconstruction scheme is whatever use_plm/use_weno3/rec_method
      ! currently select; callers set WENO3 for a production residual.

      use global_parameters
      use Conversion, only: U_to_W
      use caloric_eos, only: pressure_from_energy_density,          &
                             adiabatic_index_from_state
      use Reconstruction_step
      use RK_integration
      use viscous_conduction, only: transport_active,                   &
                                    viscous_conduction_sources
      use hydrodynamic_rows, only: generic_precision_rows_arm,          &
                                   ARM_QUADRUPLE, ARM_GENERIC_DOUBLE,   &
                                   hydrodynamic_rows_in_quadruple_precision, &
                                   hydrodynamic_rows_in_double_precision
      use ionization_equilibrium, only: ieq_sweep_state_kind,           &
                                        ieq_state_marching

      implicit none
      private
      public :: assemble_residual, residual_norms,                   &
                reconstruction_continuation_rhs,                     &
                write_residual_breakdown,                            &
                residual_row_scale, carrier_row_scale,               &
                relnorm_over_cells, state_scales_of_cell,            &
                flux_spread_of_state, steady_gates_met,              &
                flux_spread_above_radius,                            &
                n_cells_without_chemical_root,                       &
                row_terms_describe_state,                            &
                face_mass_flux_of_state

      ! THE TERMS EACH CONSERVATION ROW IS BUILT FROM, as assemble_residual
      ! last produced them, together with the state they belong to. The three
      ! row scales below are these; nothing else reads them.
      !
      !   face_mass_flux_r2(j)     F_{j+1/2} r_{j+1/2}^2     (mass, SIGNED:
      !                            the flux gate needs the sign, so a flow
      !                            that reverses inside the window is a huge
      !                            spread and not a small one; the row scale
      !                            takes the magnitude where it needs one)
      !   momentum_largest_term(j) max(|dF_2|, |S_2|)        (momentum)
      !   energy_largest_term(j)   max(|dF_3|, |S_3|, heat, cool)  (energy)
      !
      ! plus the operator-split viscous/conduction sources where those are
      ! active. The state is kept so that a caller asking for the scale of a
      ! different one is answered from ITS terms (refresh_row_terms) rather
      ! than from whatever the caller last evaluated.
      real*8, dimension(:),   allocatable :: face_mass_flux_r2
      ! The same face mass flux without the area factor, F_rho(j) at
      ! r_edg(j).  It is what every species row rides on: a transported
      ! species crosses a face carrying F_rho times its face mass fraction
      ! (species_face_flux.f90), so a stationary species row has to read the
      ! mass flux of the very state whose mass row was assembled.
      real*8, dimension(:),   allocatable :: face_mass_flux
      real*8, dimension(:),   allocatable :: momentum_largest_term
      real*8, dimension(:),   allocatable :: energy_largest_term
      real*8, dimension(:,:), allocatable :: state_of_row_terms

      contains

      ! ------------------------------------------------------!

      subroutine reconstruction_continuation_rhs(u_in, WL, WR, dF, S)
      ! Right-hand side of the homotopy between the two discretizations,
      !
      !     R_lambda(u) = (1 - lambda) R_PLM(u) + lambda R_WENO3(u),
      !
      ! returned as its flux-difference and source parts so that the marching
      ! stages and the steady residual both use it and therefore solve the
      ! same equation at every lambda. Everything PLM and WENO3 disagree about
      ! is inside the two evaluations: the face reconstruction, the
      ! extrapolated reconstructed-face boundary states of Rec_BC, and the
      ! pressure term of the momentum equation (inside the flux plus a
      ! geometric source for PLM, a face-pressure difference for WENO3). The
      ! conservative ghost fill is blended by Apply_BC_W under the same
      ! lambda.
      !
      ! lambda <= 0 and lambda >= 1 take a single evaluation and reproduce the
      ! two single-scheme operators to the bit; only the open interval costs
      ! two right-hand sides, and only for the length of the ramp.
      !
      ! For a marching stage this makes the update the convex combination of
      ! the two schemes' forward-Euler stages. The admissible set rho > 0,
      ! rho e > 0 is convex, so a stage is admissible whenever both
      ! single-scheme stages are: the continuation cannot manufacture a
      ! positivity failure that neither scheme has.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u_in
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL, WR, dF, S
      real*8, dimension(3,1-Ng:N+Ng) :: WLw, WRw, dFw, Sw, ff_plm
      real*8, dimension(1-Ng:N+Ng)   :: fp_plm
      logical :: sav_plm, sav_weno
      character(len=:), allocatable :: sav_method
      real*8  :: lam, om

      if (.not. recon_lambda_on) then
         call Reconstruct(u_in, WL, WR)
         call RK_rhs(u_in, WL, WR, dF, S)
         return
      endif

      lam        = recon_lambda
      om         = 1.0d0 - lam
      sav_plm    = use_plm
      sav_weno   = use_weno3
      sav_method = rec_method

      if (lam .le. 0.0d0) then
         rec_method = 'PLM';    use_plm = .true.;  use_weno3 = .false.
         call Reconstruct(u_in, WL, WR)
         call RK_rhs(u_in, WL, WR, dF, S)
      else if (lam .ge. 1.0d0) then
         rec_method = 'WENO3';  use_plm = .false.; use_weno3 = .true.
         call Reconstruct(u_in, WL, WR)
         call RK_rhs(u_in, WL, WR, dF, S)
      else
         rec_method = 'PLM';    use_plm = .true.;  use_weno3 = .false.
         call Reconstruct(u_in, WL, WR)
         call RK_rhs(u_in, WL, WR, dF, S)
         ff_plm = face_flux
         fp_plm = face_p

         rec_method = 'WENO3';  use_plm = .false.; use_weno3 = .true.
         call Reconstruct(u_in, WLw, WRw)
         call RK_rhs(u_in, WLw, WRw, dFw, Sw)

         dF = om*dF + lam*dFw
         S  = om*S  + lam*Sw
         WL = om*WL + lam*WLw
         WR = om*WR + lam*WRw
         ! The stored interface fluxes are what the positivity flux
         ! correction rebuilds a cell from; keep them on the same homotopy.
         face_flux = om*ff_plm + lam*face_flux
         face_p    = om*fp_plm + lam*face_p
      endif

      rec_method = sav_method
      use_plm    = sav_plm
      use_weno3  = sav_weno

      end subroutine reconstruction_continuation_rhs

      ! ------------------------------------------------------!

      ! ------------------------------------------------------!

      subroutine assemble_residual(u, n_part, heat, cool, R)
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: n_part
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: R
      ! Local scratch so callers' own WL/WR/dF/S are untouched
      real*8, dimension(3,1-Ng:N+Ng) :: WL, WR, dF, S, W
      real*8, dimension(1-Ng:N+Ng)   :: Tc, Smom, Sene
      integer :: arm

      ! THE ARITHMETIC THE HYDRODYNAMIC ROWS ARE ASSEMBLED IN.  Normally the
      ! production routines, in double.  EXHALE_RESID_QUAD=1 sends the
      ! STATIONARY evaluations -- and only those; the marching stages reach
      ! RK_rhs directly and never come here -- through the
      ! quadruple-precision instantiation of the same kind-generic text, as
      ! the control experiment of where the residual's non-smoothness floor
      ! comes from: that floor is the rounding of this flux assembly, and an
      ! arm that lowers the rounding by eighteen decades and nothing else
      ! separates the rounding from every other candidate. Default off, and
      ! nothing is adopted from it; see the header of module
      ! hydrodynamic_rows.
      arm = 0
      if (ieq_sweep_state_kind .ne. ieq_state_marching)                  &
         arm = generic_precision_rows_arm()
      if (arm .eq. ARM_QUADRUPLE) then
         call hydrodynamic_rows_in_quadruple_precision(u, WL, WR, dF, S)
      else if (arm .eq. ARM_GENERIC_DOUBLE) then
         call hydrodynamic_rows_in_double_precision(u, WL, WR, dF, S)
      else
         call reconstruction_continuation_rhs(u, WL, WR, dF, S)
      endif
      R(1,:) = dF(1,:) - S(1,:)
      R(2,:) = dF(2,:) - S(2,:)
      R(3,:) = dF(3,:) - S(3,:) - (heat - cool)

      Smom = 0.0d0;  Sene = 0.0d0
      if (transport_active()) then
         call U_to_W(u, W)
         Tc = W(3,:)/n_part
         call viscous_conduction_sources(W(2,:), Tc, Smom, Sene)
         R(2,:) = R(2,:) - Smom
         R(3,:) = R(3,:) - Sene
      endif
      ! THE LOWEST GHOST CELL CARRIES NO EQUATION: it has no lower face, so
      ! there is no balance to state there and its row is zero (RK_rhs
      ! defines the same column of dF and S as zero for the same reason).
      ! Written explicitly so that no reader of R can find a heating rate
      ! standing alone in a row that has no fluxes.
      R(:,1-Ng) = 0.0d0
      call store_row_terms(u, dF, S, heat, cool, Smom, Sene)

      end subroutine assemble_residual

      ! ------------------------------------------------------!

      subroutine state_scales_of_cell(jcell, rho_in, mom_in, ene_in,    &
                                      d1, d2, d3, vv, cs)
      ! Characteristic scale of each conserved unknown IN ITS OWN CELL, from
      ! the local state alone, together with the velocity and sound speed it
      ! is built from. THE single definition of the UNKNOWNS' scales: the
      ! steady solver's diagonal scaling (cell_state_scales, steady_newton.f90)
      ! is these three numbers. They are no longer what the convergence
      ! measure divides by -- since section 143 each residual row is divided by
      ! the largest term that row itself contains (residual_row_scale below) --
      ! and that separation is deliberate; residual_row_scale says why.
      !
      !   mass        rho
      !   momentum    |rho v| + rho c_s = rho (|v| + c_s)
      !   energy      E
      !
      ! with c_s = sqrt(gamma_eff p/rho) and p from the caloric EOS of the
      ! cell, so all of it is a function of the cell's own conserved
      ! variables and its composition.
      !
      ! rho and E are positive definite and are their own scales. The momentum
      ! density is not: it passes through zero wherever the flow reverses, and
      ! at the hand-off from marching the whole sub-sonic region below the
      ! stagnation point is still infalling, so |rho v| vanishes there and is
      ! no measure of how large a momentum residual is. The scale that does
      ! not vanish is rho times the fastest characteristic speed of the Euler
      ! system, |v| + c_s: the momentum density the cell carries when moved at
      ! its own signal speed. This is parameter-free -- no floor fraction to
      ! choose -- and it is local, so no single cell can set the scale of the
      ! whole domain.
      integer, intent(in) :: jcell
      real*8, intent(in)  :: rho_in, mom_in, ene_in
      real*8, intent(out) :: d1, d2, d3, vv, cs
      real*8 :: rho, pgas
      rho = rho_in
      if (rho .gt. 0.0d0) then
         pgas = pressure_from_energy_density(jcell, rho,               &
                   ene_in - 0.5d0*mom_in*mom_in/rho)
         cs   = sqrt(adiabatic_index_from_state(jcell, rho,            &
                     max(pgas, 0.0d0))*max(pgas, 0.0d0)/rho)
      else
         ! Not a physical state (the line search rejects these); fall back
         ! to the bare magnitudes rather than taking a root of a negative.
         rho = abs(rho);  cs = 0.0d0
      endif
      vv = mom_in/max(rho, tiny(1.0d0))
      d1 = max(rho,                        tiny(1.0d0))
      d2 = max(abs(mom_in) + rho*cs,       tiny(1.0d0))
      d3 = max(abs(ene_in),                tiny(1.0d0))
      end subroutine state_scales_of_cell

      ! ------------------------------------------------------!

      subroutine store_row_terms(u, dF, S, heat, cool, Smom, Sene)
      ! Keep the terms of each row as RK_rhs and the ionization sweep have
      ! just produced them, and the state they belong to.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u, dF, S
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: heat, cool, Smom, Sene
      integer :: j
      if (.not. allocated(face_mass_flux_r2))                            &
         allocate(face_mass_flux_r2(1-Ng:N+Ng),                          &
                  face_mass_flux(1-Ng:N+Ng),                             &
                  momentum_largest_term(1-Ng:N+Ng),                      &
                  energy_largest_term(1-Ng:N+Ng),                        &
                  state_of_row_terms(3,1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         face_mass_flux(j)        = face_flux(1,j)
         face_mass_flux_r2(j)     = face_flux(1,j)*r_edg(j)*r_edg(j)
         momentum_largest_term(j) = max(abs(dF(2,j)), abs(S(2,j)),       &
                                        abs(Smom(j)))
         energy_largest_term(j)   = max(abs(dF(3,j)), abs(S(3,j)),       &
                                        abs(heat(j)), abs(cool(j)),      &
                                        abs(Sene(j)))
      enddo
      state_of_row_terms = u
      end subroutine store_row_terms

      ! ------------------------------------------------------!

      subroutine refresh_row_terms(u)
      ! Make the stored row terms the ones of the state u. Normally a no-op:
      ! every caller of a row scale has just had assemble_residual evaluate
      ! this state, which stores them. A caller that has not gets the Riemann
      ! solve done here instead, so that a row's scale is a function of the
      ! state alone and never of what the caller happened to evaluate before.
      !
      ! WHAT THIS PATH CANNOT REBUILD, AND WHAT THAT MEANS. heat and cool come
      ! from the ionization sweep, not from u, so they are left at the values
      ! of the state the store last described and only the hydrodynamic terms
      ! are recomputed. The states that reach here are the solver's TRIAL and
      ! FINITE-DIFFERENCE PROBES -- measured: 15 of order 1e5 scale requests on
      ! the molecular hot Uranus and 45 on WASP-121b, none of them a state the
      ! run adopts -- so the energy row of a probe can be scaled by a
      ! neighboring state's radiative terms. That is acceptable for a
      ! rejected probe and would not be for an accepted state, which never
      ! takes this path because the solve evaluates its residual first.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8, dimension(3,1-Ng:N+Ng) :: WL, WR, dF, S
      integer :: j, n_limited_before
      if (allocated(state_of_row_terms)) then
         if (all(state_of_row_terms .eq. u)) return
      endif
      ! The face-limiter ledger is restored across the rebuild: this is a
      ! measurement of a state, not a step of the run, and the count the run
      ! reports must stay the count of its own updates.
      n_limited_before = n_faces_positivity_limited
      call reconstruction_continuation_rhs(u, WL, WR, dF, S)
      n_faces_positivity_limited = n_limited_before
      if (.not. allocated(face_mass_flux_r2)) then
         ! No sweep has run yet, so there are no radiative terms to keep.
         call store_row_terms(u, dF, S, 0.0d0*dF(1,:), 0.0d0*dF(1,:),    &
                              0.0d0*dF(1,:), 0.0d0*dF(1,:))
         return
      endif
      do j = 1-Ng, N+Ng
         face_mass_flux(j)        = face_flux(1,j)
         face_mass_flux_r2(j)     = face_flux(1,j)*r_edg(j)*r_edg(j)
         momentum_largest_term(j) = max(abs(dF(2,j)), abs(S(2,j)))
         energy_largest_term(j)   = max(abs(dF(3,j)), abs(S(3,j)),       &
                                        energy_largest_term(j))
      enddo
      state_of_row_terms = u
      end subroutine refresh_row_terms

      ! ------------------------------------------------------!

      subroutine face_mass_flux_of_state(rho_state, Frho)
      ! THE FACE MASS FLUX EVERY SPECIES ROW RIDES ON.  A transported species
      ! crosses a face carrying F_rho times its face mass fraction
      ! (species_face_flux.f90), so a stationary species balance and the
      ! marching stages are the same operator only if they read the same
      ! F_rho: the one the Riemann solve of THIS state returned and the mass
      ! row of this state differences.
      !
      ! IT IS THE STORED ONE, AND THE STATE IS CHECKED.  The caller supplies
      ! the density of the state it is measuring; the row terms carry the
      ! state they were built from, so a caller that has not assembled the
      ! mass row of this state is refused here rather than silently given the
      ! flux of another one.  Every path that measures a species balance --
      ! the stationary residual, the certification of a returned state, the
      ! marching stop report -- assembles that residual first.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho_state
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Frho
      integer :: j

      if (.not. allocated(face_mass_flux)) then
         write(*,'(A)') ' ERROR: a species transport balance was asked'// &
              ' for before any mass row was assembled; the face mass'//   &
              ' flux it rides on does not exist yet'
         error stop 1
      endif
      ! THE PHYSICAL CELLS IDENTIFY THE STATE.  The ghosts are derived from
      ! them by the boundary condition and a caller may hold a copy of the
      ! density written before the last ghost fill, which is a copy of the
      ! same state; the fluxes stored here were built from the ghosts that
      ! fill left, which is the state the mass row differences.
      do j = 1, N
         if (state_of_row_terms(1,j) .ne. rho_state(j)) then
            write(*,'(A,I0)') ' ERROR: the stored face mass flux belongs'//&
                 ' to another state than the one whose species balance'//  &
                 ' was asked for; first differing cell ', j
            error stop 1
         endif
      enddo
      Frho = face_mass_flux

      end subroutine face_mass_flux_of_state

      ! ------------------------------------------------------!

      double precision function mass_flux_row_scale(j, u) result(s)
      ! THE CONTINUITY ROW'S OWN LARGEST TERM:
      !
      !   s_1(j) = max( |F_{j-1/2}| r_{j-1/2}^2, |F_{j+1/2}| r_{j+1/2}^2 )
      !            / dV_j ,   dV_j = ( r_{j+1/2}^3 - r_{j-1/2}^3 ) / 3
      !
      ! dV_j is the volume RK_rhs divides the flux difference by, so R_1/s_1
      ! is exactly the fractional change of the mass flux across the cell. The
      ! mass row has no source (Source.f90, S(1) = 0), so that fraction is the
      ! whole of what the row says.
      !
      ! WHY NOT rho (|v| + c_s)/dr, WHICH IS WHAT THIS WAS. That is the state
      ! scale times the cell's signal-crossing rate, and nothing in the
      ! continuity row moves at c_s. MEASURED (docs/p54_base_layer_mass_flux.md
      ! section 3): dividing by it gives identically the fractional flux error
      ! per cell TIMES the local Mach number, verified to one percent on a
      ! converged state at r = 1.001 to 2.1. In a quasi-hydrostatic layer at
      ! Mach 5e-5 that is a factor 2e4 of blindness -- a tolerance of 1e-5 at
      ! r = 1.01 admitted a mass flux changing by 19 percent from one cell to
      ! the next, and the 1 microbar hot-Uranus rung was accepted with a flux
      ! 30 percent out at 1.03 R_p while every residual read 3e-6.
      !
      ! The face mass flux is single signed over the whole column on every
      ! state measured, so the max cannot vanish by cancellation the way the
      ! cell-centred rho v r^2 does at the base (section 2).
      integer,                        intent(in) :: j
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8 :: rp, rm
      call refresh_row_terms(u)
      rp = r_edg(j);  rm = r_edg(j-1)
      s = max(abs(face_mass_flux_r2(j)), abs(face_mass_flux_r2(j-1)))    &
          /((rp*rp*rp - rm*rm*rm)/3.0d0)
      s = max(s, tiny(1.0d0))
      end function mass_flux_row_scale

      ! ------------------------------------------------------!

      double precision function momentum_row_scale(j, u) result(s)
      ! THE MOMENTUM ROW'S OWN LARGEST TERM:
      !
      !   s_2(j) = max( |dF_2(j)| , |S_2(j)| )
      !
      ! the momentum flux divergence the row differences and the source it is
      ! balanced against -- gravity, plus the geometric pressure term under
      ! PLM, exactly as Source.f90 discretizes them, plus the operator-split
      ! viscous source where viscosity is on.
      !
      ! This REPLACES the section 62.3 gravitational bound and the
      ! signal-speed bound both: the source term IS the gravitational bound,
      ! taken from the same expression the residual subtracts rather than
      ! rebuilt from a potential difference, so there is nothing left for the
      ! max over the two to protect against. In a quasi-hydrostatic layer the
      ! two terms nearly cancel and are enormous beside the wind's mass flux
      ! -- measured, max(|dF_2|,|S_2|)/F_0 is 4.4e5 at 1.005 R_p and 1.9e2 at
      ! 1.2 on the molecular hot Uranus -- so a row that reads 1e-8 against
      ! rho(|v|+c_s)/dr can still be displacing the mass flux, and it was:
      ! docs/p55_base_mode.md section 4 predicts the marching departure rate of
      ! an accepted state from this ratio to 9 percent.
      integer,                        intent(in) :: j
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      call refresh_row_terms(u)
      s = max(momentum_largest_term(j), tiny(1.0d0))
      end function momentum_row_scale

      ! ------------------------------------------------------!

      double precision function energy_row_scale(j, u) result(s)
      ! THE ENERGY ROW'S OWN LARGEST TERM:
      !
      !   s_3(j) = max( |dF_3(j)| , |S_3(j)| , heat(j) , cool(j) )
      !
      ! the energy flux divergence (including the gravitational work term
      ! RK_rhs folds into it), the source, and the two radiative terms the row
      ! is the difference of, plus the operator-split conduction source where
      ! conduction is on. Every one of them is already computed; the cost is
      ! the max.
      !
      ! MEASURED on the states the previous scale accepted: the energy row of
      ! the molecular hot Uranus's layer reads 1.7e-3 to 4.7e-3 of its own
      ! largest term at 1.03 to 1.05 R_p, and cell 1 reads 0.99, while the
      ! same rows read 1e-8 to 1e-4 against E(|v|+c_s)/dr
      ! (docs/p54g23_row_scale_scan.md section ii). A state can be steady on
      ! the old measure and out by one part in a thousand on the physics the
      ! row is balancing.
      integer,                        intent(in) :: j
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      call refresh_row_terms(u)
      s = max(energy_largest_term(j), tiny(1.0d0))
      end function energy_row_scale

      ! ------------------------------------------------------!

      double precision function residual_row_scale(k, j, u) result(s)
      ! The scale that makes residual row k of cell j dimensionless. This is
      ! the ONLY definition of that scale in the code; the steady solver's
      ! convergence measure (steady_newton.f90) reads it from here.
      !
      ! EVERY ROW IS DIVIDED BY THE LARGEST TERM IT ITSELF CONTAINS -- the
      ! flux it differences, the source it is balanced against, and for the
      ! energy row the two radiative terms it is the difference of. The three
      ! expressions are written out and justified at mass_flux_row_scale,
      ! momentum_row_scale and energy_row_scale above; all three are already
      ! computed by RK_rhs and the ionization sweep, so the cost is the max.
      !
      ! WHAT THIS REPLACES, AND WHY.  Until section 143 the scale was
      !     s_k(j) = D_k(j) ( |v_j| + c_s,j ) / dr_j
      ! with D_k the cell's state scale: an upper bound on how fast row k
      ! could move the cell's content if its whole inventory travelled at the
      ! signal speed. That is a bound, but on the wrong quantity -- nothing in
      ! any of these rows moves at c_s -- and in a quasi-hydrostatic layer it
      ! is enormous compared with the terms the row actually holds. Measured
      ! consequences, all on states the old measure accepted at 1e-5 or better:
      !
      !   mass       R_1/s_1 was identically the fractional flux error per
      !              cell TIMES the local Mach number (verified to one percent
      !              at ten radii), so at Mach 5e-5 the gate was blind by 2e4
      !              and passed a mass flux 30 percent out at 1.03 R_p;
      !   momentum   max(|dF_2|,|S_2|)/F_0 reaches 4.4e5 in the layer, so a row
      !              reading 1e-8 on the old scale still displaced the wind's
      !              mass flux measurably per crossing time;
      !   energy     the layer's row was out by 1.7e-3 to 4.7e-3 of its own
      !              terms, and cell 1 by 99 percent, while reading 1e-8 to
      !              1e-4 on the old scale.
      !
      ! Full measurement: docs/p54_base_layer_mass_flux.md sections 3 and 10,
      ! docs/p55_base_mode.md section 4, docs/p54g23_row_scale_scan.md.
      !
      ! THE MOMENTUM ROW'S GRAVITATIONAL BOUND IS NOT A SEPARATE CASE ANY
      ! MORE. Section 62.3 added max(..., rho |dPhi/dr|) because the old scale
      ! could not see gravity; s_2 now takes gravity from S_2, the same
      ! expression the residual subtracts, so the two are one term and not two.
      !
      ! THIS SCALES THE MEASURE, NOT THE SOLVER. cell_state_scales
      ! (steady_newton.f90) keeps the state scales D_k for the Newton system
      ! and the line-search merit. Rescaling those with these quantities was
      ! tried and measured: the mass row's version varies by a factor 9 across
      ! the first two cells and left the hot-Uranus solve with no descent
      ! direction after 179 iterations against 9 for the measure-only build
      ! (docs/p54_base_layer_mass_flux.md section 10.4). So the merit and the
      ! acceptance test are deliberately no longer one expression apart.
      !
      ! There is no region switch in any row: the same expression holds at
      ! every cell, which is what makes the wind's and the layer's numbers
      ! commensurable.
      integer,                        intent(in) :: k, j
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      select case (k)
         case (1);       s = mass_flux_row_scale(j, u)
         case (2);       s = momentum_row_scale(j, u)
         case default;   s = energy_row_scale(j, u)
      end select
      end function residual_row_scale

      ! ------------------------------------------------------!

      double precision function carrier_row_scale(j, nH, vc, cs) result(s)
      ! The scale that makes the transported-carrier row of cell j
      ! dimensionless, on the same principle residual_row_scale states: the
      ! physics that balances the row where the cell is, never the unknown
      ! itself.
      !
      ! The carrier row is a number-density balance,
      !     -r^-2 d/dr [ r^2 ( n_i v + Phi_i ) ] + P_i - L_i = 0,
      ! and its unknown n_i is the one thing that CANNOT be its scale: n(H2)
      ! falls four decades across the H2 front, so |R|/n(H2) is large
      ! wherever the carrier is scarce whether or not the state is right.
      ! That is the same pathology |rho v| has in the quasi-hydrostatic
      ! layer, and it gets the same answer -- a scale built from a quantity
      ! that does not collapse there.
      !
      ! What does not collapse is the ELEMENT the carrier belongs to, moved
      ! at the cell's own fastest signal speed:
      !     s = n_H(j) ( |v_j| + c_s,j ) / dr_j
      ! the flux divergence the cell would carry if its whole hydrogen
      ! inventory travelled at |v| + c_s. It is the direct analogue of the
      ! momentum unknown's scale rho(|v| + c_s) in cell_state_scales, it is
      ! local, and it has no fitted floor.
      !
      ! nH is the hydrogen NUCLEUS density of the cell [cm^-3], vc its gas
      ! velocity and cs its sound speed [cm/s], all in physical units; the
      ! cell width is taken in physical units too (dr_j*R0), so the result is
      ! a volumetric rate [cm^-3 s^-1] -- the units the carrier residual is
      ! assembled in.
      integer, intent(in) :: j
      real*8,  intent(in) :: nH, vc, cs
      s = nH*(abs(vc) + cs)/max(dr_j(j)*R0, 1.0d0)
      end function carrier_row_scale

      ! ------------------------------------------------------!

      subroutine relnorm_over_cells(Res, u, ja, jb, rc)
      ! THE relative residual of each row over the cells [ja:jb]:
      !
      !     rc(k) = max_{j in [ja:jb]}  |R_kj| / s_kj ,
      !
      ! every cell divided by ITS OWN row scale (residual_row_scale) and the
      ! maximum taken of those ratios. "Every cell is steady on the scale its
      ! own physics sets", which is what this norm has always been described
      ! as and, until section 145, was not.
      !
      ! WHAT THIS REPLACES, AND WHY (docs/p55_base_mode.md section 11). There
      ! were two forms and neither was this one:
      !
      !   volume weighted   sum_j |R_kj| V_j / sum_j s_kj V_j
      !   ratio of maxima   max_j |R_kj|     / max_j s_kj
      !
      ! The first lets a cell be paid for by the volume of the others; the
      ! second takes the two maxima independently, and they are not in the
      ! same cell, since s ~ 1/dr is largest in the smallest cells at the
      ! base while |R| need not be. MEASURED on the two states of section
      ! 144.5, both accepted at ||R|| ~ 1e-5 under the volume-weighted form:
      ! the hot Uranus carries a mass row at 3.10e-5 and an energy row at
      ! 2.17e-3 in cell 3 -- 218 times the accepted number -- and WASP-121b
      ! carries seventy cells above 1e-5, fifty of them in the WIND between
      ! 1.30 and 2.00 R_p. Under this form the same solves converge with no
      ! cell of any row above the tolerance, in two to four more JFNK
      ! iterations, at the same mass-loss rate to the printed digit.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: Res, u
      integer,                        intent(in)  :: ja, jb
      real*8, dimension(3),           intent(out) :: rc
      integer :: k, j
      real*8  :: q
      rc = 0.0d0
      if (jb .lt. ja) return
      do k = 1,3
         do j = ja, jb
            q = abs(Res(k,j))/max(residual_row_scale(k, j, u), 1.0d-300)
            if (q .gt. rc(k)) rc(k) = q
         enddo
      enddo
      end subroutine relnorm_over_cells

      ! ------------------------------------------------------!

      subroutine residual_norms(R, u, rc)
      ! The relative residual rate of each row over the WHOLE physical column
      ! [1:N], evaluated separately over the wind [j_min:N] and the layer
      ! below the escape radius [1:j_min-1] and combined by the LARGER of the
      ! two.
      !
      ! Why the two regions are not merged, and why the combination is the
      ! maximum: their rows are made dimensionless by different physics (the
      ! wind's momentum density above r_esc, the gravitational force density
      ! below it; residual_row_scale), so the two numbers are not addends of a
      ! common quantity. Each is already a dimensionless statement whose
      ! target is "small", and the maximum says that EVERY part of the column
      ! is steady -- no region can be averaged away by another. With the
      ! cell-wise maximum inside relnorm_over_cells the split is no longer
      ! what protects the layer, but it still keeps the two regions' numbers
      ! reportable separately, which is what the JFNK's region line uses.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: R
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(3),           intent(out) :: rc
      real*8, dimension(3) :: rc_wind, rc_layer
      call relnorm_over_cells(R, u, j_min, N,       rc_wind)
      call relnorm_over_cells(R, u, 1,     j_min-1, rc_layer)
      rc = max(rc_wind, rc_layer)
      end subroutine residual_norms

      ! ------------------------------------------------------!

      subroutine write_residual_breakdown(Res, u, heat, cool, tag)
      ! BOTH NUMBERS, EVERY TIME A RESIDUAL IS REPORTED.
      !
      ! Acceptance is the cell-wise maximum (relnorm_over_cells, section 145),
      ! because that is the statement "no cell of any row is out". It is not
      ! the statement "the column's budget is closed", which is what the
      ! volume-weighted ratio of sums says, and the two are different
      ! measurements of a state: the first is a local guarantee, the second a
      ! shell-integrated one. The external review of 2026-09-03 (section 4)
      ! asks for both to be visible before either is used alone, and for the
      ! window-by-window contributions that say WHERE a number comes from.
      !
      ! Printed, for each hydrodynamic row:
      !   max_j |R_j|/s_j, with the cell index and radius where it sits
      !   sum_j |R_j| V_j / sum_j s_j V_j, the integrated norm
      !   the numerator and denominator contributions of the four windows
      !     r < 1.03, 1.03 <= r < 1.10, 1.10 <= r < 1.20, r >= 1.20
      ! and then, for the single worst cell of the whole table, the SIGNED
      ! dimensional residual and every term its row differences.
      !
      ! The terms come from one Reconstruct + RK_rhs on u, the same pair
      ! assemble_residual and refresh_row_terms run; the radiative terms are
      ! the caller's, since they are not functions of u alone.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: Res, u
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: heat, cool
      character(len=*),               intent(in) :: tag
      real*8, dimension(3,1-Ng:N+Ng) :: WL, WR, dF, S
      real*8, parameter :: redg(3) = (/1.03d0, 1.10d0, 1.20d0/)
      character(len=9), parameter :: wname(4) = (/'r<1.03   ', '1.03-1.10',  &
                                                  '1.10-1.20', 'r>=1.20  '/)
      character(len=8), parameter :: rname(3) = (/'mass    ', 'momentum',    &
                                                  'energy  '/)
      real*8  :: cmax(3), inum(3), iden(3), wnum(4,3), wden(4,3)
      integer :: cj(3), k, j, iw, n_lim
      real*8  :: rq, vol, sc
      cmax = 0.0d0;  cj = 1;  inum = 0.0d0;  iden = 0.0d0
      wnum = 0.0d0;  wden = 0.0d0
      do k = 1,3
         do j = 1, N
            sc = residual_row_scale(k, j, u)
            rq = abs(Res(k,j))/max(sc, 1.0d-300)
            if (rq .gt. cmax(k)) then
               cmax(k) = rq;  cj(k) = j
            endif
            vol = r(j)*r(j)*dr_j(j)
            inum(k) = inum(k) + abs(Res(k,j))*vol
            iden(k) = iden(k) + sc*vol
            iw = 1
            if (r(j) .ge. redg(1)) iw = 2
            if (r(j) .ge. redg(2)) iw = 3
            if (r(j) .ge. redg(3)) iw = 4
            wnum(iw,k) = wnum(iw,k) + abs(Res(k,j))*vol
            wden(iw,k) = wden(iw,k) + sc*vol
         enddo
      enddo
      write(*,'(A)') ' [diag] residual, both norms  ('//trim(tag)//')'
      write(*,'(A)') '   row       cellwise max  at cell     r      '//      &
           'integrated  (num/den)'
      do k = 1,3
         write(*,'(A,A,A,ES12.4,I8,F10.5,ES12.4,A,ES10.3,A,ES10.3,A)')       &
              '   ', rname(k), '  ', cmax(k), cj(k), r(cj(k)),               &
              inum(k)/max(iden(k), 1.0d-300), '  (', inum(k), ' /',          &
              iden(k), ')'
      enddo
      write(*,'(A)') '   window contributions, numerator / denominator:'
      do k = 1,3
         do iw = 1,4
            write(*,'(A,A,A,A,A,ES11.3,A,ES11.3)') '     ', rname(k), ' ',   &
                 wname(iw), ' :', wnum(iw,k), ' /', wden(iw,k)
         enddo
      enddo
      ! The worst cell of the whole table, with its row's terms.
      k = 1
      if (cmax(2) .gt. cmax(k)) k = 2
      if (cmax(3) .gt. cmax(k)) k = 3
      j = cj(k)
      n_lim = n_faces_positivity_limited
      call reconstruction_continuation_rhs(u, WL, WR, dF, S)
      n_faces_positivity_limited = n_lim
      write(*,'(A,A,A,I0,A,F9.5,A)') '   worst cell: row ', rname(k),        &
           ' cell ', j, ' at r =', r(j), ' -- signed terms [code units]:'
      write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)')                                &
           '     R =', Res(k,j), '   flux divergence =', dF(k,j),            &
           '   source =', S(k,j)
      if (k .eq. 3) write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)')                  &
           '     heat =', heat(j), '   cool =', cool(j),                     &
           '   heat-cool =', heat(j) - cool(j)
      write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)')                                &
           '     scale =', residual_row_scale(k, j, u),                      &
           '   rho =', u(1,j), '   rho v =', u(2,j)
      end subroutine write_residual_breakdown

      subroutine flux_spread_of_state(u, spread, fmean)
      ! THE FLUX GATE: the radial spread of the mass flux the scheme
      ! CONSERVES, over the wind window [j_flux:N], i.e. over r >= r_flux:
      !
      !     spread = ( max_f F_f - min_f F_f ) / |mean_f F_f| ,
      !     F_f    = F_{f+1/2} r_{f+1/2}^2 , the Riemann face mass flux,
      !
      ! taken over the faces bounding those cells, excluding the outermost
      ! face of the domain, whose state is a ghost extrapolation.
      !
      ! WHY THE FACE AND NOT rho v r^2 AT THE CELL CENTRE, which this was.
      ! The finite-volume update moves the face flux; the cell-centred product
      ! is a reconstruction of the state and is not the conserved quantity, and
      ! the mass row's own residual is the fractional change of the FACE flux
      ! across a cell (mass_flux_row_scale). MEASURED on the two accepted
      ! states of docs/p55_base_mode.md section 11.3: above r = 1.10 the face
      ! flux takes ONE double-precision value over three hundred cells, while
      ! the cell-centred product read 2.5e-4 over the same window -- and the
      ! sum of the mass row's own per-cell fractions bounds the face flux's
      ! total variation at 3.2e-10 there. The gate was reporting the
      ! difference between two functionals, not a failure of conservation, and
      ! it did so by six orders.
      !
      ! It also removes the cell-1 artifact for free: the centred product there
      ! is -21 times the wind flux and changes sign (section 2 of
      ! docs/p54_base_layer_mass_flux.md), while both of that cell's faces
      ! carry inflow of order the wind flux.
      !
      ! THE SIGN IS KEPT. A window in which the flux reverses is not steady,
      ! and with the magnitude alone that would read as a small spread.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8,                         intent(out) :: spread, fmean
      call flux_spread_over_faces(u, j_flux, spread, fmean)
      end subroutine flux_spread_of_state

      ! ------------------------------------------------------------- !

      subroutine flux_spread_above_radius(u, r_min, spread, fmean)
      ! The same spread of the same quantity, over an ARBITRARY inner radius,
      ! for REPORTING ONLY. The acceptance test is flux_spread_of_state above,
      ! on the r >= r_flux window; this exists because that window cannot show
      ! how far the non-flatness reaches.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8,                         intent(in)  :: r_min
      real*8,                         intent(out) :: spread, fmean
      integer :: j, ja
      ja = N
      do j = 1, N
         if (r(j) .ge. r_min) then
            ja = j;  exit
         endif
      enddo
      call flux_spread_over_faces(u, ja, spread, fmean)
      end subroutine flux_spread_above_radius

      ! ------------------------------------------------------------- !

      subroutine flux_spread_over_faces(u, ja, spread, fmean)
      ! The one implementation both of the above call: the spread of
      ! F_{f+1/2} r^2 over the faces ja-1 .. N-1, which are the interior
      ! faces bounding cells ja..N. refresh_row_terms guarantees the cached
      ! face fluxes belong to u -- it is the same guard the row scales use, so
      ! there is no second staleness rule and no second Riemann solve.
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      integer,                        intent(in)  :: ja
      real*8,                         intent(out) :: spread, fmean
      integer :: f, nf
      real*8  :: ff, fmx, fmn, fsum
      spread = 0.0d0;  fmean = 0.0d0
      if (ja .gt. N-1) return
      call refresh_row_terms(u)
      fmx = -huge(1.0d0);  fmn = huge(1.0d0);  fsum = 0.0d0;  nf = 0
      do f = max(ja-1,1), N-1
         ff   = face_mass_flux_r2(f)
         fmx  = max(fmx, ff);  fmn = min(fmn, ff)
         fsum = fsum + ff;     nf  = nf + 1
      enddo
      if (nf .le. 1) return
      fmean  = fsum/dble(nf)
      spread = (fmx - fmn)/max(abs(fmean), 1.0d-30)
      end subroutine flux_spread_over_faces

      ! ------------------------------------------------------!

      logical function row_terms_describe_state(u) result(ok)
      ! WHETHER THE STORED ROW TERMS BELONG TO THE STATE u.
      !
      ! residual_row_scale divides by terms that store_row_terms left behind,
      ! so a caller that has not assembled the residual of THIS state would
      ! be scaling its rows by another state's physics -- and before any
      ! assembly at all the arrays do not exist. store_row_terms keeps the
      ! state the terms came from beside them precisely so that the question
      ! can be asked; a certification that cannot answer it reports the
      ! hydrodynamic rows as unavailable rather than as a number.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      ok = allocated(face_mass_flux_r2) .and. allocated(state_of_row_terms)
      if (.not. ok) return
      ok = all(state_of_row_terms .eq. u)
      end function row_terms_describe_state

      ! ------------------------------------------------------!

      pure integer function n_cells_without_chemical_root(acc_n)        &
                                                            result(n_no)
      ! THE CELLS OF ONE EQUILIBRIUM SWEEP WHOSE ACCEPTED COMPOSITION IS NOT
      ! A ROOT OF THE CHEMICAL NETWORK, and the single reading of the
      ! acceptance classes every consumer of them uses.
      !
      ! The classes are defined where they are produced (ionization_equilib-
      ! rium.f90, the acceptance block of each branch). Four of the six are
      ! ROOTS of the requested equations: the state is inside the element
      ! simplex and the largest normalized reaction residual of any row is at
      ! or below ieq_res_tol.
      !
      !   1  root, and the cell solver also reported convergence
      !   2  root without a converging solver status (MINPACK info is an
      !      xtol statement about the step, so it is neither sufficient nor
      !      necessary; the residual decides)
      !   3  a state projected onto the element budget -- the molecular
      !      clamp, or the uncoupled ionization balance handed back -- whose
      !      residual, RECHECKED after the projection, still marks a root
      !   5  root of the constrained element-conserving continuation solve,
      !      certified by the same residual test as class 1
      !
      !   4  NOT a root: the largest reaction residual is above ieq_res_tol,
      !      or it is not a finite number. Such a state is kept only by the
      !      relaxation amnesty, so that a cold start may pass through it,
      !      and the chemistry does not describe it.
      !   6  NOT a root either, and the candidate was so far from one
      !      (residual above ieq_nonroot_res_cap) that it was not adopted at
      !      all: the cell kept the composition it entered the sweep with.
      !      That composition solved an EARLIER cell state, not this one, so
      !      the chemistry does not describe this state.
      !
      ! A class-5 root therefore has the validity of a class-1 root: the
      ! continuation is how the root was FOUND, not a weaker certificate.
      ! Classes 4 and 6 are the two that count here.
      integer, intent(in) :: acc_n(6)
      n_no = acc_n(4) + acc_n(6)
      end function n_cells_without_chemical_root

      ! ------------------------------------------------------!

      logical function steady_gates_met(rnorm, u, resid_tol, fspread,   &
                                        carrier_is_unknown,             &
                                        carrier_relnorm,                &
                                        n_no_chem_root) result(ok)
      ! THE ACCEPTANCE TEST of a steady state, and the ONLY one: every
      ! route that may stop a run on "the equations are satisfied" -- the
      ! JFNK and PTC solves, and the marching loop under "Resid tol:" --
      ! calls this function, so ONE definition of "accepted" exists. It
      ! lives here, with the residual and the flux functional it combines,
      ! rather than inside either solver, because neither solver owns the
      ! definition; the carrier state the third gate needs is passed in for
      ! the same reason.
      !
      !   residual gate   ||R|| < resid_tol
      !                   nothing in the column is changing, each row measured
      !                   as the fraction of its own largest term
      !                   (residual_row_scale).
      !   flux gate       spread of rho v r^2 over r >= r_flux < flux_spread_th
      !                   the wind carries ONE mass flux at every altitude
      !                   (ATES / CETIMB), which no residual norm says: a
      !                   slowly varying flux has a small local derivative, so
      !                   the two states of section 126.6 differ by a factor
      !                   25 in the spread and by only 1.9 to 5.7 in every
      !                   residual norm that was tried (section 133).
      !   chemistry gate  no cell of the state carries a composition that is
      !                   not a root of the chemical network (acceptance
      !                   class 4, n_cells_without_chemical_root). Heat and
      !                   cool in the energy row are evaluated on the
      !                   accepted composition, so on such a cell the row is
      !                   a number computed on a state the chemistry does not
      !                   describe and it says nothing about how far the
      !                   hydrodynamics is from steady. The count is passed in
      !                   because it belongs to the equilibrium sweep of the
      !                   state being judged, not to the residual.
      !                   A CALLER THAT OMITS IT MAKES NO STATEMENT ABOUT THE
      !                   CHEMISTRY of the state, and the gate then rests on
      !                   the other three; the marching stop is such a caller
      !                   today, because the marching loop does not ask
      !                   ioniz_eq for the ledger of its sweep.
      !   carrier gate    the carrier row balances to carrier_resid_th of its
      !                   own largest terms (section 139). Applied only where
      !                   n(H2) is among the unknowns, which the caller states
      !                   through carrier_is_unknown: it is a gate of its own
      !                   and not a row of ||R|| because the two are not
      !                   calibrated against each other.
      !
      ! flux_spread_th <= 0 disables the flux gate and leaves the residual
      ! alone in charge, which is what a run gets by asking for it explicitly.
      !
      ! WHY THE MARCHING LOOP PASSES carrier_is_unknown = .false. Its carriers
      ! are moved by the operator-split transport step, not solved for, so no
      ! carrier residual OF A SOLVE exists for that state; gating its stop on
      ! the last solve's number would be reporting a measurement the state
      ! never had. A coupled-carrier run reaches its answer through the JFNK
      ! finish, which passes the measurement it did make.
      real*8,                         intent(in)  :: rnorm, resid_tol
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8,                         intent(out) :: fspread
      logical,                        intent(in)  :: carrier_is_unknown
      real*8,                         intent(in)  :: carrier_relnorm
      integer, optional,              intent(in)  :: n_no_chem_root
      real*8 :: fmean
      call flux_spread_of_state(u, fspread, fmean)
      ok = (rnorm .lt. resid_tol)
      if (flux_spread_th .gt. 0.0d0) ok = ok .and.                        &
                                          (fspread .lt. flux_spread_th)
      if (carrier_is_unknown) ok = ok .and.                               &
                                   (carrier_relnorm .lt. carrier_resid_th)
      if (present(n_no_chem_root)) ok = ok .and. (n_no_chem_root .eq. 0)
      end function steady_gates_met

      ! End of module
      end module steady_residual_mod

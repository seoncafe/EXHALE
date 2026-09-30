      module steady_residual_mod
      ! Finite-volume STEADY residual  R = du/dt  (zero at a true steady
      ! state), shared by the EXHALE_RESIDUAL diagnostic, the in-loop
      ! convergence monitor, and the steady-state Newton/PTC solver.
      !
      !   R(:,1) = dF - S                          (mass)
      !   R(:,2) = dF - S - F_mu                   (momentum)
      !   R(:,3) = dF_E - S_E - (heat-cool)
      !                       - (w F_mu + q_mu + conduction)
      !                       + div q_d                        (energy)
      !
      ! q_d is the interdiffusion enthalpy flux of a mixture whose elements
      ! move relative to one another (He_diffusion; Cook 2009, Phys. Fluids
      ! 21, 055109, eqs. 11-13), formed from the composition of the state by
      ! interdiffusion_enthalpy_divergence_of_state (binary_element_
      ! diffusion) and identically zero, with no arithmetic, when the
      ! elements do not move.  The composition is therefore an argument of
      ! every assembly: the energy row of a diffusing mixture is not a
      ! function of u, n_part, heat and cool alone.
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
      use grid_construction, only: spherical_cell_volume
      use Conversion, only: U_to_W
      use caloric_eos, only: pressure_from_energy_density,          &
                             adiabatic_index_from_state
      use Reconstruction_step
      use RK_integration
      use viscous_conduction, only: transport_active,                   &
                                    viscous_conduction_sources
      use hydrodynamic_rows, only: generic_precision_rows_selected,     &
                                   ROWS_PRODUCTION, ROWS_QUADRUPLE,     &
                                   ROWS_GENERIC_DOUBLE,                 &
                                   hydrodynamic_rows_in_quadruple_precision, &
                                   hydrodynamic_rows_in_double_precision
      use ionization_equilibrium, only: ieq_sweep_state_kind,           &
                                        ieq_state_marching
      use conservation_budget, only: conservation_budget_exports_left,  &
                                     write_conservation_budget_terms
      use binary_element_diffusion, only: interdiffusion_enthalpy_active, &
                           interdiffusion_enthalpy_divergence_of_state

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
                face_mass_flux_of_state,                             &
                mass_row_rounding_floor,                             &
                energy_row_divergence_of_state,                      &
                carrier_enthalpy_divergence

      ! THE ENTHALPY THE TRANSPORTED CARRIERS CARRY (code audit of
      ! 2026-09-29, F1). The divergence is formed by the carrier operator
      ! (diffusive_photochemistry, carrier_enthalpy_divergence_of_state),
      ! which uses this module and so cannot be used by it; it registers
      ! the procedure here when the carrier set is decided
      ! (carrier_set_init), and a run without carriers leaves the pointer
      ! null and the term absent.
      abstract interface
         subroutine energy_row_divergence_of_state(rho, Tcode, f_sp,    &
                                                   divq, q_out)
            import :: N, Ng, n_species
            real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho
            real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Tcode
            real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
            real*8, dimension(1-Ng:N+Ng),           intent(out) :: divq
            real*8, dimension(0:N), optional,       intent(out) :: q_out
         end subroutine energy_row_divergence_of_state
      end interface
      procedure(energy_row_divergence_of_state), pointer, save ::      &
                carrier_enthalpy_divergence => null()

      ! THE TERMS EACH CONSERVATION ROW IS BUILT FROM, as assemble_residual
      ! last produced them, together with the state they belong to. The three
      ! row scales below are these; nothing else reads them.
      !
      !   face_mass_flux_r2(j)     F_{j+1/2} r_{j+1/2}^2     (mass, SIGNED:
      !                            the flux gate needs the sign, so a flow
      !                            that reverses inside the window is a huge
      !                            spread and not a small one; the row scale
      !                            takes the magnitude where it needs one)
      !   momentum_largest_term(j) max(|ram|, |dp/dr|, |rho dphi/dr|)
      !                            (momentum)
      !   energy_largest_term(j)   max(|dF_3|, |S_3|, heat, cool, |Sene|,
      !                            |Sidf|)  (energy: Sene the viscous and
      !                            conduction sources where those are
      !                            active, Sidf the divergence of the
      !                            enthalpy flux of the element and
      !                            carrier fluxes)
      !
      ! THE MOMENTUM ROW'S THREE TERMS ARE THE TERMS OF THE EQUATION, and
      ! they are not dF_2 and S_2. Under PLM the pressure sits partly in
      ! the momentum flux and partly in the geometric source, each of them
      ! carrying an O(2 p/r) part that cancels against the other, so
      ! max(|dF_2|, |S_2|) reads 2 p/r where the physical force is zero;
      ! under the well-balanced option the whole equilibrium pressure force
      ! and the weight cancel in the algebra and S_2 is zero, so a row
      ! scaled by its remaining terms alone would be its own scale and read
      ! one however small the imbalance became. RK_rhs gathers the pieces
      ! back into the ram divergence, the spherical pressure gradient and
      ! the weight (momentum_ram_divergence, momentum_pressure_gradient,
      ! momentum_gravity), and the weight under that option is the equilibrium
      ! pressure force, which therefore enters the max once.
      !
      ! The state is kept so that a caller asking for the scale of a
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

      ! Whether the continuity row's rounding-floor measurement is inside
      ! its own perturbed assembly, so the hook does not re-enter, and
      ! whether it has already reported: it is a property of a state and
      ! one state says it.
      logical, save :: mass_floor_scan_running = .false.
      logical, save :: mass_floor_scan_reported = .false.

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
      real*8, dimension(1-Ng:N+Ng)   :: fp_plm, epf_plm
      real*8, dimension(1-Ng:N+Ng)   :: ram_plm, pgr_plm, grv_plm
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
         ff_plm  = face_flux
         fp_plm  = face_p
         ram_plm = momentum_ram_divergence
         pgr_plm = momentum_pressure_gradient
         grv_plm = momentum_gravity
         if (well_balanced) epf_plm = equilibrium_pressure_force

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
         ! The momentum row on the homotopy is the same combination of the
         ! two schemes' rows, so each term of the equation it holds, and
         ! the equilibrium pressure force it is read against under that option,
         ! is the same combination of the two forms.
         momentum_ram_divergence    = om*ram_plm                         &
                                    + lam*momentum_ram_divergence
         momentum_pressure_gradient = om*pgr_plm                         &
                                    + lam*momentum_pressure_gradient
         momentum_gravity           = om*grv_plm                         &
                                    + lam*momentum_gravity
         if (well_balanced)                                              &
            equilibrium_pressure_force = om*epf_plm                      &
                                       + lam*equilibrium_pressure_force
      endif

      rec_method = sav_method
      use_plm    = sav_plm
      use_weno3  = sav_weno

      end subroutine reconstruction_continuation_rhs

      ! ------------------------------------------------------!

      ! ------------------------------------------------------!

      subroutine assemble_residual(u, n_part, f_sp, heat, cool, R)
      real*8, dimension(3,1-Ng:N+Ng), intent(in)  :: u
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: n_part
      ! The composition of the state u, the one n_part was formed from.
      ! The energy row reads it for the interdiffusion enthalpy flux.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: heat, cool
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: R
      ! Local scratch so callers' own WL/WR/dF/S are untouched
      real*8, dimension(3,1-Ng:N+Ng) :: WL, WR, dF, S, W
      real*8, dimension(1-Ng:N+Ng)   :: Tc, Smom, Sene
      ! The divergence of the interdiffusion enthalpy flux, in the units of
      ! the energy row (zero unless the elements move).
      real*8, dimension(1-Ng:N+Ng)   :: Sidf, Scar
      integer :: rows_kind
      ! The continuity row's scale of every cell, formed only when the
      ! conservation budget exports this assembly, and the two flux-branch
      ! counters held across that formation.
      real*8, dimension(1-Ng:N+Ng)   :: continuity_scale
      integer :: jcell_scale, n_roe_hlle_saved, n_llf_saved
      ! The reconstruction in force at entry, restored after the terms of a
      ! kind-generic row have been evaluated under the scheme that row was
      ! actually assembled with.
      logical :: sav_plm, sav_weno
      character(len=:), allocatable :: sav_method

      ! THE ARITHMETIC THE HYDRODYNAMIC ROWS ARE ASSEMBLED IN.  Normally the
      ! production routines, in double.  EXHALE_RESID_QUAD=1 sends the
      ! STATIONARY evaluations -- and only those; the marching stages reach
      ! RK_rhs directly and never come here -- through the
      ! quadruple-precision instantiation of the same kind-generic text, as
      ! the control experiment of where the residual's non-smoothness floor
      ! comes from: that floor is the rounding of this flux assembly, and an
      ! assembly that lowers the rounding by eighteen decades and nothing else
      ! separates the rounding from every other candidate. Default off, and
      ! nothing is adopted from it; see the header of module
      ! hydrodynamic_rows.
      rows_kind = ROWS_PRODUCTION
      if (ieq_sweep_state_kind .ne. ieq_state_marching)                  &
         rows_kind = generic_precision_rows_selected()
      if (rows_kind .eq. ROWS_QUADRUPLE) then
         call hydrodynamic_rows_in_quadruple_precision(u, WL, WR, dF, S)
      else if (rows_kind .eq. ROWS_GENERIC_DOUBLE) then
         call hydrodynamic_rows_in_double_precision(u, WL, WR, dF, S)
      else
         call reconstruction_continuation_rhs(u, WL, WR, dF, S)
      endif
      ! The kind-generic rows return the momentum row and store the
      ! interface fluxes, the face pressures and, under the well-balanced
      ! option, the face departures they were built from
      ! (store_the_interface_fluxes), but neither the three terms of the
      ! momentum equation nor the equilibrium pressure force the row is
      ! read against.  Both are functions of the state, the potential, the
      ! reconstruction and those stored face quantities, so the one
      ! expression of each is evaluated here for the state just assembled
      ! rather than left at whatever the last call to RK_rhs produced.  A
      ! reconstruction continuation (0 < recon_lambda < 1) never reaches
      ! the kind-generic rows: the generic text refuses it, because a blended
      ! flux carries neither scheme's pressure convention alone.  Those rows
      ! are the rounding control experiment of module hydrodynamic_rows and
      ! are default off.
      if (rows_kind .ne. ROWS_PRODUCTION) then
         ! UNDER THE SCHEME THE ROW WAS ASSEMBLED WITH, which is not always
         ! the flag pair: the kind-generic text selects the endpoint from
         ! lambda and leaves the flags alone, while the two routines below
         ! read the flags, so a run whose flags named the other endpoint
         ! scaled a row of one scheme by the terms of the other.  The
         ! production path has no such gap, because RK_rhs forms both inside
         ! the scope where reconstruction_continuation_rhs has set the flags.
         ! One rule answers it, assembled_reconstruction_is_plm, and the
         ! generic text selects on the same one.
         sav_plm    = use_plm
         sav_weno   = use_weno3
         sav_method = rec_method
         if (assembled_reconstruction_is_plm()) then
            rec_method = 'PLM';    use_plm = .true.;  use_weno3 = .false.
         else
            rec_method = 'WENO3';  use_plm = .false.; use_weno3 = .true.
         endif
         if (well_balanced) call equilibrium_pressure_force_of_state(u)
         call momentum_row_terms_of_state(WL, WR, S)
         rec_method = sav_method
         use_plm    = sav_plm
         use_weno3  = sav_weno
      endif
      R(1,:) = dF(1,:) - S(1,:)
      R(2,:) = dF(2,:) - S(2,:)
      R(3,:) = dF(3,:) - S(3,:) - (heat - cool)

      Smom = 0.0d0;  Sene = 0.0d0
      if (transport_active()) then
         call U_to_W(u, W)
         Tc = W(3,:)/n_part
         call viscous_conduction_sources(W(2,:), Tc, f_sp, Smom, Sene)
         R(2,:) = R(2,:) - Smom
         R(3,:) = R(3,:) - Sene
      endif
      ! THE ENTHALPY THE ELEMENT FLUXES CARRY.  With He_diffusion on, helium
      ! and hydrogen cross every face with different velocities, and the
      ! energy row gains the divergence of q_d = sum_s h_s J_s, positive
      ! where the moving particles take more enthalpy out of the cell than
      ! they bring in.  The element fluxes, the faces and the geometry are
      ! those of the element rows of the same state, so the energy this
      ! term moves telescopes over the column as the element flux does.
      ! The temperature is the one this assembly's own conduction term and
      ! the marching loop use, p/n_part.
      Sidf = 0.0d0
      if (interdiffusion_enthalpy_active()) then
         call U_to_W(u, W)
         Tc = W(3,:)/n_part
         call interdiffusion_enthalpy_divergence_of_state(W(1,:), Tc,     &
                                                          f_sp, Sidf)
         R(3,:) = R(3,:) + Sidf
      endif
      ! THE ENTHALPY THE CARRIER FLUXES CARRY, the same physics for the
      ! relative motion of a carrier within its own elements (diffusive_
      ! photochemistry, carrier_enthalpy_face_flux states it).  Added into
      ! Sidf, the one record of the enthalpy of relative transport that the
      ! row terms and the conservation budget report.  Its column sum is
      ! zero: both ends of the carrier column are closed.
      if (associated(carrier_enthalpy_divergence)) then
         call U_to_W(u, W)
         Tc = W(3,:)/n_part
         call carrier_enthalpy_divergence(W(1,:), Tc, f_sp, Scar)
         R(3,:) = R(3,:) + Scar
         Sidf   = Sidf + Scar
      endif
      ! THE LOWEST GHOST CELL CARRIES NO EQUATION: it has no lower face, so
      ! there is no balance to state there and its row is zero (RK_rhs
      ! defines the same column of dF and S as zero for the same reason).
      ! Written explicitly so that no reader of R can find a heating rate
      ! standing alone in a row that has no fluxes.
      R(:,1-Ng) = 0.0d0
      call store_row_terms(u, dF, S, heat, cool, Smom, Sene, Sidf)

      ! THE COMPLETE TERMS OF THE THREE ROWS OF THIS ASSEMBLY, for an
      ! independent reader of the discrete balance. Default off: nothing is
      ! written and no file is opened unless EXHALE_CONSERVATION_BUDGET is
      ! set (module conservation_budget). It writes the arrays above, the
      ! face states WL, WR the fluxes were formed from and the face data
      ! this assembly stored, and changes none of them.
      !
      ! The continuity row's scale is the one the stationary solve divides
      ! by, residual_row_scale, read here because this module owns it. The
      ! row terms of u were stored just above, so it forms no Riemann
      ! problem; the two flux-branch counters are restored all the same, so
      ! that no path of the export can move a count the run reports.
      if (conservation_budget_exports_left()) then
         n_roe_hlle_saved = n_faces_roe_hlle
         n_llf_saved      = n_faces_llf
         continuity_scale(1-Ng) = 0.0d0
         do jcell_scale = 2-Ng, N+Ng
            continuity_scale(jcell_scale) =                             &
               residual_row_scale(1, jcell_scale, u)
         enddo
         n_faces_roe_hlle = n_roe_hlle_saved
         n_faces_llf      = n_llf_saved
         call write_conservation_budget_terms(u, WL, WR, dF, S, R,      &
                                              heat, cool, Smom, Sene,   &
                                              Sidf, continuity_scale,   &
                                              rows_kind)
      endif

      ! The rounding floor of the continuity row, measured on the first
      ! state a run assembles a stationary residual for when
      ! EXHALE_MASS_FLOOR_SCAN arms it. It restores the terms stored above
      ! and decides nothing.
      if (.not. mass_floor_scan_running .and.                          &
          .not. mass_floor_scan_reported) then
         if (mass_floor_scan_armed()) then
            mass_floor_scan_reported = .true.
            call mass_row_rounding_floor_scan(u, n_part, f_sp, heat,    &
                                              cool, R)
         endif
      endif

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

      subroutine store_row_terms(u, dF, S, heat, cool, Smom, Sene, Sidf)
      ! Keep the terms of each row as RK_rhs and the ionization sweep have
      ! just produced them, and the state they belong to.  Sidf is the
      ! divergence of the interdiffusion enthalpy flux, a term of the energy
      ! row like the others.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u, dF, S
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: heat, cool, Smom, Sene
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: Sidf
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
         ! THE LARGEST TERM OF THE MOMENTUM EQUATION, not of its
         ! discretization: RK_rhs gathers the ram divergence, the whole
         ! spherical pressure gradient and the weight out of dF(2) and
         ! S(2), which split the pressure between them under PLM and
         ! cancel it against the weight under the well-balanced option.  The
         ! weight is momentum_gravity either way, so the equilibrium
         ! pressure force enters once and not twice.
         momentum_largest_term(j) =                                      &
            max(abs(momentum_ram_divergence(j)),                         &
                abs(momentum_pressure_gradient(j)),                      &
                abs(momentum_gravity(j)),                                &
                abs(Smom(j)))
         energy_largest_term(j)   = max(abs(dF(3,j)), abs(S(3,j)),       &
                                        abs(heat(j)), abs(cool(j)),      &
                                        abs(Sene(j)), abs(Sidf(j)))
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
                              0.0d0*dF(1,:), 0.0d0*dF(1,:),          &
                              0.0d0*dF(1,:))
         return
      endif
      do j = 1-Ng, N+Ng
         face_mass_flux(j)        = face_flux(1,j)
         face_mass_flux_r2(j)     = face_flux(1,j)*r_edg(j)*r_edg(j)
         momentum_largest_term(j) =                                      &
            max(abs(momentum_ram_divergence(j)),                         &
                abs(momentum_pressure_gradient(j)),                      &
                abs(momentum_gravity(j)))
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
      ! continuity row moves at c_s. MEASURED: dividing by it gives
      ! identically the fractional flux error
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
      call refresh_row_terms(u)
      s = max(abs(face_mass_flux_r2(j)), abs(face_mass_flux_r2(j-1)))    &
          /spherical_cell_volume(j)
      s = max(s, tiny(1.0d0))
      end function mass_flux_row_scale

      ! ------------------------------------------------------!

      double precision function mass_row_rounding_floor(j, u) result(f)
      ! THE ROUNDING FLOOR OF THE CONTINUITY ROW OF CELL j, in the units the
      ! row is judged in: |R_1|/s_1, the fractional change of the mass flux
      ! across the cell.
      !
      ! The row is a difference of two face fluxes over the cell volume,
      !
      !     R_1(j) = ( A_+ F_+  -  A_- F_- ) / dV_j ,
      !
      ! and each interface flux is assembled from the two reconstructed
      ! states and from their JUMPS. Where the column is nearly hydrostatic
      ! the jumps stand many decades below the states themselves, so a
      ! correctly rounded interface flux carries the last bit of the
      ! O(rho c_s) quantities the Riemann solve is built from and not the
      ! last bit of its own value rho v (MEASURED, N33: the first quantity
      ! of the assembly that steps is the interface flux, by 4.8e-16 of an
      ! O(1) state, where its own value is smaller by 3.2e5 to 5.9e6). The
      ! smallest flux difference the arithmetic can resolve at that face is
      ! therefore
      !
      !     eps * rho ( |v| + c_s ) * A / dV_j ,
      !
      ! the momentum the cell would carry at its own fastest signal speed,
      ! and dividing by the row's own scale s_1 = max_faces|A F| / dV_j
      ! leaves
      !
      !     f(j) = eps * max_faces[ rho (|v| + c_s) A ] / max_faces| A F | ,
      !
      ! which is one ulp in the wind, where the flux IS the momentum of the
      ! cell, and about eps/Mach in a layer whose flux is a tiny residue of
      ! it. Nothing in the continuity row moves at c_s -- that is why the
      ! row's SCALE is the flux and not this quantity (mass_flux_row_scale)
      ! -- but the rounding of the row does, because the cancellation that
      ! forms the flux happens between quantities of that size.
      !
      ! rho(|v| + c_s) is the momentum unknown's own state scale
      ! (state_scales_of_cell), taken at the two cells the faces separate.
      integer,                        intent(in) :: j
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u
      real*8 :: d1, d3, vv, cs, dm, dp, ap, am, sc
      call state_scales_of_cell(j-1, u(1,j-1), u(2,j-1), u(3,j-1),        &
                                d1, dm, d3, vv, cs)
      call state_scales_of_cell(j,   u(1,j),   u(2,j),   u(3,j),          &
                                d1, dp, d3, vv, cs)
      ap = r_edg(j)*r_edg(j);  am = r_edg(j-1)*r_edg(j-1)
      sc = mass_flux_row_scale(j, u)
      f  = epsilon(1.0d0)*max(dm*am, dp*ap)                              &
           /spherical_cell_volume(j)/max(sc, tiny(1.0d0))
      end function mass_row_rounding_floor

      ! ------------------------------------------------------!

      logical function mass_floor_scan_armed() result(on)
      ! EXHALE_MASS_FLOOR_SCAN=1 arms the measurement of the continuity
      ! row's rounding floor. Default off, and it changes nothing a run
      ! decides: it prints one block and restores the terms it found.
      character(len=32) :: env
      call get_environment_variable('EXHALE_MASS_FLOOR_SCAN', env)
      on = (trim(env) .eq. '1')
      end function mass_floor_scan_armed

      ! ------------------------------------------------------!

      subroutine mass_row_rounding_floor_scan(u, n_part, f_sp, heat, cool, &
                                              Rgiven)
      ! WHAT THE CONTINUITY ROW OF THIS STATE CANNOT GO BELOW, MEASURED.
      !
      ! One ulp is added to the density of every physical cell and the whole
      ! flux assembly is run again. A smooth response to that perturbation
      ! is one ulp of the row's own scale, because the mass flux is linear
      ! in rho; what is measured instead is the STEP the assembly takes,
      ! which is the rounding of the interface fluxes amplified by the
      ! near-hydrostatic cancellation (N33). The step is reported against
      ! two estimates of it:
      !
      !   signal  mass_row_rounding_floor above, eps rho(|v|+c_s) A / dV
      !           over the row's scale;
      !   flux    the textbook rounding of a difference of two correctly
      !           rounded addends, eps (|A_+ F_+| + |A_- F_-|) / dV over the
      !           same scale, which is between one and two ulps everywhere
      !           because the scale is the larger addend. It is reported so
      !           that the amplification is a measured ratio and not an
      !           assumption.
      !
      ! The state is put back: the perturbed assembly leaves the module's
      ! row terms describing the perturbed state, so the state given is
      ! assembled once more at the end and its rows are asserted to
      ! reproduce bitwise.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u, Rgiven
      real*8, dimension(1-Ng:N+Ng),   intent(in) :: n_part, heat, cool
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(3,1-Ng:N+Ng) :: up, Rp, Rb
      real*8, dimension(1:N) :: sc, fsg, ffl, jmp
      real*8 :: dV, mach, d1, d2, d3, vv, cs, worst_sg, worst_fl, back
      integer :: j, jw_sg, jw_fl, nshow
      ! The scales and both estimates, on the state as given: the stored
      ! terms describe it here and will not after the perturbed assembly.
      do j = 1, N
         sc(j)  = mass_flux_row_scale(j, u)
         fsg(j) = mass_row_rounding_floor(j, u)
         dV     = spherical_cell_volume(j)
         ffl(j) = epsilon(1.0d0)                                          &
                  *(abs(face_mass_flux_r2(j)) + abs(face_mass_flux_r2(j-1)))&
                  /dV/max(sc(j), tiny(1.0d0))
      enddo
      up = u
      do j = 1, N
         up(1,j) = nearest(u(1,j), 1.0d0)
      enddo
      mass_floor_scan_running = .true.
      call assemble_residual(up, n_part, f_sp, heat, cool, Rp)
      call assemble_residual(u,  n_part, f_sp, heat, cool, Rb)
      mass_floor_scan_running = .false.
      back = 0.0d0
      do j = 1, N
         back   = max(back, abs(Rb(1,j) - Rgiven(1,j)))
         jmp(j) = abs(Rp(1,j) - Rgiven(1,j))/max(sc(j), tiny(1.0d0))
      enddo
      worst_sg = 0.0d0;  worst_fl = 0.0d0;  jw_sg = 0;  jw_fl = 0
      do j = 1, N
         if (jmp(j)/max(fsg(j), tiny(1.0d0)) .gt. worst_sg) then
            worst_sg = jmp(j)/max(fsg(j), tiny(1.0d0));  jw_sg = j
         endif
         if (jmp(j)/max(ffl(j), tiny(1.0d0)) .gt. worst_fl) then
            worst_fl = jmp(j)/max(ffl(j), tiny(1.0d0));  jw_fl = j
         endif
      enddo
      write(*,'(A)') ' '
      write(*,'(A)') ' [mass floor scan] the continuity row under one ulp'&
           //' of every cell density'
      write(*,'(A,ES9.2)') '   the given state reproduces its own mass'// &
           ' rows to ', back
      write(*,'(A)') '   cell      r        |R_1|/s_1      step/s_1'//    &
           '     floor signal   floor flux   step/signal  step/flux'//    &
           '   Mach'
      nshow = min(N, 30)
      do j = 1, nshow
         call state_scales_of_cell(j, u(1,j), u(2,j), u(3,j),             &
                                   d1, d2, d3, vv, cs)
         mach = abs(vv)/max(cs, tiny(1.0d0))
         write(*,'(A,I5,F10.5,6ES14.4,ES11.3)') '   ', j, r(j),           &
              abs(Rgiven(1,j))/max(sc(j), tiny(1.0d0)), jmp(j), fsg(j),        &
              ffl(j), jmp(j)/max(fsg(j), tiny(1.0d0)),                    &
              jmp(j)/max(ffl(j), tiny(1.0d0)), mach
      enddo
      write(*,'(A,ES11.4,A,I0,A,ES11.4,A,ES11.4,A)')                      &
           '   largest step over the signal estimate ', worst_sg,         &
           ' at cell ', jw_sg, ' (step ', jmp(max(jw_sg,1)),              &
           ', estimate ', fsg(max(jw_sg,1)), ')'
      write(*,'(A,ES11.4,A,I0,A,ES11.4,A,ES11.4,A)')                      &
           '   largest step over the flux estimate   ', worst_fl,         &
           ' at cell ', jw_fl, ' (step ', jmp(max(jw_fl,1)),              &
           ', estimate ', ffl(max(jw_fl,1)), ')'
      flush(6)
      end subroutine mass_row_rounding_floor_scan

      ! ------------------------------------------------------!

      double precision function momentum_row_scale(j, u) result(s)
      ! THE LARGEST PHYSICAL TERM OF THE MOMENTUM EQUATION IN THIS CELL:
      !
      !   s_2(j) = max( |ram| , |dp/dr| , |rho dphi/dr| , |S_visc| )
      !
      ! the three left-hand terms of
      !
      !   d(rho v)/dt + div(rho v v) + dp/dr + rho dphi/dr = S_visc ,
      !
      ! each as the evaluation that produced the row assembled it
      ! (momentum_ram_divergence, momentum_pressure_gradient,
      ! momentum_gravity in RK_rhs), plus the operator-split viscous source
      ! where viscosity is on.
      !
      ! WHY NOT max(|dF_2|, |S_2|), WHICH IS WHAT THIS WAS.  Those are the
      ! pieces the discretization splits the equation into, and the split
      ! does not follow the terms:
      !
      !   PLM    Phys_flux gives the momentum flux the pressure, so dF_2 is
      !          the ram divergence plus (A+ p_up - A- p_dn)/dV, and S_2 is
      !          the weight MINUS the geometric term (A+ - A-) p_c/dV.  The
      !          two pressure pieces are one term, the spherical pressure
      !          gradient, and each holds an O(2 p/r) part that cancels
      !          against the other.  MEASURED (grid_and_gates
      !          hydrostatic_residual, zero-gravity uniform-pressure state,
      !          where the physical force is zero): the old scale read
      !          2 p/r, r s_2/p = 2.000 for both Riemann solvers at every N,
      !          and stood 0.944 above what the row carries; the WENO3 scale
      !          of the same state agreed with the row to 2e-13.
      !   WENO3  dF_2 is the ram divergence plus (p_R - p_L)/dr in ONE
      !          number, so a near-hydrostatic cell reads |dF_2| ~ |dp/dr|
      !          ~ |rho g| and the old scale was right by the accident of
      !          that balance; where ram and pressure gradient cancel each
      !          other, as they do at a sonic point, it is smaller than
      !          either term.
      !   WB     Under "Well balanced:" the equilibrium pressure force and
      !          the weight cancel in the algebra before the row is formed,
      !          so S_2 is zero and dF_2 holds the DEPARTURE alone: the
      !          remaining terms ARE the numerator and a row divided by them
      !          reads one whatever the imbalance is.  MEASURED before the
      !          weight was put back: every normalized momentum residual of
      !          the carrier reload was exactly 1.000000E+00, against 1.954
      !          for the same state without it (Update_EXHALE N37).
      !
      ! Under that option the weight is the pressure force of the cell's own
      ! hydrostatic equilibrium,
      !
      !   |rho_j [A+ (phi_i(j) - phi_c(j))
      !         + A- (phi_c(j) - phi_i(j-1))]|/dV       (PLM form)
      !   |rho_j (phi_i(j) - phi_i(j-1))|/dr            (WENO3 form)
      !
      ! which is the same physics in that option's own discretization and is
      ! what momentum_gravity carries there; the pressure gradient is then
      ! the departure the row holds.  The numerator is untouched either way
      ! and is never re-formed from a cancelled term.
      !
      ! WHY A TERM OF THE EQUATION AND NOT A SIGNAL-SPEED BOUND.  In a
      ! quasi-hydrostatic layer the pressure gradient and the weight nearly
      ! cancel and are enormous beside the wind's mass flux -- measured,
      ! the scale over F_0 is 4.4e5 at 1.005 R_p and 1.9e2 at 1.2 on the
      ! molecular hot Uranus -- so a row that reads 1e-8 against
      ! rho(|v|+c_s)/dr can still be displacing the mass flux, and it was:
      ! the marching departure rate of an accepted state is predicted from
      ! this ratio to 9 percent.  With no gravity
      ! the weight is zero and the max falls back on the dynamic terms; on a
      ! state at rest with a uniform pressure all three terms vanish and the
      ! scale is the tiny floor below, which is the only place a fully zero
      ! row is divided by anything but a term of its own.
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
      ! same rows read 1e-8 to 1e-4 against E(|v|+c_s)/dr. A state can be
      ! steady on
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
      !   momentum   the largest term of the equation over F_0 reaches 4.4e5
      !              in the layer, so a row reading 1e-8 on the old scale
      !              still displaced the wind's mass flux measurably in a
      !              crossing time;
      !   energy     the layer's row was out by 1.7e-3 to 4.7e-3 of its own
      !              terms, and cell 1 by 99 percent, while reading 1e-8 to
      !              1e-4 on the old scale.
      !
      ! THE MOMENTUM ROW'S GRAVITATIONAL BOUND IS NOT A SEPARATE CASE ANY
      ! MORE. Section 62.3 added max(..., rho |dPhi/dr|) because the old
      ! scale could not see gravity; the weight is now one of the three
      ! terms s_2 is the max over, taken from the discretization the row was
      ! assembled by, so the two are one term and not two.
      !
      ! THIS SCALES THE MEASURE, NOT THE SOLVER. cell_state_scales
      ! (steady_newton.f90) keeps the state scales D_k for the Newton system
      ! and the line-search merit. Rescaling those with these quantities was
      ! tried and measured: the mass row's version varies by a factor 9 across
      ! the first two cells and left the hot-Uranus solve with no descent
      ! direction after 179 iterations against 9 for the measure-only build.
      ! So the merit and the
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
      ! WHAT THIS REPLACES, AND WHY. There
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
      ! Every row of every cell is divided by the same expression, its own
      ! largest term (residual_row_scale, which has no region switch), and
      ! relnorm_over_cells takes the cell-wise maximum, so the maximum of the
      ! two regions is the maximum over the whole column: the split changes
      ! no number. It is kept because it makes the two regions' numbers
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
      ! states: above r = 1.10 the face
      ! flux takes ONE double-precision value over three hundred cells, while
      ! the cell-centred product read 2.5e-4 over the same window -- and the
      ! sum of the mass row's own fractions bounds the face flux's
      ! total variation at 3.2e-10 there. The gate was reporting the
      ! difference between two functionals, not a failure of conservation, and
      ! it did so by six orders.
      !
      ! It also removes the cell-1 artifact for free: the centred product there
      ! is -21 times the wind flux and changes sign, while both of that
      ! cell's faces
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

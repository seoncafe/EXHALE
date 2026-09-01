      module diffusive_photochemistry
      ! Vertical transport of the molecular carriers, solved together with
      ! the chemistry that makes and destroys them:
      !
      !   d n_i/dt + (1/r^2) d/dr [ r^2 ( n_i v + Phi_i ) ] = P_i - L_i
      !
      !   Phi_i = - n_tot (D_i + K_zz) d f_i/dr
      !           - n_i D_i [ 1/H_i - 1/H_atm ]        (f_i = n_i/n_tot)
      !
      ! for i = H2, OH, H2O, CO.  This is milestone M3 of
      ! docs/a2_oxygen_option_design.md; sections 3 and 4 of that document
      ! are the specification and this header states what was built and what
      ! was measured, not what was intended.
      !
      ! ---------------------------------------------------------------
      ! 1. WHY THE CARRIERS CANNOT BE A LOCAL STEADY STATE
      !
      ! Everywhere else in EXHALE the composition of a cell is the root of a
      ! local algebraic system: f_sp is not advected by the hydro update, it
      ! is recomputed from scratch each step by ioniz_eq.  That is correct
      ! wherever the chemistry is fast compared with the flow, and it is the
      ! reason the code can carry thirty ion stages cheaply.
      !
      ! The measurement that motivated this module says it is not correct
      ! for H2 at a cool base (design sec. 3.1, from the P4 block of
      ! docs/oxygen_chemistry_new_plan.md): at the HD 189733 b 1-microbar
      ! level tau_chem(H2)/tau_adv is 0.20 to 1.67 on H/v, and 16.4 when the
      ! two are compared at the same pressure, against 1e-3 on HD 209458 b.
      !
      ! WHAT M3 THEN MEASURED IN EXHALE'S OWN STRUCTURE, AND IT IS NOT THE
      ! SAME NUMBER.  That ratio was built from the PHOTOCHEMICAL model's
      ! scale height and velocity at that pressure.  In an EXHALE run of the
      ! same planet the base cell has v = 0 EXACTLY -- the lower boundary
      ! condition puts it there -- and the flow time r/|v| in the cells above
      ! it is 1e8-2e8 s against tau_chem(H2) = 178-552 s.  The diffusive time
      ! of the base cell is 8.0e5 s with K_zz = 0, i.e. 4500 chemical times,
      ! and 634 s with K_zz = 1e9, i.e. 3.6.  So in this code's structure the
      ! base partition is a LOCAL quantity, and transport is not the term
      ! that decides it; what transport does decide is the profile of a
      ! species with no chemistry of its own, which here is CO.  The M3
      ! result block of the design document carries the numbers and what
      ! follows from them.
      !
      ! That is why every run writes tau_chem, tau_adv AND the diffusive time
      ! of each cell to output/Oxygen_chemistry.txt: the regime is a
      ! measurement of the run, not a property of the option.
      !
      ! ---------------------------------------------------------------
      ! 2. WHICH SPECIES ARE TRANSPORTED, AND WHY THE OTHERS ARE NOT
      !
      ! Transported: H2, OH, H2O, CO.  Between them they carry every H
      ! nucleus that is not atomic or in a molecular ion, and every O and C
      ! nucleus that is not a free atom.  OH is in the set even though its
      ! lifetime is short: measured over the H2 layer of the HD 189733 b
      ! run, OH lives 4.0e-3 to 0.40 s against O1 plus its photolysis while
      ! H2 lives 33 to 553 s against the oxygen cycle, i.e. four decades
      ! longer.  It is carried because it holds up to 14% of the oxygen
      ! where the carriers live (against 43% for H2O and 55% for CO), so
      ! leaving it out would open a hole in the oxygen closure, and because
      ! an implicit step carries a short-lived species at no cost: it simply
      ! relaxes to its own local balance.
      !
      ! NOT transported, and closed by their local balance instead:
      !  - H2+, H3+, HeH+.  Measured on the HD 209458 b molecular example
      !    over the cells that hold H2 (r = 1.000-1.018 R_p), the three of
      !    them together carry at most 1.2e-9 of the H nuclei.  Leaving them
      !    local is therefore a closure statement at that level and not a
      !    leak: the H-nucleus budget is closed against the ELEMENT total,
      !    not against a sum of transported species, so whatever they hold
      !    is simply not available to the carriers.
      !  - O(1D).  Its lifetime is 1.6e-3 s at the HD 189733 b base and it
      !    has one sink; the design's sec. 2.3 states the case, and it is the
      !    same argument as above rather than a different one.
      !  - H I / H II and the ion stages of every element.  Their balance is
      !    local, and gate G2 requires that limit to be reproduced.  They are
      !    where the transported carriers' nuclei come from and go to: the
      !    write-back below moves nuclei between a carrier and the closure
      !    species of its element, so no element can drift.
      !
      ! The one term this drops that is not small everywhere is the AMBIPOLAR
      ! field.  All four transported species are neutral, so eq. (2) of the
      ! design has no (Z e E)/(k T) term for any of them; the field is not
      ! neglected, it is absent.  Had the molecular ions been transported it
      ! would not be, and settling_coefficient in binary_element_diffusion
      ! already builds it.
      !
      ! ---------------------------------------------------------------
      ! 3. THE COEFFICIENTS
      !
      ! MOLECULAR DIFFUSION D_i.  Blanc's law over the background carriers
      ! (H I, H II, H2, H2+, H3+, He I, He II, He III, HeH+), each pair from
      ! the rigid-sphere coefficient of Banks & Kockarts (1973) that
      ! binary_element_diffusion already carries
      ! (hard_sphere_pair_diffusion, reused rather than restated).
      !
      ! What that leaves out is the ION-NEUTRAL polarization enhancement.
      ! For an ion colliding with one of these molecules the induced-dipole
      ! (Langevin) cross section is larger than the rigid-sphere one, so the
      ! rigid-sphere D is an UPPER bound on that pair.  Two reasons it is
      ! used anyway, and both are measurements rather than preferences.
      ! First, the transported species only exist where the gas is neutral:
      ! over the cells that hold H2 on the HD 209458 b molecular example the
      ! electron fraction per H nucleus is at most 1.0e-7 (3.1e-11 at the
      ! base cell), so the ion pairs carry that share of the Blanc friction
      ! sum and the choice of coefficient for them is numerically
      ! irrelevant.  Second, for H2O the induced-dipole formula would be the
      ! wrong one in any case: a permanent dipole of 1.85 D dominates the
      ! ion-molecule capture rate, and the locked-dipole treatment that needs
      ! is not in this code.  Putting an induced polarizability there would
      ! look more careful and be less honest.
      !
      ! EDDY DIFFUSION K_zz.  Read from kzz_cell and nowhere else.  This
      ! module defines no eddy coefficient and no key for one (design
      ! sec. 3.5): kzz_cell is filled once by eddy_diffusion_on_grid, from
      ! the scalar He_Kzz or from the lower-atmosphere profile's own K_zz
      ! column, and every diffusive term in the code reads that one array.
      ! With no profile and no He_Kzz it is zero, and the transport is then
      ! pure molecular diffusion -- which on a lower atmosphere is the wrong
      ! limit, since eddy mixing is what holds the composition well mixed
      ! below the homopause.  The run says so at startup rather than
      ! refusing, so that the K_zz = 0 limit stays available as a test.
      !
      ! SETTLING.  G_i = (m_i - m_bar) g/(k T), with m_bar the mean mass per
      ! gas particle (rho/n_tot) and g from the same Dphi the hydro uses.
      ! Thermal diffusion is not carried (alpha_T = 0): there is no evaluated
      ! alpha_T for H2O, OH or CO in an H2/He bath, and the element operator
      ! makes the same choice for the same reason.
      !
      ! ---------------------------------------------------------------
      ! 4. THE CHEMISTRY IS THE SAME CHEMISTRY
      !
      ! P_i - L_i is not written here.  It is obtained by calling the very
      ! rows the local equilibrium solve uses -- mol_heh_rows row 4 for H2
      ! and oxygen_carrier_rows for OH and H2O -- at the trial densities,
      ! with the background frozen.  A rate that changes there changes here,
      ! and the two cannot drift apart the way two transcriptions of one
      ! coefficient do (the R16-R20 lesson recorded in mol_rates.f90).
      ! The Jacobian is a forward difference of those same calls.
      !
      ! CO has no chemical source: it is frozen chemically by decision D4 and
      ! only its transport moves it.  What that decision MEANS changes when
      ! CO is transported, and the change is not cosmetic.  In the local
      ! solve "inert" meant "at its own CO <-> C + O chemical equilibrium",
      ! which dissociates CO above about 4000 K.  Transported, "inert" would
      ! mean INDESTRUCTIBLE, and it was: measured on HD 189733 b before the
      ! ceiling of limit_to_element_budget was added, the wind carried CO out
      ! to 1.67 R_p and 2e4 K and it held 55% of the oxygen there, switching
      ! off the O I and C II cooling of the entire wind.  The ceiling is the
      ! statement the thermodynamics supports and the audited set has no rate
      ! for; its own limits are written at that code site.
      !
      ! WHAT IS FROZEN INSIDE ONE STEP: the ion stages, the molecular ions,
      ! n_e, n_tot, T, and the photoionization rates.  The photolysis rates
      ! are NOT frozen -- they are rebuilt from the current H2O, OH and H2
      ! columns at the top of every step, so the self-shielding of the water
      ! layer responds within the relaxation instead of only between Picard
      ! passes.
      !
      ! ---------------------------------------------------------------
      ! 5. BOUNDARY CONDITIONS (design decision D6)
      !
      ! ELEMENT RESERVOIRS ARE DIRICHLET; THE PARTITION AMONG CARRIERS IS
      ! NOT.  The H, He, O and C totals at the base are boundary data the
      ! handoff legitimately supplies, and they are imposed where they always
      ! were -- through melem_ab, the base density and the EOS.  The split of
      ! each element among its carriers is what this option exists to
      ! compute, so it is not imposed anywhere: the base face carries zero
      ! diffusive flux and the advective inflow there carries the base cell's
      ! own partition (df/dr = 0 at the face), which is the same statement.
      !
      ! That is what makes the A/B gate of design sec. 6.2 a test.  Had the
      ! base H2 fraction been pinned from the same handoff the gate compares
      ! against, the gate would measure nothing.
      !
      ! OUTER: zero gradient, as the element operator does.  The carriers are
      ! negligible there by many decades in every configuration this option
      ! targets; a run in which they are not is outside its validity range
      ! and output/Oxygen_chemistry.txt says so.
      !
      ! ---------------------------------------------------------------
      ! 6. DISCRETIZATION AND THE SIMPLEX
      !
      ! Cell-centred f_i, face-centred fluxes on r_edg, one backward-Euler
      ! step per call, Newton with a backtracking line search, block
      ! tridiagonal in space (4x4 blocks: the transport is diagonal in
      ! species and the chemistry is dense within a cell).  Faces 0 and N
      ! carry zero diffusive flux.  The drift term switches from central to
      ! donor-cell on the same Peclet test the element operator uses.
      !
      ! ADVECTION IS A ONE-SIDED UPWIND DIFFERENCE ON THE CELL VELOCITY, not
      ! on face-averaged velocities, and that is a lesson inherited rather
      ! than rediscovered: binary_element_diffusion's header records that the
      ! face-averaged form lets a cell whose two face velocities straddle
      ! zero lose its advective term entirely, which froze an isolated helium
      ! hole at the breathing HD 209458 b base.
      !
      ! POSITIVITY AND THE ELEMENT SIMPLEX.  f_i >= 0 is enforced on every
      ! Newton update.  After the step the carriers are limited so that the
      ! nuclei they hold do not exceed the element's own total, in H, O and
      ! C alike; every limiter is counted and reported, never silent.  The
      ! write-back then shares the remaining nuclei out over the element's
      ! other species in proportion to what they already held.
      !
      ! WHAT THAT MAKES THIS OPERATOR, STATED PLAINLY.  Each cell's ELEMENT
      ! TOTALS are unchanged by the step: only the partition among an
      ! element's carriers moves.  So element conservation is exact by
      ! construction rather than by cancellation -- which is what gate G5
      ! measures -- and rho is untouched, because bsp_mass of a carrier is
      ! exactly the sum of its nuclei's masses.
      !
      ! It also fixes the approximation.  A water molecule that diffuses out
      ! of a cell really takes its oxygen with it; here the oxygen stays and
      ! reappears as O I.  That is the same statement EXHALE already makes
      ! about every metal when the trace-metal diffusion arm is off -- the
      ! element is slaved to hydrogen at melem_ab -- and it is why this
      ! operator transports SPECIATION and not elements.  Elemental
      ! separation of oxygen is he_metal_diffusion's subject, and the two
      ! are refused together until CO, which carries an oxygen AND a carbon
      ! nucleus, has a rule for following two element factors at once
      ! (input_read states the refusal).  For hydrogen the distinction does
      ! not arise: the H total is the element, so transporting the partition
      ! 2 n_H2/n_H IS transporting what the gate measures.
      !
      ! ---------------------------------------------------------------
      ! 7. VALIDITY
      !
      ! The same range the reaction set has (oxygen_rates header): this is
      ! valid where the base H2/H partition is kinetic and set by the oxygen
      ! cycle.  Transport moves the boundary of validity outward, it does not
      ! remove it -- design sec. 9 -- and the Damkohler profile in
      ! output/Oxygen_chemistry.txt is where each run states its own.
      !
      ! References: Banks & Kockarts (1973), Aeronomy, Part B, ch. 15 (the
      ! diffusive flux and the rigid-sphere coefficient); Draine & Bertoldi
      ! (1996) for the H2 band; the reaction set is
      ! docs/a2_reaction_audit.md.

      use global_parameters
      use grav_func,     only: Dphi
      use species_table, only: n_bsp, bsp_fsp, bsp_nH, bsp_nO, bsp_nC,   &
                               bsp_is_excited_level, bsp_mass,           &
                               isp_HI, isp_HII, isp_HeI, isp_HeII,       &
                               isp_HeIII, isp_HeTR,                      &
                               isp_H2, isp_H2p, isp_H3p, isp_HeHp,       &
                               isp_OH, isp_H2O, isp_CO,                  &
                               n_mion, mion_fsp, melem_i0, melem_top,    &
                               iel_O, iel_C
      use binary_element_diffusion, only: hard_sphere_pair_diffusion
      use ion_cell_state, only: ieq_cell, ion_rates
      use oxygen_rates, only: co_equilibrium_density
      use System_HeH_mol, only: set_mol_coeffs, set_oxygen_coeffs,       &
                                mol_heh_rows, oxygen_carrier_rows,       &
                                oj3, oj4, oj5, oj7
      use water_photolysis, only: n_fuv_band
      use utils_ion_eq, only: fuv_lw_photon_field
      use utils, only: calc_ne
      use ionization_equilibrium, only: bg_cell, bg_ready

      implicit none
      private

      public :: photochemical_transport_step
      public :: relax_photochemical_composition
      public :: carrier_transport_diagnostics, carrier_co_ceiling_cells
      public :: carrier_diffusion_coefficient
      public :: n_carrier, carrier_name, ic_H2, ic_OH, ic_H2O, ic_CO

      integer, parameter :: dp = kind(1.0d0)

      ! The transported set, in the order the solve uses.
      integer, parameter :: n_carrier = 4
      integer, parameter :: ic_H2 = 1, ic_OH = 2, ic_H2O = 3, ic_CO = 4
      character(len=3), parameter :: carrier_name(n_carrier) =           &
           (/ 'H2 ', 'OH ', 'H2O', 'CO ' /)

      ! Background collision partners for Blanc's law.  Masses in m_H, the
      ! same convention bsp_mass uses.
      integer, parameter :: n_bkg = 9
      real(dp), parameter :: bkg_mass(n_bkg) =                           &
           (/ 1.0d0, 1.0d0, 2.0d0, 2.0d0, 3.0d0,                         &
              4.0d0, 4.0d0, 4.0d0, 5.0d0 /)

      ! Atomic mass unit [g], the value binary_element_diffusion uses.
      real(dp), parameter :: m_amu_g = 1.67353284d-24

      ! Newton controls, mirroring solve_mass_fraction.
      integer,  parameter :: newton_maxit  = 30
      real(dp), parameter :: newton_tol    = 1.0d-12
      real(dp), parameter :: newton_floor  = 1.0d-8
      real(dp), parameter :: newton_drop   = 1.0d-6
      integer,  parameter :: newton_halves = 8
      ! Relaxation controls, mirroring relax_element_composition.
      integer,  parameter :: relax_maxstep = 400
      real(dp), parameter :: relax_tol     = 1.0d-10

      ! The background this step holds frozen: the ion stages, the molecular
      ! ions and the free electrons (sec. 4).  Filled once per step by
      ! carrier_state from the f_sp the step was handed, so the chemistry
      ! rows see a state and not a moving target.
      real(dp), dimension(:), allocatable :: cbg_nhii, cbg_nh2p
      real(dp), dimension(:), allocatable :: cbg_nh3p, cbg_nhehp
      real(dp), dimension(:), allocatable :: cbg_nhei, cbg_nheii
      real(dp), dimension(:), allocatable :: cbg_nheiii, cbg_nheiTR
      real(dp), dimension(:), allocatable :: cbg_ne
      ! Ionized stages of oxygen and carbon, frozen with the rest of the
      ! background: they are what separates the element headroom (which the
      ! limiter uses) from the FREE ATOMIC density the reaction rows need.
      real(dp), dimension(:), allocatable :: cbg_nOion, cbg_nCion
      ! Photolysis rates rebuilt from the CURRENT carrier columns at the top
      ! of each step, so the water layer shields itself inside the
      ! relaxation rather than only between Picard passes.
      real(dp), dimension(:),   allocatable :: cph_klw
      real(dp), dimension(:,:), allocatable :: cph_jh2o, cph_joh

      ! Diagnostics of the last step, read by write_output.
      integer,  save :: pct_newton_steps = 0
      real(dp), save :: pct_newton_resid = 0.0d0
      integer,  save :: pct_limited      = 0
      integer,  save :: pct_co_ceiling    = 0
      integer,  save :: pct_worst_j       = 0
      integer,  save :: pct_worst_ic      = 1
      real(dp), save :: pct_worst_limit  = 0.0d0
      ! Molecular diffusion coefficient of each carrier on the grid, kept
      ! from the last step so the run can print the transport time scale
      ! beside the chemical one it already prints.
      real(dp), dimension(:,:), allocatable :: pct_Dco

      contains

      ! Diagnostics of the last transport step: Newton iterations, the
      ! relative residual it stopped at, how many cells needed the element
      ! limiter, and the worst relative overshoot it removed.
      subroutine carrier_transport_diagnostics(nsteps, resid, nlim, worst)
      integer,  intent(out) :: nsteps, nlim
      real(dp), intent(out) :: resid, worst
      nsteps = pct_newton_steps
      resid  = pct_newton_resid
      nlim   = pct_limited
      worst  = pct_worst_limit
      end subroutine carrier_transport_diagnostics

      ! Cells in which the CO thermal ceiling was applied at the last step.
      integer function carrier_co_ceiling_cells() result(n)
      n = pct_co_ceiling
      end function carrier_co_ceiling_cells

      ! Molecular diffusion coefficient of carrier ic in cell j [cm^2/s],
      ! from the last transport step; zero before the first one.
      double precision function carrier_diffusion_coefficient(j, ic)      &
                                result(D)
      integer, intent(in) :: j, ic
      if (allocated(pct_Dco)) then
         D = pct_Dco(j,ic)
      else
         D = 0.0d0
      endif
      end function carrier_diffusion_coefficient

      ! ------------------------------------------------------------- !

      ! One backward-Euler step of the coupled transport-chemistry system,
      ! over the step dt_code of each cell (adimensional, as the hydro's own).
      ! Updates f_sp in place; a no-op unless the oxygen chemistry is on with
      ! its transport, and unless ioniz_eq has run at least once (the frozen
      ! background comes from there).
      subroutine photochemical_transport_step(rho, v, Tcode, f_sp, dt_code)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho, v
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: Tcode
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: dt_code
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp

      real(dp), dimension(1-Ng:N+Ng,n_carrier) :: fc, fc_old
      real(dp), dimension(1-Ng:N+Ng,n_carrier) :: Dco
      real(dp), dimension(1-Ng:N+Ng)           :: ntot, TK, gphys, mbar
      real(dp), dimension(1-Ng:N+Ng)           :: rp, rep, ntv, dt_phys
      real(dp), dimension(1-Ng:N+Ng)           :: nH_free, nO_free, nC_free
      real(dp), dimension(0:N,n_carrier)       :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier)       :: updrf

      if (.not. thereis_oxychem)   return
      if (.not. oxygen_transport)  return
      if (.not. bg_ready)          return

      call carrier_state(rho, Tcode, f_sp, fc, ntot, TK, mbar,           &
                         nH_free, nO_free, nC_free)
      fc_old = fc

      call carrier_geometry(v, dt_code, ntot, rp, rep, ntv, dt_phys,     &
                            gphys)
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      if (.not. allocated(pct_Dco)) allocate(pct_Dco(1-Ng:N+Ng,n_carrier))
      pct_Dco = Dco
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,     &
                                     Agrd, Bdrf, updrf)
      call carrier_photolysis(rho, TK, f_sp)

      call solve_carriers(fc, fc_old, ntot, dt_phys, rp, rep, ntv,       &
                          Agrd, Bdrf, updrf, TK,                         &
                          nH_free, nO_free, nC_free)

      call carrier_write_back(rho, f_sp, fc, ntot, nH_free, nO_free,     &
                              nC_free)

      if (carrier_debug_on()) then
         write(*,'(a,i4,a,es9.2,a,f6.0,a,i0,a,es9.2,a,f6.0,a,es9.2,'//  &
              'a,es9.2,a,es9.2)')                                        &
            ' (carrier_transport) newton ', pct_newton_steps,            &
            ' resid ', pct_newton_resid, ' at cell '//                  &
            trim(carrier_name(pct_worst_ic))//' j=', dble(pct_worst_j),  &
            ' limited ', pct_limited,                                    &
            ' cells, worst ', pct_worst_limit, ', CO ceiling ',          &
            dble(pct_co_ceiling), '; base x_H2 ',                        &
            2.0d0*fc(1,ic_H2)*ntot(1)/max(nH_free(1), 1.0d-99),          &
            ' <- ', 2.0d0*fc_old(1,ic_H2)*ntot(1)                        &
                    /max(nH_free(1), 1.0d-99),                           &
            ' dt ', dt_phys(1)
      endif

      end subroutine photochemical_transport_step

      ! EXHALE_CARRIER_DEBUG=1 turns on one line per step from the transport
      ! operator: Newton iterations, the relative residual, how many cells
      ! the element limiter touched, and how far the base partition moved.
      logical function carrier_debug_on()
      character(len=8) :: env
      call get_environment_variable('EXHALE_CARRIER_DEBUG', env)
      carrier_debug_on = (trim(env) .eq. '1')
      end function carrier_debug_on

      ! ------------------------------------------------------------- !

      ! The carrier mixing ratios of the current state, and the element
      ! headroom each of them is limited against.
      !
      ! nH_free is the H nuclei NOT held by the ion stages and the molecular
      ! ions, i.e. what H I and the transported carriers share; nO_free and
      ! nC_free are the same for oxygen and carbon against their ion stages.
      ! They are the element totals minus the parts this step cannot move, so
      ! a carrier that stays inside them cannot break an element budget.
      subroutine carrier_state(rho, Tcode, f_sp, fc, ntot, TK, mbar,     &
                               nH_free, nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(out) :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: ntot, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: mbar
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nC_free
      real(dp), dimension(1-Ng:N+Ng) :: nd
      real(dp), dimension(1-Ng:N+Ng,4) :: nmol_l
      real(dp), dimension(1-Ng:N+Ng,n_mion) :: nm_l
      integer :: j, ib, i0, k, im

      if (.not. allocated(cbg_nhii)) then
         allocate(cbg_nhii(1-Ng:N+Ng),  cbg_nh2p(1-Ng:N+Ng),             &
                  cbg_nh3p(1-Ng:N+Ng),  cbg_nhehp(1-Ng:N+Ng),            &
                  cbg_nhei(1-Ng:N+Ng),  cbg_nheii(1-Ng:N+Ng),            &
                  cbg_nheiii(1-Ng:N+Ng),cbg_nheiTR(1-Ng:N+Ng),           &
                  cbg_ne(1-Ng:N+Ng),                                     &
                  cbg_nOion(1-Ng:N+Ng), cbg_nCion(1-Ng:N+Ng))
         allocate(cph_klw(1-Ng:N+Ng),                                    &
                  cph_jh2o(1-Ng:N+Ng,n_fuv_band),                        &
                  cph_joh(1-Ng:N+Ng,n_fuv_band))
         cph_klw  = 0.0d0
         cph_jh2o = 0.0d0
         cph_joh  = 0.0d0
      endif

      nd = rho*n0                      ! cm^-3 per unit of f_sp

      ! The frozen background of this step.
      cbg_nhii   = f_sp(:,isp_HII)  *nd
      cbg_nh2p   = f_sp(:,isp_H2p)  *nd
      cbg_nh3p   = f_sp(:,isp_H3p)  *nd
      cbg_nhehp  = f_sp(:,isp_HeHp) *nd
      cbg_nhei   = f_sp(:,isp_HeI)  *nd
      cbg_nheii  = f_sp(:,isp_HeII) *nd
      cbg_nheiii = f_sp(:,isp_HeIII)*nd
      cbg_nheiTR = 0.0d0
      if (thereis_HeITR) cbg_nheiTR = f_sp(:,isp_HeTR)*nd
      do im = 1, n_mion
         nm_l(:,im) = f_sp(:,mion_fsp(im))
      enddo
      nmol_l(:,1) = f_sp(:,isp_H2)
      nmol_l(:,2) = f_sp(:,isp_H2p)
      nmol_l(:,3) = f_sp(:,isp_H3p)
      nmol_l(:,4) = f_sp(:,isp_HeHp)
      call calc_ne(f_sp(:,isp_HII), f_sp(:,isp_HeII), f_sp(:,isp_HeIII), &
                   cbg_ne, nm_l, nmol_l)
      cbg_ne = cbg_ne*nd
      cbg_nOion = 0.0d0
      i0 = melem_i0(iel_O)
      do k = 1, melem_top(iel_O)
         cbg_nOion = cbg_nOion + f_sp(:,mion_fsp(i0+k))*nd
      enddo
      cbg_nCion = 0.0d0
      i0 = melem_i0(iel_C)
      do k = 1, melem_top(iel_C)
         cbg_nCion = cbg_nCion + f_sp(:,mion_fsp(i0+k))*nd
      enddo

      ! Gas-particle density and the mean mass per particle.  Both come from
      ! the cell state ioniz_eq last solved, so they are the same numbers the
      ! EOS is using; freezing them inside one step is stated in sec. 4.
      do j = 1-Ng, N+Ng
         ntot(j) = bg_cell(j)%ntot
         TK(j)   = bg_cell(j)%T_K
      enddo
      do j = 1-Ng, N+Ng
         if (ntot(j) .gt. 0.0d0) then
            mbar(j) = rho(j)*n0*m_amu_g/ntot(j)
         else
            mbar(j) = m_amu_g
         endif
      enddo

      fc(:,ic_H2)  = f_sp(:,isp_H2) *nd/max(ntot, 1.0d-99)
      fc(:,ic_OH)  = f_sp(:,isp_OH) *nd/max(ntot, 1.0d-99)
      fc(:,ic_H2O) = f_sp(:,isp_H2O)*nd/max(ntot, 1.0d-99)
      fc(:,ic_CO)  = f_sp(:,isp_CO) *nd/max(ntot, 1.0d-99)
      where (fc .lt. 0.0d0) fc = 0.0d0

      ! ---- element headroom -------------------------------------------
      ! What the transported carriers of each element may hold.  It is the
      ! ELEMENT TOTAL, not the neutral stage: transport moves nuclei, and
      ! which stage a nucleus sits in is re-solved from scratch by the very
      ! next ionization sweep.  Limiting the carriers to the neutral stage
      ! instead would cap them wherever the element happens to be ionized --
      ! measured on the HD 209458 b molecular gate, that clamped 150-240 of
      ! 503 cells by up to 2% every step, which is a boundary of the
      ! bookkeeping and not of the physics.
      !
      ! The one nucleus this step genuinely cannot move is the H in HeH+,
      ! because that molecule also carries a helium nucleus and rescaling it
      ! would move helium.  It is therefore taken out of the hydrogen
      ! headroom rather than shared.
      nH_free = 0.0d0
      nO_free = 0.0d0
      nC_free = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nH(ib) .gt. 0)                                          &
            nH_free = nH_free + dble(bsp_nH(ib))*f_sp(:,bsp_fsp(ib))*nd
         if (bsp_nO(ib) .gt. 0)                                          &
            nO_free = nO_free + dble(bsp_nO(ib))*f_sp(:,bsp_fsp(ib))*nd
         if (bsp_nC(ib) .gt. 0)                                          &
            nC_free = nC_free + dble(bsp_nC(ib))*f_sp(:,bsp_fsp(ib))*nd
      enddo
      ! The metal ion stages hold the rest of the oxygen and the carbon.
      i0 = melem_i0(iel_O)
      do k = 0, melem_top(iel_O)
         nO_free = nO_free + f_sp(:,mion_fsp(i0+k))*nd
      enddo
      i0 = melem_i0(iel_C)
      do k = 0, melem_top(iel_C)
         nC_free = nC_free + f_sp(:,mion_fsp(i0+k))*nd
      enddo
      nH_free = nH_free - f_sp(:,isp_HeHp)*nd
      where (nH_free .lt. 0.0d0) nH_free = 0.0d0
      where (nO_free .lt. 0.0d0) nO_free = 0.0d0
      where (nC_free .lt. 0.0d0) nC_free = 0.0d0

      end subroutine carrier_state

      ! ------------------------------------------------------------- !

      ! Grid, step and gravity in physical units, and the advective
      ! coefficient n_tot v of the one-sided upwind term.
      subroutine carrier_geometry(v, dt_code, ntot, rp, rep, ntv,        &
                                  dt_phys, gphys)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: v
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: dt_code, ntot
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: rp, rep, ntv
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: dt_phys, gphys
      integer :: j
      rp      = r*R0
      rep     = r_edg*R0
      ntv     = ntot*v*v0
      dt_phys = dt_code*R0/v0
      where (dt_phys .lt. 1.0d-30) dt_phys = 1.0d-30
      do j = 1-Ng, N+Ng
         gphys(j) = Dphi(r(j))*v0*v0/R0
      enddo
      end subroutine carrier_geometry

      ! ------------------------------------------------------------- !

      ! Molecular diffusion coefficient of each transported carrier in the
      ! background mixture, by Blanc's law over the rigid-sphere pairs
      ! (module header sec. 3).
      subroutine carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: ntot
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(out) :: Dco
      real(dp) :: y(n_bkg), ysum, fric, Dpair, nd, nc
      integer  :: j, ic, is
      integer  :: bkg_isp(n_bkg)
      bkg_isp = (/ isp_HI, isp_HII, isp_H2, isp_H2p, isp_H3p,            &
                   isp_HeI, isp_HeII, isp_HeIII, isp_HeHp /)
      do j = 1-Ng, N+Ng
         nd = rho(j)*n0
         ysum = 0.0d0
         do is = 1, n_bkg
            y(is) = max(f_sp(j,bkg_isp(is))*nd, 0.0d0)
            ysum  = ysum + y(is)
         enddo
         if (ysum .le. 0.0d0) then
            Dco(j,:) = 0.0d0
            cycle
         endif
         y  = y/ysum
         nc = max(ntot(j), 1.0d0)
         do ic = 1, n_carrier
            fric = 0.0d0
            do is = 1, n_bkg
               if (y(is) .le. 0.0d0) cycle
               Dpair = hard_sphere_pair_diffusion(max(TK(j), 1.0d0), nc, &
                            carrier_mass_amu(ic), bkg_mass(is))
               fric  = fric + y(is)/max(Dpair, 1.0d-99)
            enddo
            if (fric .gt. 0.0d0) then
               Dco(j,ic) = 1.0d0/fric
            else
               Dco(j,ic) = 0.0d0
            endif
         enddo
      enddo
      end subroutine carrier_diffusivities

      ! Mass of a transported carrier [m_H], read from the species table so
      ! that this module and calc_rho cannot disagree about it.
      double precision function carrier_mass_amu(ic) result(m)
      integer, intent(in) :: ic
      select case (ic)
      case (ic_H2)
         m = bsp_mass(7)
      case (ic_OH)
         m = bsp_mass(11)
      case (ic_H2O)
         m = bsp_mass(12)
      case default
         m = bsp_mass(13)
      end select
      end function carrier_mass_amu

      ! ------------------------------------------------------------- !

      ! Face gradient and drift coefficients of each carrier, and the
      ! Peclet-hybrid donor switch.  Faces 0 and N are left at zero, which
      ! IS the zero-diffusive-flux boundary of sec. 5.
      subroutine carrier_face_coefficients(ntot, TK, mbar, gphys, Dco,   &
                                           rp, Agrd, Bdrf, updrf)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: ntot, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: mbar, gphys
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(in)  :: Dco
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rp
      real(dp), dimension(0:N,n_carrier),       intent(out) :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier),       intent(out) :: updrf
      real(dp), dimension(1-Ng:N+Ng) :: Gco
      real(dp) :: dr_f, ntf, Df, Gf, Kf
      integer  :: j, ic

      Agrd  = 0.0d0
      Bdrf  = 0.0d0
      updrf = 0

      do ic = 1, n_carrier
         ! Settling coefficient of this carrier, G = (m_i - m_bar) g/(k T).
         do j = 1-Ng, N+Ng
            Gco(j) = (carrier_mass_amu(ic)*m_amu_g - mbar(j))*gphys(j)   &
                     /(kb_erg*max(TK(j), 1.0d0))
         enddo
         do j = 1, N-1
            dr_f = max(rp(j+1) - rp(j), 1.0d0)
            ntf  = 0.5d0*(ntot(j) + ntot(j+1))
            Df   = 0.5d0*(Dco(j,ic) + Dco(j+1,ic))
            Gf   = 0.5d0*(Gco(j) + Gco(j+1))
            Kf   = 0.5d0*(kzz_cell(j) + kzz_cell(j+1))
            Agrd(j,ic) = ntf*(Df + Kf)/dr_f
            Bdrf(j,ic) = ntf*Df*Gf
            if (abs(Bdrf(j,ic))*dr_f .le. 2.0d0*ntf*(Df + Kf)) then
               updrf(j,ic) = 0        ! central: the drift is resolved
            else if (Bdrf(j,ic) .ge. 0.0d0) then
               updrf(j,ic) = -1       ! settling inward: donor is cell j+1
            else
               updrf(j,ic) = +1       ! rising outward: donor is cell j
            endif
         enddo
      enddo
      end subroutine carrier_face_coefficients

      ! ------------------------------------------------------------- !

      ! Rebuild the photolysis and Lyman-Werner rates from the CURRENT
      ! carrier columns and store them on the frozen cell state, so the
      ! water layer shields itself within the relaxation (sec. 4).  Uses the
      ! opa_pf weighting the last ioniz_eq left, like every other column.
      subroutine carrier_photolysis(rho, TK, f_sp)
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho, TK
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real(dp), dimension(1-Ng:N+Ng) :: nH2, nH2O, nOH
      real(dp), dimension(1-Ng:N+Ng) :: c1, c2, c3, fsh, trl, klw
      real(dp), dimension(1-Ng:N+Ng,n_fuv_band) :: tau_b, jh2o, joh
      real(dp) :: amax
      integer  :: j
      nH2  = rho*n0*f_sp(:,isp_H2)
      nH2O = rho*n0*f_sp(:,isp_H2O)
      nOH  = rho*n0*f_sp(:,isp_OH)
      call fuv_lw_photon_field(nH2, nH2O, nOH, TK, c1, c2, c3,           &
                               fsh, trl, tau_b, klw, jh2o, joh, amax)
      do j = 1-Ng, N+Ng
         cph_klw(j)      = klw(j)
         cph_jh2o(j,:)   = jh2o(j,:)
         cph_joh(j,:)    = joh(j,:)
      enddo
      end subroutine carrier_photolysis

      ! ------------------------------------------------------------- !

      ! Chemical production minus loss of the four carriers in cell j, at the
      ! trial densities nc(1..4) [cm^-3].  The rows are the SAME rows the
      ! local equilibrium solve uses (sec. 4); nothing is rewritten here.
      subroutine carrier_source(j, nc, nH_free, nO_free, src)
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: nc(n_carrier), nH_free, nO_free
      real(dp), intent(out) :: src(n_carrier)
      real(dp) :: fv(10)
      real(dp) :: n_hi, n_hii, n_h2p, n_h3p, n_hehp
      real(dp) :: n_hei, n_heii, n_heiii, n_heiTR, n_heiSI, n_e, n_o0

      n_hii   = cbg_nhii(j)
      n_h2p   = cbg_nh2p(j)
      n_h3p   = cbg_nh3p(j)
      n_hehp  = cbg_nhehp(j)
      n_heii  = cbg_nheii(j)
      n_heiii = cbg_nheiii(j)
      n_heiTR = cbg_nheiTR(j)
      n_hei   = cbg_nhei(j)
      n_heiSI = max(n_hei - n_heiTR, 0.0d0)
      n_e     = cbg_ne(j)

      ! Atomic H closes the hydrogen budget: what the element has left after
      ! the ions, the molecular ions and the trial carriers.  This is the
      ! coupling that makes the H2 source a function of the unknown.
      ! nH_free is the element headroom, so the frozen ions come off here.
      n_hi = max(nH_free - n_hii - 2.0d0*n_h2p - 3.0d0*n_h3p             &
                         - 2.0d0*nc(ic_H2) - nc(ic_OH)                   &
                         - 2.0d0*nc(ic_H2O), 0.0d0)
      ! Free atomic oxygen closes the oxygen budget the same way: the
      ! element total less its ionized stages and less the carriers.
      n_o0 = max(nO_free - cbg_nOion(j) - nc(ic_OH) - nc(ic_H2O)         &
                         - nc(ic_CO), 0.0d0)

      ! The cell state of the last sweep, with only the photolysis rates
      ! refreshed: n_ofam and n_co are carried over unchanged because the
      ! rows called below take the densities they need as arguments.
      ieq_cell       = bg_cell(j)
      ieq_cell%k_LW  = cph_klw(j)
      call set_mol_coeffs(ieq_cell%T_K, ieq_cell%ntot)
      call set_oxygen_coeffs(ieq_cell%T_K, cph_jh2o(j,:), cph_joh(j,:))

      fv = 0.0d0
      call mol_heh_rows(fv, n_hi, n_hii, nc(ic_H2), n_h2p, n_h3p,        &
                        n_hehp, n_heiSI, n_heiTR, n_heii, n_heiii, n_e,  &
                        ieq_cell%ntot,                                   &
                        ieq_cell%P_HI, ieq_cell%P_HeI, ieq_cell%P_HeII,  &
                        ieq_cell%P_HeITR, ieq_cell%P_H2, ieq_cell%k_LW,  &
                        ieq_cell%rchiiB, ieq_cell%rcheiiB,               &
                        ieq_cell%rcheiiiB, ieq_cell%rcheiTR,             &
                        ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,           &
                        ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR,       &
                        ieq_cell%q13, ieq_cell%q31a, ieq_cell%q31b,      &
                        ieq_cell%Q31, ieq_cell%A31)
      ! The photolysis channels set_oxygen_coeffs built for this cell: the
      ! transport operator evaluates the chemistry at the cell's own,
      ! unscaled radiation field.
      call oxygen_carrier_rows(fv, 9, n_hi, nc(ic_H2), nc(ic_OH),        &
                               nc(ic_H2O), n_o0, oj3, oj4, oj5, oj7)

      src(ic_H2)  = fv(4)
      src(ic_OH)  = fv(9)
      src(ic_H2O) = fv(10)
      src(ic_CO)  = 0.0d0        ! chemically frozen (decision D4)

      end subroutine carrier_source

      ! ------------------------------------------------------------- !

      ! Face flux of one carrier and its two derivatives.  The drift term is
      ! linear in f (not the X(1-X) form of the element operator, whose two
      ! factors exist because that variable is one half of a binary pair).
      subroutine carrier_face_flux(fl, fr, Agr, Bst, upw, Jf, dJl, dJr)
      real(dp), intent(in)  :: fl, fr, Agr, Bst
      integer,  intent(in)  :: upw
      real(dp), intent(out) :: Jf, dJl, dJr
      real(dp) :: wl, wr
      if (upw .eq. 0) then
         wl = 0.5d0
         wr = 0.5d0
      else if (upw .lt. 0) then
         wl = 0.0d0
         wr = 1.0d0
      else
         wl = 1.0d0
         wr = 0.0d0
      endif
      Jf  = -Agr*(fr - fl) - Bst*(wl*fl + wr*fr)
      dJl =  Agr - Bst*wl
      dJr = -Agr - Bst*wr
      end subroutine carrier_face_flux

      ! ------------------------------------------------------------- !

      ! Residual of the backward-Euler transport-chemistry system, and the
      ! face-flux derivatives the Jacobian needs.  rnorm is the largest
      ! RELATIVE row imbalance, built the way composition_residual builds
      ! its own: the row divided by the sum of the magnitudes of its terms.
      subroutine carrier_residual(fc, fc_old, ntot, dt_phys, rp, rep,    &
                                  ntv, Agrd, Bdrf, updrf,                &
                                  nH_free, nO_free, nC_free,             &
                                  res, Jf, dJl, dJr, rnorm)
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(in)  :: fc, fc_old
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: ntot
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: dt_phys
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rp, rep, ntv
      real(dp), dimension(0:N,n_carrier),       intent(in)  :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier),       intent(in)  :: updrf
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nC_free
      real(dp), dimension(1:N,n_carrier),       intent(out) :: res
      real(dp), dimension(0:N,n_carrier),       intent(out) :: Jf, dJl, dJr
      real(dp),                                 intent(out) :: rnorm
      ! Cell and species carrying the worst relative row imbalance, for the
      ! debug line: a Newton that stops short says WHERE it stopped.
      real(dp) :: nc(n_carrier), src(n_carrier)
      real(dp) :: Kj, sL, sR, cadv, dsc, rr, nfl
      integer  :: j, ic

      Jf  = 0.0d0
      dJl = 0.0d0
      dJr = 0.0d0
      do ic = 1, n_carrier
         do j = 1, N-1
            call carrier_face_flux(fc(j,ic), fc(j+1,ic), Agrd(j,ic),     &
                                   Bdrf(j,ic), updrf(j,ic),              &
                                   Jf(j,ic), dJl(j,ic), dJr(j,ic))
         enddo
      enddo

      rnorm = 0.0d0
      do j = 1, N
         nc = fc(j,:)*ntot(j)
         call carrier_source(j, nc, nH_free(j), nO_free(j), src)
         Kj = 1.0d0/(rp(j)**2*max(rep(j) - rep(j-1), 1.0d0))
         sL = rep(j-1)**2
         sR = rep(j)**2
         do ic = 1, n_carrier
            rr  = ntot(j)*(fc(j,ic) - fc_old(j,ic))/dt_phys(j)
            dsc = abs(ntot(j)*fc(j,ic)/dt_phys(j))                       &
                + abs(ntot(j)*fc_old(j,ic)/dt_phys(j))
            rr  = rr + Kj*(sR*Jf(j,ic) - sL*Jf(j-1,ic))
            dsc = dsc + Kj*(sR*abs(Jf(j,ic)) + sL*abs(Jf(j-1,ic)))
            ! One-sided upwind on the CELL velocity; the base face carries
            ! no advective term, which is the zero-flux partition boundary.
            if (ntv(j) .ge. 0.0d0) then
               if (j .gt. 1) then
                  cadv = ntv(j)/max(rp(j) - rp(j-1), 1.0d0)
                  rr   = rr + cadv*(fc(j,ic) - fc(j-1,ic))
                  dsc  = dsc + abs(cadv)*max(abs(fc(j,ic)),              &
                                             abs(fc(j-1,ic)))
               endif
            else
               if (j .lt. N) then
                  cadv = ntv(j)/max(rp(j+1) - rp(j), 1.0d0)
                  rr   = rr + cadv*(fc(j+1,ic) - fc(j,ic))
                  dsc  = dsc + abs(cadv)*max(abs(fc(j,ic)),              &
                                             abs(fc(j+1,ic)))
               endif
            endif
            rr  = rr - src(ic)
            dsc = dsc + abs(src(ic))
            ! ABSOLUTE FLOOR ON THE ROW SCALE, and it is not cosmetic. A
            ! relative residual with no floor asks a row whose species is
            ! 1e-28 of its element to balance to 1e-12 of ITSELF, and out in
            ! the wind every term of such a row is round-off: measured on
            ! the HD 209458 b example, H2O jumped four decades from cell to
            ! cell at 1e-28 and the Newton spent its whole iteration budget
            ! there at a relative residual of 1 while every cell that
            ! carries a molecule was already at 1e-13. The floor is 1e-20 of
            ! the element the carrier belongs to, converted to the row's own
            ! units by dividing by the step: a density that small cannot
            ! change any observable, so a row below it IS converged.
            if (ic .eq. ic_H2) then
               nfl = nH_free(j)
            else if (ic .eq. ic_CO) then
               nfl = min(nO_free(j), nC_free(j))
            else
               nfl = nO_free(j)
            endif
            dsc = dsc + 1.0d-20*nfl/dt_phys(j)
            res(j,ic) = rr
            if (abs(rr)/max(dsc, 1.0d-300) .gt. rnorm) then
               rnorm      = abs(rr)/max(dsc, 1.0d-300)
               pct_worst_j  = j
               pct_worst_ic = ic
            endif
         enddo
      enddo
      end subroutine carrier_residual

      ! ------------------------------------------------------------- !

      ! Newton solve of the coupled system, block-tridiagonal in space with
      ! 4x4 blocks.  The chemistry Jacobian is a forward difference of the
      ! same rows the residual calls, so a change to the network reaches the
      ! Jacobian without a second edit.
      subroutine solve_carriers(fc, fc_old, ntot, dt_phys, rp, rep, ntv, &
                                Agrd, Bdrf, updrf, TK,                   &
                                nH_free, nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(inout) :: fc
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(in)    :: fc_old
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: ntot
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: dt_phys
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rp, rep
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: ntv, TK
      real(dp), dimension(0:N,n_carrier),       intent(in)    :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier),       intent(in)    :: updrf
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free

      real(dp), dimension(1:N,n_carrier)  :: res, rhs
      real(dp), dimension(0:N,n_carrier)  :: Jf, dJl, dJr
      real(dp), dimension(1:N,n_carrier,n_carrier) :: aa, bb, cc
      real(dp), dimension(1:N,n_carrier)  :: dfc
      real(dp), dimension(1-Ng:N+Ng,n_carrier) :: ftry
      real(dp) :: nc(n_carrier), s0(n_carrier), s1(n_carrier)
      real(dp) :: Kj, sL, sR, cadv, dn, nref
      real(dp) :: rnorm, rprev, rstart, rtry, damp
      integer  :: j, ic, kc, it, ihalf

      call carrier_residual(fc, fc_old, ntot, dt_phys, rp, rep, ntv,     &
                            Agrd, Bdrf, updrf, nH_free, nO_free, nC_free,&
                            res, Jf, dJl, dJr, rnorm)
      rstart = rnorm
      rprev  = rnorm
      pct_newton_steps = 0

      do it = 1, newton_maxit
         if (rnorm .le. newton_tol) exit

         ! ---- assemble the block tridiagonal
         aa = 0.0d0
         bb = 0.0d0
         cc = 0.0d0
         do j = 1, N
            Kj = 1.0d0/(rp(j)**2*max(rep(j) - rep(j-1), 1.0d0))
            sL = rep(j-1)**2
            sR = rep(j)**2
            do ic = 1, n_carrier
               bb(j,ic,ic) = ntot(j)/dt_phys(j)
               if (j .lt. N) then
                  bb(j,ic,ic) = bb(j,ic,ic) + Kj*sR*dJl(j,ic)
                  cc(j,ic,ic) = cc(j,ic,ic) + Kj*sR*dJr(j,ic)
               endif
               if (j .gt. 1) then
                  bb(j,ic,ic) = bb(j,ic,ic) - Kj*sL*dJr(j-1,ic)
                  aa(j,ic,ic) = aa(j,ic,ic) - Kj*sL*dJl(j-1,ic)
               endif
               if (ntv(j) .ge. 0.0d0) then
                  if (j .gt. 1) then
                     cadv = ntv(j)/max(rp(j) - rp(j-1), 1.0d0)
                     bb(j,ic,ic) = bb(j,ic,ic) + cadv
                     aa(j,ic,ic) = aa(j,ic,ic) - cadv
                  endif
               else
                  if (j .lt. N) then
                     cadv = ntv(j)/max(rp(j+1) - rp(j), 1.0d0)
                     bb(j,ic,ic) = bb(j,ic,ic) - cadv
                     cc(j,ic,ic) = cc(j,ic,ic) + cadv
                  endif
               endif
            enddo
            ! ---- chemistry block, by forward difference of the same rows
            !
            ! THE STEP NEEDS A FLOOR TIED TO THE ELEMENT, and that was a
            ! real defect rather than a precaution. Every source term is
            ! built from densities of order the element total, so a
            ! perturbation far below the round-off of THOSE terms returns
            ! noise for a derivative. Where a carrier is empty -- OH high in
            ! the wind, which photolysis keeps at 1e-30 of its element --
            ! a step proportional to the carrier itself is exactly that,
            ! and the Newton then churned to its iteration cap with a
            ! relative residual of 1 in that one cell (measured on
            ! HD 189733 b). The floor below is 1e-12 of the element the
            ! carrier belongs to, which is resolvable against a double's
            ! 1e-16 while still being a small perturbation of the row.
            nc = fc(j,:)*ntot(j)
            call carrier_source(j, nc, nH_free(j), nO_free(j), s0)
            do kc = 1, n_carrier
               if (kc .eq. ic_H2) then
                  nref = nH_free(j)
               else if (kc .eq. ic_CO) then
                  nref = min(nO_free(j), nC_free(j))
               else
                  nref = nO_free(j)
               endif
               dn = max(1.0d-6*abs(nc(kc)), 1.0d-12*nref, 1.0d-300)
               nc(kc) = nc(kc) + dn
               call carrier_source(j, nc, nH_free(j), nO_free(j), s1)
               nc(kc) = nc(kc) - dn
               do ic = 1, n_carrier
                  bb(j,ic,kc) = bb(j,ic,kc)                              &
                              - (s1(ic) - s0(ic))/dn*ntot(j)
               enddo
            enddo
         enddo

         rhs = -res
         call block_thomas(aa, bb, cc, rhs, dfc)

         ! ---- damped update, with f >= 0 enforced on every try
         damp = 1.0d0
         do ihalf = 0, newton_halves
            ftry = fc
            do ic = 1, n_carrier
               do j = 1, N
                  ftry(j,ic) = max(fc(j,ic) + damp*dfc(j,ic), 0.0d0)
               enddo
               ftry(1-Ng:0,ic)   = ftry(1,ic)
               ftry(N+1:N+Ng,ic) = ftry(N,ic)
            enddo
            call carrier_residual(ftry, fc_old, ntot, dt_phys, rp, rep,  &
                                  ntv, Agrd, Bdrf, updrf, nH_free,       &
                                  nO_free, nC_free, res, Jf, dJl, dJr,   &
                                  rtry)
            if (rtry .lt. rnorm .or. ihalf .eq. newton_halves) exit
            damp = 0.5d0*damp
         enddo
         fc    = ftry
         rprev = rnorm
         rnorm = rtry
         pct_newton_steps = it
         if (it .ge. 3 .and. rnorm .gt. 0.5d0*rprev .and.                &
             (rnorm .le. newton_floor .or.                               &
              rnorm .le. newton_drop*rstart)) exit
      enddo
      pct_newton_resid = rnorm

      call limit_to_element_budget(fc, ntot, TK, nH_free, nO_free,        &
                                   nC_free)

      end subroutine solve_carriers

      ! ------------------------------------------------------------- !

      ! Block-tridiagonal Thomas sweep with dense n_carrier x n_carrier
      ! blocks.  Written here because the element operator's scalar Thomas
      ! sweep does not generalize: the species system is diagonal in space
      ! but dense in species through the chemistry.
      subroutine block_thomas(aa, bb, cc, dd, xx)
      real(dp), dimension(1:N,n_carrier,n_carrier), intent(in)  :: aa,bb,cc
      real(dp), dimension(1:N,n_carrier),           intent(in)  :: dd
      real(dp), dimension(1:N,n_carrier),           intent(out) :: xx
      real(dp), dimension(1:N,n_carrier,n_carrier) :: cp
      real(dp), dimension(1:N,n_carrier)           :: dp_
      real(dp) :: mm(n_carrier,n_carrier), rhs(n_carrier,n_carrier+1)
      real(dp) :: piv, fac
      integer  :: j, i, k, l, ip

      cp  = 0.0d0
      dp_ = 0.0d0
      do j = 1, N
         ! mm = bb(j) - aa(j) cp(j-1)
         mm = bb(j,:,:)
         if (j .gt. 1) then
            do i = 1, n_carrier
               do k = 1, n_carrier
                  do l = 1, n_carrier
                     mm(i,k) = mm(i,k) - aa(j,i,l)*cp(j-1,l,k)
                  enddo
               enddo
            enddo
         endif
         ! Solve mm * [cp(j) | dp(j)] = [cc(j) | dd(j) - aa(j) dp(j-1)]
         rhs(:,1:n_carrier) = cc(j,:,:)
         rhs(:,n_carrier+1) = dd(j,:)
         if (j .gt. 1) then
            do i = 1, n_carrier
               do l = 1, n_carrier
                  rhs(i,n_carrier+1) = rhs(i,n_carrier+1)                &
                                     - aa(j,i,l)*dp_(j-1,l)
               enddo
            enddo
         endif
         ! Gaussian elimination with partial pivoting on the small block.
         do i = 1, n_carrier
            ip  = i
            piv = abs(mm(i,i))
            do k = i+1, n_carrier
               if (abs(mm(k,i)) .gt. piv) then
                  piv = abs(mm(k,i))
                  ip  = k
               endif
            enddo
            if (ip .ne. i) then
               do k = 1, n_carrier
                  fac      = mm(i,k)
                  mm(i,k)  = mm(ip,k)
                  mm(ip,k) = fac
               enddo
               do k = 1, n_carrier+1
                  fac       = rhs(i,k)
                  rhs(i,k)  = rhs(ip,k)
                  rhs(ip,k) = fac
               enddo
            endif
            if (mm(i,i) .eq. 0.0d0) mm(i,i) = 1.0d-300
            do k = i+1, n_carrier
               fac = mm(k,i)/mm(i,i)
               if (fac .eq. 0.0d0) cycle
               do l = i, n_carrier
                  mm(k,l) = mm(k,l) - fac*mm(i,l)
               enddo
               do l = 1, n_carrier+1
                  rhs(k,l) = rhs(k,l) - fac*rhs(i,l)
               enddo
            enddo
         enddo
         do l = 1, n_carrier+1
            do i = n_carrier, 1, -1
               fac = rhs(i,l)
               do k = i+1, n_carrier
                  fac = fac - mm(i,k)*rhs(k,l)
               enddo
               rhs(i,l) = fac/mm(i,i)
            enddo
         enddo
         cp(j,:,:) = rhs(:,1:n_carrier)
         dp_(j,:)  = rhs(:,n_carrier+1)
      enddo

      xx(N,:) = dp_(N,:)
      do j = N-1, 1, -1
         do i = 1, n_carrier
            xx(j,i) = dp_(j,i)
            do k = 1, n_carrier
               xx(j,i) = xx(j,i) - cp(j,i,k)*xx(j+1,k)
            enddo
         enddo
      enddo
      end subroutine block_thomas

      ! ------------------------------------------------------------- !

      ! Hold the carriers inside the element simplex: the H nuclei they
      ! carry cannot exceed what the cell has left for them, and the same for
      ! O and C.  Every limiter is counted and the worst overshoot is
      ! reported (design sec. 3.6: no silent clamp).
      subroutine limit_to_element_budget(fc, ntot, TK, nH_free, nO_free,  &
                                         nC_free)
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(inout) :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: ntot, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      real(dp) :: nc(n_carrier), got, sc, over, nco_eq
      logical  :: hit
      integer  :: j, ic
      ! A carrier that holds the WHOLE of its element sits exactly on the
      ! face of the simplex -- all the carbon is CO in a molecular base --
      ! so a limiter that counted every round-off would report half the grid
      ! and mean nothing. Only an overshoot above this relative size is
      ! counted; the clamp itself is applied whatever the size.
      real(dp), parameter :: limit_report = 1.0d-10

      pct_limited     = 0
      pct_worst_limit = 0.0d0
      pct_co_ceiling  = 0
      ! Ghosts take the partition of the cell they mirror -- that is the
      ! zero-gradient outer boundary and the zero-flux inner one -- and then
      ! go through the same limiter as an interior cell.  They HAVE to: a
      ! ghost is denser than the cell it copies, so the same mixing ratio is
      ! a larger density there and can ask for more nuclei than the ghost's
      ! element holds.  Leaving them out of the limiter was measured, and it
      ! broke the oxygen and carbon budgets of the base ghost by 2% in a
      ! 400-step run while every interior cell closed at round-off.
      do ic = 1, n_carrier
         fc(1-Ng:0,ic)   = fc(1,ic)
         fc(N+1:N+Ng,ic) = fc(N,ic)
      enddo
      do j = 1-Ng, N+Ng
         nc  = max(fc(j,:)*ntot(j), 0.0d0)
         hit = .false.
         ! ---- CO cannot be carried where it cannot exist ----------------
         ! Decision D4 makes CO chemically inert, and with the LOCAL solve
         ! that meant "at its own CO <-> C + O chemical equilibrium", which
         ! dissociates it above ~4000 K. Transported, "inert" would mean
         ! INDESTRUCTIBLE: measured on HD 189733 b before this ceiling, the
         ! wind advected CO out to 1.6 R_p and 28 000 K and it took the
         ! WHOLE oxygen and carbon inventory there, switching off the O I
         ! and C II coolants of the entire wind rather than of the
         ! molecular layer alone.
         !
         ! The network has no CO kinetics to destroy it with, so the
         ! statement made here is the one the thermodynamics supports: where
         ! the equilibrium constant says CO cannot exist, it does not, and
         ! the collisional and radiative processes that take it apart at
         ! 2e4 K run far faster than the flow. Where the equilibrium DOES
         ! allow CO the ceiling is inactive and the transported value
         ! stands, which is the whole molecular layer -- so this is a
         ! one-sided constraint, not a return to equilibrium.
         !
         ! WHERE IT OVER-SUPPRESSES. A parcel carried out of the molecular
         ! layer faster than CO can dissociate is QUENCHED above
         ! equilibrium, which is exactly what a transport operator exists to
         ! capture, and this ceiling removes it -- in the temperature
         ! interval where CO dissociation is neither fast nor negligible,
         ! around the 3000-4000 K turnover of co_equilibrium_density. Far
         ! above that interval the destruction is instantaneous and the
         ! ceiling IS the answer; far below it the ceiling is inactive. The
         ! width of the interval it gets wrong cannot be measured without
         ! the rate the audited set does not contain. A real CO network
         ! would replace the ceiling with that rate, and would also have to
         ! carry CO's own photodissociation (water_photolysis.f90 sec. 3).
         nco_eq = co_equilibrium_density(nC_free(j), nO_free(j),          &
                                         max(TK(j), 1.0d0))
         if (nc(ic_CO) .gt. nco_eq) then
            over = nc(ic_CO)/max(nco_eq, 1.0d-300) - 1.0d0
            nc(ic_CO) = nco_eq
            if (over .gt. limit_report) pct_co_ceiling = pct_co_ceiling + 1
         endif
         ! carbon: CO is the only carbon carrier
         if (nc(ic_CO) .gt. nC_free(j)) then
            over = nc(ic_CO)/max(nC_free(j), 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            nc(ic_CO) = nC_free(j)
            if (over .gt. limit_report) hit = .true.
         endif
         ! oxygen: CO already fits, scale the water family into what is left
         got = nc(ic_OH) + nc(ic_H2O) + nc(ic_CO)
         if (got .gt. nO_free(j)) then
            over = got/max(nO_free(j), 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            sc = max(nO_free(j) - nc(ic_CO), 0.0d0)                      &
                 /max(nc(ic_OH) + nc(ic_H2O), 1.0d-300)
            nc(ic_OH)  = nc(ic_OH) *sc
            nc(ic_H2O) = nc(ic_H2O)*sc
            if (over .gt. limit_report) hit = .true.
         endif
         ! hydrogen: scale the three H-bearing carriers together
         got = 2.0d0*nc(ic_H2) + nc(ic_OH) + 2.0d0*nc(ic_H2O)
         if (got .gt. nH_free(j)) then
            over = got/max(nH_free(j), 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            sc = nH_free(j)/max(got, 1.0d-300)
            nc(ic_H2)  = nc(ic_H2) *sc
            nc(ic_OH)  = nc(ic_OH) *sc
            nc(ic_H2O) = nc(ic_H2O)*sc
            if (over .gt. limit_report) hit = .true.
         endif
         if (hit) pct_limited = pct_limited + 1
         fc(j,:) = nc/max(ntot(j), 1.0d-99)
      enddo
      end subroutine limit_to_element_budget

      ! ------------------------------------------------------------- !

      ! Write the solved carriers back into f_sp, moving the nucleus
      ! difference into the closure species of each element so that no
      ! element total and no mass density is changed by the step.  bsp_mass
      ! of a carrier is exactly the sum of its nuclei's masses, so the mass
      ! moved out of H I, O I and C I is the mass moved into the carriers.
      subroutine carrier_write_back(rho, f_sp, fc, ntot, nH_free,        &
                                    nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho
      real(dp), dimension(1-Ng:N+Ng,n_carrier), intent(in)    :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: ntot
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp) :: nd, nH2n, nOHn, nH2On, nCOn, rest, held, sc
      integer  :: j, k, i0

      do j = 1-Ng, N+Ng
         nd = rho(j)*n0
         if (nd .le. 0.0d0) cycle
         nH2n  = fc(j,ic_H2) *ntot(j)
         nOHn  = fc(j,ic_OH) *ntot(j)
         nH2On = fc(j,ic_H2O)*ntot(j)
         nCOn  = fc(j,ic_CO) *ntot(j)
         f_sp(j,isp_H2)  = nH2n /nd
         f_sp(j,isp_OH)  = nOHn /nd
         f_sp(j,isp_H2O) = nH2On/nd
         f_sp(j,isp_CO)  = nCOn /nd

         ! The nuclei the carriers did NOT take are shared out over the
         ! element's other species IN PROPORTION TO WHAT THEY ALREADY HELD.
         ! Two reasons for the proportional rule rather than putting the
         ! whole remainder in the neutral stage: it cannot drive a stage
         ! negative, and it leaves the ionization fraction of the cell where
         ! it was, which is the right starting point for the sweep that
         ! re-solves it a moment later.  HeH+ is excluded from the hydrogen
         ! share because rescaling it would move a helium nucleus too.
         rest = nH_free(j) - 2.0d0*nH2n - nOHn - 2.0d0*nH2On
         held = f_sp(j,isp_HI)*nd + f_sp(j,isp_HII)*nd                   &
              + 2.0d0*f_sp(j,isp_H2p)*nd + 3.0d0*f_sp(j,isp_H3p)*nd
         if (held .gt. 0.0d0 .and. rest .gt. 0.0d0) then
            sc = rest/held
            f_sp(j,isp_HI)  = f_sp(j,isp_HI) *sc
            f_sp(j,isp_HII) = f_sp(j,isp_HII)*sc
            f_sp(j,isp_H2p) = f_sp(j,isp_H2p)*sc
            f_sp(j,isp_H3p) = f_sp(j,isp_H3p)*sc
         else if (held .le. 0.0d0) then
            f_sp(j,isp_HI) = max(rest, 0.0d0)/nd
         else
            f_sp(j,isp_HI)  = 0.0d0
            f_sp(j,isp_HII) = 0.0d0
            f_sp(j,isp_H2p) = 0.0d0
            f_sp(j,isp_H3p) = 0.0d0
         endif

         rest = nO_free(j) - nOHn - nH2On - nCOn
         i0   = melem_i0(iel_O)
         held = 0.0d0
         do k = 0, melem_top(iel_O)
            held = held + f_sp(j,mion_fsp(i0+k))*nd
         enddo
         if (held .gt. 0.0d0) then
            sc = max(rest, 0.0d0)/held
            do k = 0, melem_top(iel_O)
               f_sp(j,mion_fsp(i0+k)) = f_sp(j,mion_fsp(i0+k))*sc
            enddo
         else
            f_sp(j,mion_fsp(i0)) = max(rest, 0.0d0)/nd
         endif

         rest = nC_free(j) - nCOn
         i0   = melem_i0(iel_C)
         held = 0.0d0
         do k = 0, melem_top(iel_C)
            held = held + f_sp(j,mion_fsp(i0+k))*nd
         enddo
         if (held .gt. 0.0d0) then
            sc = max(rest, 0.0d0)/held
            do k = 0, melem_top(iel_C)
               f_sp(j,mion_fsp(i0+k)) = f_sp(j,mion_fsp(i0+k))*sc
            enddo
         else
            f_sp(j,mion_fsp(i0)) = max(rest, 0.0d0)/nd
         endif
      enddo
      end subroutine carrier_write_back

      ! ------------------------------------------------------------- !

      ! Relax the carriers to their steady state at a fixed wind, for the
      ! Picard loop of the steady solver.  Same shape as
      ! relax_element_composition: a time scale for each cell, grown
      ! geometrically until the state stops moving, with the steady mass
      ! flux in place of the instantaneous velocity so that a breathing base
      ! does not drive the relaxation.
      subroutine relax_photochemical_composition(rho, v, Tcode, f_sp,    &
                                                 omega, drift, nstep)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho, v
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp),                                 intent(in)    :: omega
      real(dp),                                 intent(out)   :: drift
      integer,                                  intent(out)   :: nstep

      real(dp), dimension(1-Ng:N+Ng,n_species) :: f_in
      real(dp), dimension(1-Ng:N+Ng) :: dt_code
      real(dp), dimension(1-Ng:N+Ng,n_carrier) :: fprev, fnow, fentry
      real(dp), dimension(1-Ng:N+Ng,n_carrier) :: Dco
      real(dp), dimension(1-Ng:N+Ng) :: ntot, TK, mbar
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp) :: drj, tdiff, tadv, grow, tscale, dmax, x_ref
      integer  :: j, k, ic

      drift = 0.0d0
      nstep = 0
      if (.not. thereis_oxychem)  return
      if (.not. oxygen_transport) return
      if (.not. bg_ready)         return

      f_in   = f_sp
      tscale = R0/v0
      call carrier_state(rho, Tcode, f_sp, fprev, ntot, TK, mbar,        &
                         nH_free, nO_free, nC_free)
      fentry = fprev
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      do j = 1, N
         drj   = max((r_edg(j) - r_edg(j-1))*R0, 1.0d0)
         tdiff = drj*drj/max(maxval(Dco(j,:)) + kzz_cell(j), 1.0d-30)
         tadv  = drj/max(abs(v(j))*v0, 1.0d-30)
         dt_code(j) = min(tdiff, tadv)/tscale
      enddo
      dt_code(1-Ng:0)   = dt_code(1)
      dt_code(N+1:N+Ng) = dt_code(N)

      grow = 1.0d0
      do k = 1, relax_maxstep
         call photochemical_transport_step(rho, v, Tcode, f_sp,          &
                                           dt_code*grow)
         call carrier_state(rho, Tcode, f_sp, fnow, ntot, TK, mbar,      &
                            nH_free, nO_free, nC_free)
         nstep = k
         dmax  = 0.0d0
         do ic = 1, n_carrier
            do j = 1, N
               dmax = max(dmax, abs(fnow(j,ic) - fprev(j,ic)))
            enddo
         enddo
         x_ref = max(maxval(fnow(1:N,ic_H2)), 1.0d-30)
         if (dmax/x_ref .lt. relax_tol) exit
         fprev = fnow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo

      ! Undamped drift over the WHOLE pass, on the same measure the element
      ! relaxation reports (its own drift is max|X_relaxed - X_old|/X_base);
      ! then the damped update the Picard loop asked for.
      drift = 0.0d0
      do ic = 1, n_carrier
         do j = 1, N
            drift = max(drift, abs(fnow(j,ic) - fentry(j,ic)))
         enddo
      enddo
      drift = drift/max(maxval(fentry(1:N,ic_H2)),                       &
                        maxval(fnow(1:N,ic_H2)), 1.0d-30)
      if (omega .lt. 1.0d0) then
         f_sp(:,isp_H2)  = f_in(:,isp_H2)                                &
                         + omega*(f_sp(:,isp_H2)  - f_in(:,isp_H2))
         f_sp(:,isp_OH)  = f_in(:,isp_OH)                                &
                         + omega*(f_sp(:,isp_OH)  - f_in(:,isp_OH))
         f_sp(:,isp_H2O) = f_in(:,isp_H2O)                               &
                         + omega*(f_sp(:,isp_H2O) - f_in(:,isp_H2O))
         f_sp(:,isp_CO)  = f_in(:,isp_CO)                                &
                         + omega*(f_sp(:,isp_CO)  - f_in(:,isp_CO))
         call carrier_state(rho, Tcode, f_sp, fnow, ntot, TK, mbar,      &
                            nH_free, nO_free, nC_free)
         call limit_to_element_budget(fnow, ntot, TK, nH_free, nO_free,   &
                                      nC_free)
         call carrier_write_back(rho, f_sp, fnow, ntot, nH_free,         &
                                 nO_free, nC_free)
      endif

      end subroutine relax_photochemical_composition

      ! End of module
      end module diffusive_photochemistry

      module binary_element_diffusion
      ! Binary (two-component) H/He element transport for the single-fluid
      ! EXHALE wind.
      !
      ! COMPONENTS.  The gas is treated as two components moving through each
      ! other with a single mass-averaged velocity v (the hydro's velocity):
      !   component 1 : hydrogen carriers together with the trace metals and
      !                 the heavy nuclei bound in the molecular carriers.
      !                 Its mass rho_1 is read from the cell's own species
      !                 (mixture_mass_split), and its mass per hydrogen
      !                 nucleus is rho_1/n_H, which equals the reservoir
      !                 m_1 = mass_per_H_nucleus_without_He() only while the
      !                 cell sits at the reservoir metal/H;
      !   component He: helium, all stages.
      ! With that split rho_1 + rho_He = rho exactly (the metal mass is in rho
      ! under the same eos_metals policy that puts it into rho_1), so the two
      ! diffusive mass fluxes close, J_1 = -J_He, and the composition update
      ! creates or destroys no mass.
      !
      ! STATE.  The transported variable is the helium MASS FRACTION
      !
      !    X = rho_He/rho = m_He n_He / (rho_1 + m_He n_He),   0 <= X <= 1
      !
      ! with n_H, n_He the ELEMENT (nucleus) densities counted over every
      ! species through the bsp_nH / bsp_nHe weights of species_table, so the
      ! count is the same quantity in the atomic, triplet-helium and molecular
      ! regions.  X is bounded at both ends of the composition axis, unlike the
      ! ratio f = n_He/n_H of the earlier trace kernel, which diverges exactly
      ! in the helium-rich limit this module exists for.
      !
      ! TRANSPORT.  Equation (1) of the memo,
      !
      !    d(rho X)/dt + (1/r^2) d/dr [ r^2 ( rho X v + J ) ] = 0
      !
      ! SPLIT IN TWO, AND WHERE EACH HALF LIVES.  The advection of the element
      ! mass fraction is a divergence of the SAME face mass flux the density
      ! rides on, so it is carried inside the Runge-Kutta stages by
      ! species_advection_stage below, on the hydro's own faces, volumes and
      ! time step.  Each element's nucleus total over the domain then changes
      ! only by its boundary fluxes, to round-off, and not to the order of the
      ! hydro's truncation.  What is left for the implicit step here is the
      ! diffusive half,
      !
      !    rho (X^new - X^adv)/dt = -(1/r^2) d/dr [ r^2 J ] ,
      !
      ! taken at the composition the stages advected.  The operator therefore
      ! advects only when it is GIVEN the face mass flux (Frho_in), which is
      ! the fixed-wind relaxation of relax_element_composition, a steady solve
      ! that has no Runge-Kutta stages to ride on; the marching path omits it,
      ! and the stages carry the advection instead.  When it does advect it
      ! forms the term from that face mass flux through the same face
      ! routines the stages use (species_face_fraction, species_face_flux)
      ! and the same spherical geometry, so the relaxation's fixed point is
      ! the zero of the stationary elemental row and not of a second
      ! discretization.
      !
      ! THIS MODULE KNOWS NOTHING OF THE STEADY RESIDUAL: IT IS GIVEN THE
      ! FLUX.  F_rho is a required argument of the fixed-wind relaxation and
      ! of the stationary elemental row, and every caller fetches it from the
      ! mass row of the state it is measuring (face_mass_flux_of_state,
      ! steady_residual.f90) at the call site.  Element transport is a
      ! closure of the mixture and stands below the hydrodynamic solver in
      ! the module order; reaching up to that solver from here inverted the
      ! layering and dragged the whole hydrodynamic tree, LAPACK included,
      ! into every program that links this module, so the standalone
      ! acceptance binary (diffusion_tests.x) could not be linked at all.
      !
      ! The face-flux form removes two failures the earlier cell-velocity
      ! upwind difference had to work around one at a time, because a cell can
      ! now only lose what the mass row says it loses: a cell whose two face
      ! velocities straddle zero keeps an advective term of its own, and a
      ! cell that is outflowing at both faces cannot be evacuated of an
      ! element by a divergence the hydro's own density does not have.
      !
      ! DIFFUSIVE MASS FLUX (memo eqs. 2-4, 6).  For a binary mixture only one
      ! diffusive flux is independent; the helium one is
      !
      !    J = -rho (D_12 + K_zz) dX/dr  -  rho D_12 G X (1-X)
      !    G = [ (m_He - m_c1) m_H g - (Zbar_He - Zbar_c1) eE ] / (k T)
      !        + alpha_T dlnT/dr  -  dln(psi)/dr                      [1/cm]
      !
      ! and the hydrogen-component flux is its negative.  The masses and
      ! charges in G are CARRIER quantities: m_c1 is the mass of component 1
      ! per collision partner and Zbar_c1 its mean charge per collision
      ! partner (see CARRIERS below), which in the atomic region are the mass
      ! and charge per nucleus, m_1 and Zbar_1.
      !
      ! The gradient coefficient of the memo,
      ! A = (m_c1 m_He n_car/mbar_car) D_12 (dx/dX), collapses to rho D_12,
      ! and it does so with the carriers as well as with the nuclei.  Write
      ! n_car = n_1c + n_He for the collision partners, psi = n_1c/n_H for the
      ! carriers each hydrogen nucleus is spread over (1 atomic, 1/2 fully
      ! H2), m_c1 = m_1/psi and mbar_car = rho/n_car.  Then
      ! x = n_He/n_car obeys x = (X/m_He)/(X/m_He + psi (1-X)/m_1), so at
      ! frozen psi -- the chemistry is lagged over the step, as every other
      ! coefficient is -- dx/dX = psi mbar_car^2/(m_1 m_He), while
      ! rho_1 rho_He/(rho x (1-x)) = m_c1 m_He n_car/mbar_car.  The product is
      ! rho D_12 (m_c1 psi/m_1) = rho D_12 exactly, for any psi.  So the
      ! gradient term of the molecular region is the same Fickian rho D dX/dr
      ! as the atomic one.
      !
      ! What the frozen-psi derivative leaves out is the piece of dx/dr the
      ! CHEMISTRY carries, and that piece is a genuine driver: the
      ! Chapman-Cowling force is the gradient of the MOLE fraction x, and
      !
      !    logit(x) = logit(X) + ln(m_1/(m_He psi)),
      !
      ! so d logit(x)/dr = d logit(X)/dr - dln(psi)/dr.  Where the hydrogen
      ! turns molecular going down, each nucleus is spread over fewer
      ! collision partners, helium's mole fraction rises, and helium diffuses
      ! down that gradient -- outward across the molecular front -- even at a
      ! uniform mass fraction.  The term therefore rides in G as
      ! -dln(psi)/dr, which leaves the discretization, the Peclet hybrid and
      ! the M-matrix argument untouched and vanishes identically wherever psi
      ! is constant, i.e. over the whole atomic region.  The eddy term does
      ! not carry it: eddy mixing transports the mixture as a whole and has no
      ! preferred species, so it acts on dX/dr alone.
      !
      ! Nothing divides by a vanishing element density anywhere -- that
      ! division is what made the trace kernel's self-consistent coupling
      ! diverge.  The settling prefactor rho X (1-X) switches the settling off
      ! as either element is exhausted, while the gradient flux does not vanish
      ! there, so a helium-free cell next to a helium-bearing one still
      ! receives helium.
      !
      ! Sign check (memo 2.3): neutral gas at uniform composition has G > 0 and
      ! therefore J < 0 -- helium drifts inward, i.e. settles.
      !
      ! AMBIPOLAR FIELD.  eE is computed, not assumed:
      !
      !    eE = -(1/n_e) d(n_e k T)/dr = -k T dln(n_e T)/dr        (memo 3a)
      !
      ! (Koskinen et al. 2013, section 2.1), evaluated from the solved electron
      ! density and temperature of the cell.  The logarithmic form is used
      ! because it is exact for the barometric stratification of the base while
      ! being the same second-order central difference elsewhere.  Limits: a
      ! neutral gas gives Zbar_He - Zbar_1 = 0 (relative settling mass 3 m_H),
      ! an H+ plasma eE = m_H g/2 (2.5), a fully ionized He++ plasma
      ! eE = 4/3 m_H g (5/3).  he_ambipolar = .false. sets eE = 0.
      !
      ! DIFFUSION COEFFICIENT -- STAGE-RESOLVED (memo 2.6, decision D3).  The
      ! friction between the two elements is built pair by pair over their
      ! IONIZATION STAGES, because a helium ion moving through a proton gas is
      ! held by Coulomb collisions whose momentum-transfer cross section at
      ! T ~ 1e4 K is ~1e-12 cm^2 against the ~1e-15 cm^2 of a hard sphere.  A
      ! single neutral coefficient lets helium and the metals settle in the
      ! ionized wind by two to three orders of magnitude too fast.  Each pair
      ! coefficient is the Chapman-Enskog first approximation
      !
      !    D_st = 3 k T / (16 n mu_st Omega_st^(1,1))
      !
      ! with n the total CARRIER density (the Chapman-Enskog n is a density of
      ! colliding particles; in the atomic region it is the nucleus density,
      ! and where the hydrogen is molecular it is smaller by the nuclei bound
      ! into each molecule) and mu_st the reduced mass, evaluated
      ! in the limit the pair belongs to (routines below):
      !   neutral-neutral   hard_sphere_pair_diffusion    ~ T^(1/2)/n
      !   ion-neutral       ion_neutral_pair_diffusion    (polarization and
      !                                                    rigid core added)
      !   ion-ion           coulomb_pair_diffusion        ~ T^(5/2)/n
      ! and the element-element coefficient is the stage-fraction weighted
      ! harmonic sum (frictions add; Blanc's law across stages)
      !
      !    1/D_eff = sum_{s in He} sum_{t in comp.1} y_s y_t / D_st .
      !
      ! Limits: an all-neutral gas leaves one term, the Banks & Kockarts
      ! hard-sphere coefficient the earlier milestones used everywhere; a
      ! fully ionized gas leaves the He++/H+ Coulomb term alone (test T12).
      ! Resonant charge exchange H+ + H is NOT a term here: both partners
      ! carry the same element, and a binary diffusion coefficient is driven
      ! by the friction BETWEEN the elements, so every ion-neutral pair that
      ! appears is non-resonant.  HeH+ carries both elements and is left out
      ! of both carrier lists, and the metal nuclei are left out of n and of
      ! the reduced mass, both as trace approximations.
      !
      ! CARRIERS, AND THE MOLECULAR REGION (memo section 5).
      ! Equation (1) transports ELEMENTS, and the chemistry moves hydrogen
      ! between H, H+, H2, H2+, H3+ and helium between its stages without
      ! changing either element's mass, so the transport equation is the same
      ! above and below the molecular front.  What the species set changes is
      ! the CLOSURE, because the collision partners of helium are then a
      ! mixture of carriers:
      !
      !   friction   the stage sum above already runs over the carrier lists
      !              hcar_isp = {HI, HII, H2, H2+, H3+} and hecar_isp, each
      !              carrier entering with its own mass (an H2 is one partner
      !              of mass 2 m_H), charge and polarizability.  Blanc's law
      !              across carriers is the same harmonic sum as across
      !              stages, so it reduces to D(He,H) in the atomic region and
      !              to D(He,H2) in a fully molecular one, and interpolates
      !              between (test T7);
      !   forces     gravity and the ambipolar field act on the particles that
      !              actually collide, so G uses the mean carrier mass m_c1
      !              and mean carrier charge Zbar_c1 of the hydrogen component
      !              (per nucleus the weight of H2 equals that of H; per
      !              collision partner it is twice), and the trace-metal loop
      !              settles a metal against the same mean carrier mass;
      !   density    n in every pair coefficient is the carrier density.
      !
      ! The molecular columns of f_sp are zero unless the molecular chemistry
      ! is on, so no flag is tested anywhere: with it off every carrier list,
      ! mean carrier mass and carrier density collapses onto the nuclei and
      ! the atomic-region behavior is recovered identically.
      !
      ! DISCRETIZATION (memo section 3).  Finite volume on the existing grid,
      ! with ONE geometry for every term of one row: the face areas
      ! A(f) = r_edg(f)^2 and the exact shell volumes
      ! V(j) = (r_+^3 - r_-^3)/3 of spherical_face_area_and_cell_volume
      ! (grid_construction), which this module forms no copy of.  Both halves
      ! of the row -- the diffusive divergence and the divergence of the face
      ! mass flux -- are [A_+ F_+ - A_- F_-]/V_j on those faces, so an
      ! internal face cancels between the two cells that share it and the
      ! column sum of V_j times the row is the difference of the two boundary
      ! face fluxes.  A term weighted by r_j^2 (r_+ - r_-) instead carries the
      ! extra factor V_j/(r_j^2 dr_j), which differs between neighbours on a
      ! stretched grid, and the mixed operator is then the divergence of no
      ! single flux.
      !
      ! BOUNDARY NUCLEUS FLUX.  The column's two ends carry no diffusive
      ! flux: faces 0 and N are left at Agrd = Bdrf = 0, so an element
      ! crosses them only with the gas.  At the base the advective face flux
      ! is F_rho(1) times the reservoir composition wherever the face mass
      ! flux flows inward, which is what makes the Dirichlet cell a
      ! reservoir; at the outer face the ghost continues the interior,
      ! X_ghost = X_N, in BOTH directions of the face mass flux, so gas
      ! leaving carries the column's own composition and gas entering brings
      ! back the same, and no element flux is imposed at the outflow that the
      ! elemental face flux does not already carry.
      ! faces carrying r^2 areas, one implicit (backward-Euler) step per call
      ! with the transport COEFFICIENTS frozen within the step.  The step is
      ! nonlinear in X because the drift flux carries the product X(1-X), and
      ! it is solved by Newton on the residual with a tridiagonal Jacobian
      ! (solve_mass_fraction).  The gradient term is central; the drift term
      ! keeps the Peclet-based central/upwind hybrid: central where
      ! |B| dr <= 2 (A + E), donor-cell upwind otherwise.
      !
      ! THE DRIFT FLUX VANISHES AT BOTH ENDS OF THE COMPOSITION AXIS, IN THE
      ! DISCRETE OPERATOR AS IN THE CONTINUUM.  This is the property that
      ! bounds X, and it is why the two factors of X(1-X) are taken from
      ! OPPOSITE sides of the face.  The drift is a counter-flow: the helium
      ! mass flux -B X(1-X) is matched by an equal and opposite hydrogen flux,
      ! so a donor-cell rule has to take each element's mass fraction from the
      ! cell that element leaves.  With B >= 0 the helium drifts inward, so
      ! helium is donated by cell j+1 and hydrogen by cell j, and the face
      ! flux is
      !
      !    J_drift(f) = -B_f X(j+1) (1 - X(j))          [B_f >= 0, inward]
      !    J_drift(f) = -B_f X(j)   (1 - X(j+1))        [B_f <  0, outward]
      !
      ! Consequences, which the earlier lagged form did not have: a donor at
      ! X = 0 sends nothing and an ACCEPTOR at X = 1 receives nothing, so a
      ! cell sitting on either end of the axis cannot be pushed past it.  In
      ! the central branch the same holds because the Peclet condition makes
      ! the gradient flux dominate there: with X(j) = 1 the face value
      ! X_f (1 - X_f) = (1 + X(j+1))(1 - X(j+1))/4 <= (1 - X(j+1))/2, and
      ! |B_f| <= 2 A_f/dr bounds the drift flux by the gradient flux
      ! A_f (1 - X(j+1))/dr that carries helium OUT of that cell.
      !
      ! BOUNDS 0 <= X <= 1.  Let X be the solution of the implicit step and
      ! suppose it first touches 1 in cell m, every other cell still inside
      ! [0,1].  Then (i) the time term rho (X_m - X_m^old)/dt >= 0, (ii) the
      ! gradient flux leaves m at both faces because X_m is the maximum,
      ! (iii) the upwind advection contributes rho|v|(X_m - X_donor)/dr >= 0,
      ! and (iv) the drift flux is an outflow or zero at both faces by the
      ! rule above.  Every term of the row has the same sign and their sum is
      ! the row residual, which is zero -- so no term can be strictly
      ! positive, and X_m > 1 is impossible.  X <= 0 is the same statement:
      ! the scheme is exactly symmetric under X -> 1 - X, B -> -B (the
      ! donor/acceptor pair exchanges roles and the central branch is even in
      ! X_f - 1/2), which is the binary symmetry of the transport equation
      ! itself.  There is no cap f <= HeH and no base pile-up limiter: a
      ! pile-up, if one appears, is physics for the base boundary condition to
      ! answer, and it saturates at pure helium of its own accord.
      !
      ! The row sum is NOT rho/dt and the bounds do not need it to be: the
      ! argument above is made cell by cell on the residual, not from a
      ! comparison principle for the linear system.  (The lagged linearization
      ! this replaced entered the drift as [rho D G (1-X_lag)] X, whose flux
      ! does not vanish at X = 1 unless the lag is already there; measured, it
      ! let X reach 1.52.  Sections 84 and 86 of docs/Update_EXHALE_stage1.pdf.)
      !
      ! M-MATRIX CONDITION.  The Newton Jacobian has
      !   off-diagonals   aa = -K r^2 dJ/dX(j-1) <= 0,  cc = K r^2 dJ/dX(j+1) <= 0
      !   diagonal        bb = rho/dt + (outflow face terms) > 0
      ! for any iterate in [0,1], because (i) the gradient coefficient
      ! A + E >= 0 contributes +A/dr and -A/dr to the two slopes, (ii) the
      ! drift slopes are +B X or +|B|(1-X) on the left and -B(1-X) or -|B| X
      ! on the right in the two upwind branches, and never exceed A/dr in the
      ! central one (that is exactly the Peclet switch), and (iii) the upwind
      ! advection contributes +rho|v|/dr to the diagonal and the same amount,
      ! negated, to the donor neighbor and to nothing else.
      !
      ! What is measured, before the clip: the excursion outside [0,1] is
      ! NEGATIVE in every case tested -- the solve stays strictly inside the
      ! range and does not reach the clip at all.  -9.2e-4 on the strong-drift
      ! column of test T14, -8.6e-3 on the pure-helium band of T13, -1.05e-1
      ! on the LHS 1140 b wind; the lagged form overshot the same T14 column
      ! by +0.514.  So the clip is an assertion at both ends, and the
      ! excursions are measured every step and exposed
      ! (he_fraction_over_one / he_fraction_under_zero) rather than asserted.
      !
      ! Gated on he_diffusion (default .false.); when off the module is never
      ! entered, so flag-off runs are byte-identical to before.

      use global_parameters
      use grav_func,     only: Dphi
      ! The one spherical geometry of this grid: face areas r_edg^2 and the
      ! exact shell volumes.  Every divergence below divides by these, so
      ! the advective and the diffusive halves of one element row are the
      ! divergence of one flux (see DISCRETIZATION in the header).
      use grid_construction, only: spherical_face_area_and_cell_volume
      use species_table, only: n_bsp, bsp_fsp, bsp_nH, bsp_nHe,           &
                               bsp_is_excited_level, bsp_mass,            &
                               bsp_charge, isp_HI, isp_HeI,               &
                               isp_HII, isp_HeII, isp_HeIII, isp_HeTR,    &
                               isp_H2, isp_H2p, isp_H3p,                  &
                               bsp_nO, bsp_nC,                            &
                               n_mion, mion_fsp, mion_stage, mion_elem,   &
                               n_melem, melem_i0, melem_top, melem_A,      &
                               melem_name, iel_C, iel_O
      use composition,   only: mass_per_H_nucleus_without_He
      ! The rovibrational energy of H2 of the caloric equation of state, the
      ! one internal energy the sensible enthalpy of the components adds to
      ! (5/2) kT per particle (component_specific_enthalpies).
      use caloric_eos,   only: h2_rovibrational_energy_and_heat_capacity
      use species_advective_transport, only: species_advective_update,   &
                                    species_face_fraction,              &
                                    species_face_flux,                  &
                                    n_species_faces_bounded,            &
                                    species_face_excursion
      use lower_atmosphere_profile, only: lap_in_use, lap_r_top_RJ,       &
                          lap_flux_measured, lap_flux_window_empty,       &
                          lap_flux_nface, lap_FH_median, lap_FH_spread,   &
                          lap_FHe_median, lap_FHe_spread,                 &
                          lap_Mdot_median, lap_Mdot_spread,               &
                          lap_flux_r_lo_Rp, lap_flux_r_hi_Rp,             &
                          lap_flux_r_lo_measured,                         &
                          lap_steady_r_lo_Rp, lap_steady_nface,           &
                          lap_steady_FH_median,  lap_steady_FH_spread,    &
                          lap_steady_FHe_median, lap_steady_FHe_spread,   &
                          lap_steady_Mdot_median, lap_steady_Mdot_spread

      implicit none
      private
      public :: element_diffusion_step, relax_element_composition
      public :: interdiffusion_enthalpy_scale
      ! The advective half of the transport of every species row carried on
      ! the hydro's own face mass fluxes inside the Runge-Kutta stages: the
      ! element mass fractions and the declared molecular and proton
      ! carriers alike.
      public :: species_advection_active, species_advection_begin_step
      public :: species_advection_stage, species_advection_project
      ! The molecular and proton carriers of the photochemical transport
      ! operator ride on the same face mass fluxes as the elements.  That
      ! operator declares them here once the input keys are parsed, and reads
      ! their advected fractions back as the state one transport step starts
      ! from.
      public :: advected_carrier_reset, advected_carrier_register
      public :: advected_carrier_count
      public :: advected_carrier_fractions, advected_carrier_state_ready
      public :: advected_carrier_state_consumed
      ! The mass fractions the faces carry, and the mixture mass they are
      ! fractions of.  A stationary species balance forms its advective term
      ! from the SAME face quantities the Runge-Kutta stages do, so it reads
      ! them here instead of building a second copy of the conversion.
      public :: carrier_mass_fractions, element_mass_fractions
      public :: project_element_mass_fractions
      public :: mixture_mass_sum
      public :: element_transport_residual
      ! THE ELEMENT NUCLEUS FLUX THROUGH THE FACES, the single object the
      ! element row is the divergence of: the advective face element mass
      ! flux, the diffusive face flux, and the face element nucleus density
      ! and mass per nucleus a stage row is written against.  Exposed so
      ! that a caller reads the operator's own flux instead of rebuilding
      ! it, which would make an identity between two operators out of an
      ! identity within one.
      public :: element_nucleus_face_flux
      ! THE INTERDIFFUSION ENTHALPY FLUX of the energy equation of the
      ! diffusing mixture, the sensible enthalpy the element fluxes carry
      ! (Cook 2009, eqs. 11-13): its face value from given element fluxes,
      ! its divergence in the code units of the energy row, and both from a
      ! state in hand.  The marching update and the stationary energy row
      ! are written with these and nothing else.
      public :: interdiffusion_enthalpy_active
      public :: interdiffusion_enthalpy_face_flux
      public :: interdiffusion_enthalpy_divergence
      public :: interdiffusion_enthalpy_divergence_of_state
      public :: helium_diffusive_face_flux, trace_element_diffusive_face_flux
      public :: component_specific_enthalpies

      ! THE BUDGET ENTRY OF THIS OPERATOR'S BASE BOUNDARY CONDITION: the
      ! helium element mass flux [g cm^-2 s^-1] the base carries, ADVECTIVE
      ! AND DIFFUSIVE HALF SEPARATELY, positive OUTWARD (toward increasing
      ! r), at the two faces the boundary condition is written between.
      !
      !   the base face, r_edg(0), the contact with the reservoir and the
      !   level the lower atmosphere hands the column over at.  Its
      !   diffusive half is zero BY CONSTRUCTION and not by the state: the
      !   gradient and drift coefficients of faces 0 and N are left at zero
      !   (drift_and_gradient_face_coefficients), so an element crosses
      !   either end of the column with the gas alone.  It is carried in the
      !   record as its own number rather than left to be assumed.
      !
      !   the first solved face, r_edg(first_solved_face), which bounds the
      !   region the operator actually solves (rows 2..N against the
      !   Dirichlet reservoir, so face 1; face 0 when a closed base makes
      !   cell 1 an unknown).  This is the face through which the reservoir
      !   composition feeds the solved column, and the diffusive half here
      !   is the flux the Dirichlet condition drives.
      !
      ! THE ADVECTIVE HALF RIDES ON THE RIEMANN FACE MASS FLUX the element
      ! transport of the same evaluation rode on: the flux the caller handed
      ! the step (the fixed-wind relaxation), or on the marching path the
      ! face mass flux of the step the Runge-Kutta stages moved the mass and
      ! the elements with (Frho_step below).  An evaluation that had neither
      ! -- a direct call with no stages before it, which only the acceptance
      ! tests make -- has no face mass flux to state, and advective_measured
      ! says so instead of a zero standing in for the flux.
      !
      ! A REPORT READS THE EVALUATION THE NUMBERS BELONG TO.  evaluation
      ! counts the operator evaluations that produced a composition, so a
      ! caller that notes the count before its own step can tell this
      ! record from an older one, and evaluation = 0 means NO OPERATOR HAS
      ! RUN, which is not a flux of zero.  element_base_flux_report writes
      ! the record with that distinction made.
      type, public :: base_element_flux_record
         integer :: evaluation = 0
         integer :: first_solved_face = 0
         logical :: advective_measured = .false.
         real*8  :: advective_base_face = 0.0d0
         real*8  :: diffusive_base_face = 0.0d0
         real*8  :: advective_solved_face = 0.0d0
         real*8  :: diffusive_solved_face = 0.0d0
      end type base_element_flux_record
      type(base_element_flux_record), public, save :: element_base_flux
      public :: element_base_flux_report
      ! The conversion of a face mass flux from the code units of the mass
      ! row to g cm^-2 s^-1, the one rule every face record of this module
      ! states its fluxes in.
      public :: face_mass_flux_cgs

      ! THE ELEMENTAL FACE FLUXES OF A STATE, as a radial profile
      ! (output/element_flux_profile.txt) and, with a lower-atmosphere
      ! profile in use, reduced over the two windows the flux closure reads.
      ! A statement about ONE state, written by the evaluation that measures
      ! that state (the certification) on that state's own Riemann face mass
      ! flux, and not by a transport step.
      public :: write_element_flux_profile, element_flux_profile_on

      ! THE HELIUM ELEMENT ROW OF A STATE, TERM BY TERM, as the stationary
      ! elemental transport balance forms it (element_transport_residual),
      ! so that a row the certification reports can be read by the terms it
      ! balances rather than inferred from its measure.  Filled only when a
      ! caller asks for it (the he_terms argument of
      ! element_transport_residual), written by element_row_terms_write on
      ! request (EXHALE_ELEMENT_ROW_TERMS=1, default off).  Nothing in the
      ! solution reads it.
      !
      ! Cells 1..N, every rate in g cm^-3 s^-1:
      !   dif_in, dif_out  the diffusive (gradient, eddy and settling drift)
      !                    face flux of the inner face j-1 and of the outer
      !                    face j, each with its face area and divided by the
      !                    cell volume, signed as they enter the row, so
      !                    dif_in + dif_out is the diffusive divergence;
      !   adv_in, adv_out  the same for the advective face element mass flux
      !                    F_rho Y_face, converted with the row's own factor
      !                    n0 mu msum(j) v0/R0;
      !   res, scale       the row and the scale the certification divides it
      !                    by (the sum of the magnitudes of the row's own
      !                    terms, floored at floor);
      !   rho_phys, X      the mass density and the helium mass fraction of
      !                    the cell, D12 and Kzz its binary and eddy
      !                    coefficients [cm^2 s^-1].
      ! Faces 0..N: F_rho the Riemann face mass flux and J the diffusive
      ! helium flux [g cm^-2 s^-1], Y_face the face helium mass fraction.
      type, public :: element_row_terms_record
         real*8, allocatable :: rho_phys(:), X(:), D12(:), Kzz(:)
         real*8, allocatable :: dif_in(:), dif_out(:)
         real*8, allocatable :: adv_in(:), adv_out(:)
         real*8, allocatable :: res(:), scale(:), floor(:)
         real*8, allocatable :: F_rho(:), Y_face(:), J(:)
      end type element_row_terms_record
      public :: element_row_terms_on, element_row_terms_write
      ! The nucleus counts per unit mass of the two elements, from the one
      ! stoichiometric map of the species table.  A row written per nucleus
      ! of an element divides by that element's nucleus DENSITY, and the
      ! face density the flux above returns is the arithmetic mean of that
      ! same cell quantity, so the cell and the face value have to come
      ! from one count.
      public :: element_nucleus_counts
      ! The elemental transport balance of a state reduced to ONE number,
      ! and the admissibility test a residual norm has to meet before a
      ! progress control may read it.
      public :: element_transport_residual_norm
      public :: residual_norm_is_admissible
      ! The floor a row measure divides by, so that a cell in which every
      ! term of a row vanishes reads zero instead of dividing by zero.  It
      ! is the value certification_row_measure uses (cert_scale_floor),
      ! repeated here because the certification stands above this module and
      ! the two must give one number for one row.
      real*8, parameter, public :: element_row_scale_floor = 1.0d-300
      ! Whether a species vector is a vector of ordinary reals: the outer
      ! iteration asks it of the composition a relaxation hands back.
      public :: every_species_is_finite
      ! THE TERMS OF ONE TRACE-ELEMENT ROW AT THE OUTERMOST CELLS, printed
      ! and nothing else: a row that stands away from zero is read by the
      ! sizes of the two fluxes it balances, and those sizes exist only
      ! inside the residual that forms them.  Off unless a caller asks, and
      ! the caller states the first cell it wants (PLAN_20260909_rev1 item
      ! N25, the sodium row of the outermost cell).  The arithmetic of the
      ! row is untouched by it.
      public :: trace_row_terms_diag, trace_row_terms_from
      public :: relative_settling_mass
      ! Exposed so the acceptance tests read the same coefficients the
      ! operator uses -- there is no second copy of the friction anywhere.
      public :: helium_hydrogen_diffusion
      public :: hard_sphere_pair_diffusion, polarization_pair_diffusion
      public :: ion_neutral_pair_diffusion
      public :: coulomb_pair_diffusion, coulomb_logarithm
      public :: alpha_HI, alpha_HeI
      ! Exposed so an acceptance test reads the masses this operator closes
      ! its mixture with, and not a second copy of them.
      public :: m_He_amu, m_H_amu
      ! How far the last solve left [0,1] BEFORE the range clip, signed so
      ! that a negative value means it stayed inside.  Exposed because the
      ! clipped X cannot tell an overshoot from an exact 1 (both read 1.0), so
      ! an acceptance test that means to check the bounds has to read the
      ! solve and not its clip.
      public :: he_fraction_over_one, he_fraction_under_zero
      real*8, protected :: he_fraction_over_one   = -1.0d0
      real*8, protected :: he_fraction_under_zero = -1.0d0
      ! Newton iterations and the residual the last solve stopped at, in units
      ! of X (the residual is scaled by dt/rho).  Reported by report_step.
      integer, protected :: he_fraction_newton_steps = 0
      real*8,  protected :: he_fraction_newton_resid = 0.0d0
      public :: he_fraction_newton_steps, he_fraction_newton_resid
      ! The same measurement for the trace-metal mixing ratio, whose clip
      ! carried the same unmeasured assertion the helium one did.  Reset at
      ! the start of each step and maximized over the elements.
      real*8, protected :: trace_ratio_under_zero = -1.0d0
      public :: trace_ratio_under_zero

      ! WHETHER THE COMPOSITION ONE TRANSPORT STEP HANDS BACK IS ADMISSIBLE.
      ! The step solves a discrete transport equation at a FIXED density, so
      ! the composition it returns is admissible only if it still carries
      ! that density as well as the one it was handed did: sum_i m_i n_i =
      ! rho, the closure calc_rho reads out of the species vector.  Finite
      ! and nonnegative populations are necessary and not sufficient, and a
      ! helium fraction clipped into [0,1] says nothing about the row it was
      ! meant to zero.  The named outcomes below are what a caller reads
      ! instead of the diagnostics above, which measure and do not judge.
      integer, parameter :: element_step_accepted            = 0
      integer, parameter :: element_step_solve_failed        = 1
      integer, parameter :: element_step_nonfinite           = 2
      integer, parameter :: element_step_out_of_bounds       = 3
      integer, parameter :: element_step_mass_closure_failed = 4
      public :: element_step_accepted, element_step_solve_failed
      public :: element_step_nonfinite, element_step_out_of_bounds
      public :: element_step_mass_closure_failed
      ! THE OUTCOME OF THE MOST RECENT CALL, whether or not that call asked
      ! for status.  The candidate a step produces is judged against the
      ! same four tests every time: acceptance of a numerical update cannot
      ! depend on whether the caller asked to be told the result.  A caller
      ! that does not read status still gets the composition it was handed
      ! back unmoved on a refusal, and can read this variable afterward to
      ! find out that happened (the marching path at EXHALE_main.f90 reads
      ! it to decide whether to shorten the step).
      integer, protected :: element_step_last_status = element_step_accepted
      public :: element_step_last_status
      ! TEST KNOB: the next this many calls report the nonlinear solve
      ! unsolved whatever it did, so a caller's handling of a refused
      ! step can be exercised on a column that would not refuse by
      ! itself (EXHALE_ELEMENT_REFUSE_STEPS, read by the marching loop;
      ! zero in every production run).
      integer, public :: element_refuse_leading_steps_for_test = 0
      ! TEST KNOB: with a value n >= 0, n steps are let through and every
      ! later one is reported unsolved, so a relaxation that keeps steps and
      ! then runs out of retries can be exercised; -1 (every production
      ! run) leaves it off.
      integer, public :: element_refuse_after_accepted_for_test = -1
      public :: element_step_outcome_text
      ! THE OUTCOME OF A RELAXATION, which is not the outcome of one step.
      ! A relaxation whose composition stopped moving has reached the fixed
      ! point of the transport operator; one that ran out of steps has not,
      ! and its returned state is a partial advance whatever its movement
      ! over the last step was.  A relaxation that could take no admissible
      ! step at all returns the composition it was given.
      integer, parameter :: element_relaxation_converged   = 0
      integer, parameter :: element_relaxation_step_budget = 1
      integer, parameter :: element_relaxation_failed      = 2
      public :: element_relaxation_converged
      public :: element_relaxation_step_budget
      public :: element_relaxation_failed
      ! max_j |sum_i m_i n_i - rho|/rho over the physical cells, in the mass
      ! policy of mixture_mass_split (which is calc_rho's own, term for
      ! term).  It is the departure of the composition the last judged step
      ! or relaxation PRODUCED: on a refusal that is the candidate that was
      ! refused and not the state handed back, which is what says why the
      ! refusal happened.  Negative where no candidate was produced (a step
      ! whose nonlinear solve failed) and until the first judged call.
      real*8, protected :: element_mass_closure_departure = -1.0d0
      public :: element_mass_closure_departure
      ! HOW MUCH FARTHER FROM THE DENSITY THE SPECIES MAY END THAN THEY
      ! BEGAN.  The transport moves element amounts and the projection meets
      ! both element totals exactly while closing the mixture around the mass
      ! the entry species carried, so the departure from rho is an invariant
      ! of the step and what is left in it is the round-off of that
      ! inversion: 1.0e-15 at the entry of the atomic element fixture and
      ! 1.0e-14 after a relaxation of the synthetic columns (both MEASURED).
      ! Over the he_relax_maxstep
      ! steps of a relaxation those accumulate at worst to ~4e-13, so this
      ! bound stands two and a half decades above the arithmetic and ten
      ! decades below a mass error with any physical meaning.  It is stated
      ! as an increase and not as an absolute departure because establishing
      ! the closure belongs to whoever builds the composition (the ionization
      ! solve normalizes it: 1.0e-15 MEASURED on the atomic fixture) and
      ! preserving it belongs to the transport; on an admissible entry state
      ! the two readings are the same number.
      real*8, parameter :: element_mass_closure_tol = 1.0d-10
      ! The drift flux vanishes at both ends of the composition axis in the
      ! discrete operator as in the continuum (module header), so a converged
      ! step leaves X inside [0,1] up to the round-off of the tridiagonal
      ! solve.
      real*8, parameter :: element_fraction_bound_tol = 1.0d-12

      ! THE ADVECTED ELEMENT MASS FRACTIONS, between the beginning of a
      ! marching step and the projection that ends it.  Column 1 is helium,
      ! the member of the normalized set; columns 1+im are the trace metals,
      ! which ride on the same faces without entering the normalization
      ! because their mass is already inside m_1 at the fixed reservoir
      ! abundance the trace closure assumes (the same closure that lets a
      ! metal element be solved against a frozen hydrogen background).
      ! Yetr_n is the composition the step began with, which the second and
      ! third Runge-Kutta stages combine with.
      real*8, allocatable :: Yetr(:,:), Yetr_n(:,:)
      integer :: n_etr = 0

      ! THE FACE MASS FLUX OF ONE MARCHING STEP.  The three stages of the
      ! SSP-RK3 update advance the mass row by
      !
      !    rho^(n+1) = rho^n - dt [ D F^(0)/6 + D F^(1)/6 + 2 D F^(2)/3 ]
      !
      ! with D the spherical divergence and F^(k) the Riemann face mass flux
      ! of stage k (species_advective_update combines the element rows with
      ! the same weights), so Frho_step = F^(0)/6 + F^(1)/6 + 2 F^(2)/3 is
      ! the face mass flux whose divergence is the density change of the
      ! step, and at a stationary state, where every stage returns the flux
      ! of the state, it is that flux.  species_advection_stage accumulates
      ! it; the element step of the marching path, which is handed no face
      ! mass flux of its own, reads it once for its base budget entry.  It
      ! exists from the third stage of an attempt to that element step:
      ! species_advection_begin_step and the reading both clear
      ! step_flux_ready, so no attempt reads another's flux.
      real*8, allocatable :: Frho_step(:)
      logical :: step_flux_ready = .false.

      ! WHETHER THE TERMS OF A TRACE-ELEMENT ROW ARE PRINTED, and the first
      ! cell they are printed from (trace_composition_residual).  Read only
      ! by a print; the row and its boundary are the same either way.
      logical :: trace_row_terms_diag = .false.
      integer :: trace_row_terms_from = 0

      ! THE ADVECTED CARRIER MASS FRACTIONS.  A carrier is one species of an
      ! element, so it rides on the same faces WITHOUT entering the
      ! normalized set: the element it belongs to is already carried by the
      ! closing member, and normalizing the carrier beside it would count
      ! that mass twice.  What the face identity states for a carrier is
      ! therefore a statement about its element, and it is the element
      ! totals, not the carriers, that the closing member closes.
      !
      ! THE INNER GHOSTS ARE THE INFLOW COMPOSITION.  Where a lower-
      ! atmosphere handoff states a carrier's base partition the ghost holds
      ! that value, pinned by the ionization sweep, and the base face carries
      ! it in; where nothing states it the ghost takes the base cell's own
      ! partition, which is the zero-gradient condition written on the face
      ! the gas actually crosses.  Which of the two a carrier gets is
      ! declared with it.
      integer, parameter :: n_car_max = 8
      integer :: car_isp(n_car_max)  = 0
      real*8  :: car_mass(n_car_max) = 0.0d0
      logical :: car_base_imposed(n_car_max) = .false.
      integer :: n_car = 0
      real*8, allocatable :: Ycar(:,:), Ycar_n(:,:)
      ! Whether Ycar holds a composition this attempt advected.  False from
      ! the top of every attempt until the projection, so an attempt that is
      ! discarded before the projection leaves no advected state behind.
      logical :: car_ready = .false.

      ! Species masses in the m_H units the code counts f_sp in, READ FROM
      ! THE SPECIES TABLE rather than restated, because the two-component
      ! closure below (msum = rho_1 + m_He_amu*nucHe = rho/(n0 mu)) is exact
      ! only while this helium mass is the one calc_rho weighs helium with.
      ! Indexing note: bsp_mass is indexed by bsp POSITION and bsp_fsp(1:6)
      ! = 1..6, so isp_HI/isp_HeI are also the positions of the two atomic
      ! rows; the shortcut does NOT extend to the molecular species
      ! (isp_H2 = 34 but bsp position 7).  The gram value of the H = 1 unit
      ! is the hydrogen ATOM as in parameters.f90, not the atomic mass unit u.
      real*8, parameter :: m_He_amu = bsp_mass(isp_HeI)
      real*8, parameter :: m_H_amu  = bsp_mass(isp_HI)
      ! The mass unit of the density normalization is the hydrogen ATOM mass mu
      ! of parameters.f90 (one definition); this module reads it and keeps no copy.
      ! Elementary charge in electrostatic units, from the CODATA 2018 exact
      ! coulomb value: e = 1.602176634e-19 C / (10/c) = 4.803204713e-10 esu.
      real*8, parameter :: e_esu    = 4.803204713d-10
      ! Hard-sphere prefactor of Banks & Kockarts (1973): the neutral binary
      ! coefficient is D = 1.52e18 (1/A_s + 1/A_t)^(1/2) sqrt(T)/n.
      real*8, parameter :: bk_hs_pref = 1.52d18

      ! Static dipole polarizabilities of the NEUTRAL collision partners, the
      ! only atomic datum the ion-neutral (polarization) coefficient needs.
      ! Quoted in atomic units a_0^3 and converted with
      ! a_0^3 = 1.481847e-25 cm^3.  alpha(H) = 4.5 a_0^3 is the exact
      ! nonrelativistic ground-state value; the rest are the recommended
      ! static dipole polarizabilities of Schwerdtfeger & Nagle (2019,
      ! Mol. Phys. 117, 1200).  As a cross-check H and He come out at 0.667
      ! and 0.205 Angstrom^3, the values tabulated for aeronomy by
      ! Schunk & Nagy (Ionospheres, Table 4.1).
      real*8, parameter :: a0cub_cm3 = 1.481847d-25
      real*8, parameter :: alpha_HI  =   4.500d0*a0cub_cm3
      real*8, parameter :: alpha_HeI =   1.383d0*a0cub_cm3
      real*8, parameter :: alpha_H2  =   5.315d0*a0cub_cm3
      ! neutral metals, in the melem order C O N Mg Si Ca Na K S Fe
      real*8, parameter :: alpha_melem(n_melem) = a0cub_cm3*                &
           [  11.3d0,   5.3d0,   7.4d0,  71.2d0,  37.3d0,                   &
             160.8d0, 162.7d0, 289.7d0,  19.4d0,  62.0d0 ]

      ! Carriers of component 1 (hydrogen).  H2/H2+/H3+ are zero unless the
      ! molecular chemistry is on, so no flag has to be tested; with it off
      ! the list collapses to HI/HII.  Masses are per CARRIER (an H2 is one
      ! collision partner of mass 2), which is what the friction counts.
      integer, parameter :: n_hcar = 5
      integer, parameter :: hcar_isp(n_hcar) =                              &
           [ isp_HI, isp_HII, isp_H2, isp_H2p, isp_H3p ]
      real*8,  parameter :: hcar_Z(n_hcar) =                                &
           [ 0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0 ]
      real*8,  parameter :: hcar_m(n_hcar) =                                &
           [ 1.0d0, 1.0d0, 2.0d0, 2.0d0, 3.0d0 ]
      real*8,  parameter :: hcar_alpha(n_hcar) =                            &
           [ alpha_HI, 0.0d0, alpha_H2, 0.0d0, 0.0d0 ]
      character(len=5), parameter :: hcar_name(n_hcar) =                    &
           [ 'HI   ', 'HII  ', 'H2   ', 'H2+  ', 'H3+  ' ]

      ! Carriers of the helium component.  The 2^3S triplet is a neutral
      ! helium atom for the friction -- same mass, same charge, same
      ! polarizability as a ground-state He I atom -- and it is an EXCITED
      ! LEVEL of He I, so the HeI column already counts it
      ! (bsp_is_excited_level).  It is therefore not a list entry of its own:
      ! a separate HeTR row made every triplet atom two collision partners.
      integer, parameter :: n_hecar = 3
      integer, parameter :: hecar_isp(n_hecar) =                            &
           [ isp_HeI, isp_HeII, isp_HeIII ]
      real*8,  parameter :: hecar_Z(n_hecar) =                              &
           [ 0.0d0, 1.0d0, 2.0d0 ]
      ! One helium nucleus with its electrons weighs the same in every
      ! ionization stage to the digits the mass table keeps, so the three
      ! rows are the single table mass and not three literals.
      real*8,  parameter :: hecar_m(n_hecar) =                              &
           [ m_He_amu, m_He_amu, m_He_amu ]
      real*8,  parameter :: hecar_alpha(n_hecar) =                          &
           [ alpha_HeI, 0.0d0, 0.0d0 ]
      character(len=5), parameter :: hecar_name(n_hecar) =                  &
           [ 'HeI  ', 'HeII ', 'HeIII' ]

      ! Largest number of ionization stages carried for one metal element.
      integer, parameter :: n_mstage = 3
      ! Composition relaxation at a fixed wind (relax_element_composition):
      ! step budget, and the absolute convergence measure max|dX|/X_base.
      integer, parameter :: he_relax_maxstep = 400
      real*8,  parameter :: he_relax_tol     = 1.0d-12
      ! A step the operator refuses is discarded and retried at half its
      ! length.  The retry is bounded twice: by the number of discarded steps
      ! and by the step length itself, which at 1e-6 of the shortest
      ! composition time scale of the grid moves the composition by less than
      ! the movement measure he_relax_tol can see, so a shorter one carries
      ! no advance whether it solves or not.
      integer, parameter :: he_relax_retry_max  = 20
      real*8,  parameter :: he_relax_grow_floor = 1.0d-6

      contains

      ! ------------------------------------------------------------------ !

      subroutine element_diffusion_step(rho, Tcode, f_sp, dt_code,        &
                                        closed_base, Jface_out, Frho_in,  &
                                        status, JXface_out)
      ! Advance the helium mass fraction X one relaxation step and project the
      ! new element totals back into f_sp.  rho, Tcode are the current
      ! adimensional primitives, dt_code the adimensional relaxation timestep;
      ! f_sp is modified in place.  No cell velocity is taken: the material
      ! transport is the divergence of the FACE mass fluxes (below).
      !
      ! closed_base (optional, default .false.) replaces the Dirichlet
      ! reservoir base by a zero-flux inner boundary, which is what the closed
      ! column of tests T1a/T4/T6 needs; production runs never set it.
      ! Jface_out (optional) returns the diffusive helium mass flux
      ! [g cm^-2 s^-1] at the faces r_edg(0:N) evaluated with the coefficients
      ! the step actually used and the NEW X, so that a discrete elemental
      ! budget closes exactly against it (test T1b).
      ! JXface_out (optional) returns, in the same units and at the same
      ! faces, the diffusive MASS flux of every trace metal relative to the
      ! hydrogen it diffuses through, formed with the coefficients the metal
      ! step used and its new mixing ratio (trace_face_coefficients); zero
      ! unless He_metal_diffusion moves the metals.  Both are zero on a call
      ! that moves nothing (flag off, no helium, or an unsolved step), and a
      ! caller uses them only for a step whose composition was accepted.
      ! Frho_in (optional) is the FACE MASS FLUX F_rho(j) at r_edg(j) in code
      ! units, the one the Riemann solve of this state returned and the mass
      ! row of this state differences, and its presence is what turns the
      ! advective term on.  The marching path omits it: there the advection is
      ! the divergence of those same face fluxes taken inside the Runge-Kutta
      ! stages (species_advection_stage), so repeating it here would advect
      ! the elements twice.  The relaxation at a fixed wind has no stages to
      ! ride on and supplies the flux itself.  The base budget entry of the
      ! step is written on the face mass flux the element transport of THIS
      ! evaluation rode on: Frho_in where it is given, otherwise the step's
      ! flux of the Runge-Kutta stages that preceded the call (Frho_step).
      !
      ! THE STEP HANDS BACK A NEW COMPOSITION ONLY IF IT IS ADMISSIBLE, EVERY
      ! CALL, WHETHER OR NOT status IS ASKED FOR.  A composition that does
      ! not reconstruct the density it was advanced under, or whose helium
      ! mass fraction left [0,1], is not a state of the gas whichever loop
      ! asked for the step: acceptance of a numerical update cannot depend
      ! on a diagnostic argument.  The candidate is accepted only if the
      ! nonlinear solve reported the step solved and the result is finite
      ! everywhere, X inside [0,1] to element_fraction_bound_tol before the
      ! range clip, and standing no farther from the density the step held
      ! fixed, by more than element_mass_closure_tol, than the composition
      ! it was handed did; otherwise the entry composition is restored.  The
      ! outcome is always left in element_step_last_status, and in status
      ! too where the caller supplied it.

      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: dt_code
      logical, optional,                      intent(in)    :: closed_base
      real*8, dimension(0:N), optional,       intent(out)   :: Jface_out
      real*8, dimension(1-Ng:N+Ng), optional, intent(in)    :: Frho_in
      integer, optional,                      intent(out)   :: status
      real*8, dimension(0:N,n_melem), optional, intent(out) :: JXface_out

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, mass1, Xhe, Xold
      real*8, dimension(1-Ng:N+Ng) :: msum_out
      real*8, dimension(1-Ng:N+Ng) :: carH, carHe, mcarH, dlnpsi
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, TK, Dco, Gco, dt_phys
      real*8, dimension(1-Ng:N+Ng) :: rp, nH_phys, ntot_phys, dmeff
      real*8, dimension(1-Ng:N+Ng) :: nX, nXold, DcoX, GcoX, zb1, Frho
      real*8, dimension(1-Ng:N+Ng) :: cadvf, wYtr, cadvX
      real*8, dimension(1-Ng:N+Ng) :: eEf, ne_phys, zbHe, ne_rel
      real*8, dimension(1-Ng:N+Ng,n_hcar)    :: yH
      real*8, dimension(0:N)       :: Agrd, Bdrf, Jf, PLX, PRX
      real*8, dimension(1-Ng:N+Ng) :: fXnew
      integer, dimension(0:N)      :: updrf
      ! The face mass flux the element transport of this evaluation rode on
      ! (code units, and in g cm^-2 s^-1), whether one exists, and the face
      ! helium mass fraction it carries: the advective half of the base
      ! budget entry below.
      real*8, dimension(1-Ng:N+Ng) :: Frho_eval, Frho_face, Yface
      logical :: have_face_flux
      real*8  :: sv_face_excursion
      integer :: sv_faces_bounded
      ! NB: local scalars are checked against global_parameters case-
      ! insensitively.  In particular the time scale is tscale, NOT t0 (a local
      ! t0 would alias the global temperature normalization T0), and nothing
      ! here is named N, Ng, r, g, mu, info, count or du.
      real*8, dimension(:,:), allocatable :: f_entry
      real*8 :: tscale, X_base, m_1, fXbase, rXsc, Xover, Xunder, qdep
      real*8 :: dJl, dJr, closure
      integer :: j, jlo, im, i0m, top, k, n_vanished, outcome
      logical :: shut_base, advect, solved

      if (present(status)) status = element_step_accepted
      element_step_last_status = element_step_accepted
      if (present(Jface_out))  Jface_out  = 0.0d0
      if (present(JXface_out)) JXface_out = 0.0d0
      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      ! The composition to restore if the candidate turns out inadmissible.
      ! Kept on every call, not only where a caller reads status: whether
      ! the update is a state of the gas does not depend on who is asking.
      allocate(f_entry(1-Ng:N+Ng,n_species))
      f_entry = f_sp

      shut_base = .false.
      if (present(closed_base)) shut_base = closed_base
      jlo = 2
      if (shut_base) jlo = 1

      m_1 = mass_per_H_nucleus_without_He()

      ! --- element nucleus counts per unit mass, and the helium mass fraction
      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      Xhe  = m_He_amu*nucHe/msum
      Xold = Xhe

      ! --- dimensional fields
      TK       = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      rp       = r*R0
      rho_phys = rho*n0*mu*msum                       ! [g/cm^3]
      ! Carrier (collision-partner) density, which is what the Chapman-Enskog
      ! coefficients divide by.  Equal to the nucleus density in the atomic
      ! region; smaller where the hydrogen is bound into molecules.
      call carrier_counts(f_sp, carH, carHe, mcarH)
      ntot_phys = (carH + carHe)*rho*n0                    ! [cm^-3] carriers
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      ! D_12 [cm^2/s], stage-resolved (module header / memo 2.6)
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dco)
      tscale  = R0/v0
      dt_phys = dt_code*tscale
      where (dt_phys .lt. 1.0d-30) dt_phys = 1.0d-30

      ! --- face mass flux, and with it the advective term.  Absent, the
      !     Runge-Kutta stages already carried it.  cadvf converts the code
      !     divergence of the helium mass flux into the [g cm^-3 s^-1] the
      !     row is written in: a code density times a code velocity over a
      !     code length is a rate in units of v0/R0, and n0 mu msum is the
      !     mass one unit of the code density carries.
      advect = present(Frho_in)
      if (advect) then
         Frho = Frho_in
      else
         Frho = 0.0d0
      endif
      cadvf = n0*mu*msum*v0/R0

      ! --- the face mass flux this evaluation's element transport rode on,
      !     for the base budget entry: the one handed in, or the step's flux
      !     of the stages that preceded this call, read once.  Taken here,
      !     before the solve, so that an unsolved step consumes it as well
      !     and no later call can read a flux that belongs to this attempt.
      have_face_flux = .true.
      if (advect) then
         Frho_eval = Frho_in
      else if (step_flux_ready) then
         Frho_eval = Frho_step
         step_flux_ready = .false.
      else
         Frho_eval = 0.0d0
         have_face_flux = .false.
      endif

      ! --- settling coefficient G [1/cm] (gravity + ambipolar field + thermal)
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1, eEf,   &
                                dlnpsi)

      ! --- THE BASE BOUNDARY CONDITION OF THIS OPERATOR, stated here and
      ! independent of what the advective boundary does with the contact.
      !
      ! The elemental composition of the lower atmosphere at the base level
      ! is prescribed data (the He/H of the input, or of the lower-atmosphere
      ! profile at the matching level), so the base cell and the inner ghosts
      ! carry it as a DIRICHLET value and the operator solves rows 2..N
      ! against it.  It is not an inflow condition: a Dirichlet value drives
      ! a diffusive flux through the base face whichever way the bulk gas
      ! moves and whether or not it moves at all, and the flux it drives is
      ! recorded below as this boundary condition's budget entry.  The
      ! closed-column tests set closed_base and get a zero-flux inner
      ! boundary instead; production runs never do.
      X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)

      ! --- implicit step.  The face coefficients do not depend on X (the
      ! composition enters only through the drift product X(1-X), which the
      ! solve carries at the new level), so there is no Picard sweep to make:
      ! one Newton solve of the nonlinear step is the whole update.
      call drift_and_gradient_face_coefficients(rho_phys, Dco, Gco, rp,   &
                                                Agrd, Bdrf, updrf)
      call solve_mass_fraction(Xold, Xhe, rho_phys, dt_phys, Frho,        &
                               cadvf, Agrd, Bdrf, updrf, X_base, jlo,      &
                               advect, solved)

      ! Range clip.  Both lines are assertions: the drift flux vanishes at
      ! both ends of the composition axis in the discrete operator as in the
      ! continuum (header), so only round-off can leave [0,1].  How far the
      ! solve left it is measured and not asserted -- the clipped X cannot
      ! tell an assertion from a limiter, because both read 1.0 -- and the
      ! excursions are kept in module variables the acceptance tests read.
      ! NB a solve returning X = 1 to the last bit needs no clip at all
      ! (1 - 1e-17 rounds to 1 in double), so X = 1 in the report is not by
      ! itself evidence that these lines fired.
      Xover  = maxval(Xhe(1:N)) - 1.0d0
      Xunder = -minval(Xhe(1:N))
      he_fraction_over_one   = Xover
      he_fraction_under_zero = Xunder

      ! An unsolved step has no composition to project: the row it was meant
      ! to zero is not zeroed and the bounds on X are the converged step's.
      ! The entry composition stands, and the trace metals are not advanced
      ! against a hydrogen background that was never moved.
      if (element_refuse_leading_steps_for_test .gt. 0) then
         element_refuse_leading_steps_for_test =                           &
            element_refuse_leading_steps_for_test - 1
         solved = .false.
      endif
      if (element_refuse_after_accepted_for_test .eq. 0) then
         solved = .false.
      else if (element_refuse_after_accepted_for_test .gt. 0 .and. solved) &
         then
         element_refuse_after_accepted_for_test =                          &
            element_refuse_after_accepted_for_test - 1
      endif
      if (.not. solved) then
         f_sp   = f_entry
         element_step_last_status = element_step_solve_failed
         if (present(status)) status = element_step_solve_failed
         element_mass_closure_departure = -1.0d0     ! no candidate to weigh
         deallocate(f_entry)
         return
      endif

      where (Xhe .lt. 0.0d0) Xhe = 0.0d0
      where (Xhe .gt. 1.0d0) Xhe = 1.0d0
      if (.not. shut_base) Xhe(1-Ng:1) = X_base    ! base + inner ghosts
      Xhe(N+1:N+Ng) = Xhe(N)                       ! zero-gradient outer ghost

      ! --- diffusive face flux actually carried by the step, with the
      ! coefficients of the matrix that was solved and the new X (test T1b).
      ! Face 0 is always evaluated: it is the flux the base boundary
      ! condition above drives, the budget entry that says how much helium
      ! the lower atmosphere put through the level in this step, and a
      ! quantity a zero bulk velocity does not make zero.
      do j = 0, N
         call element_face_flux(Xhe(j), Xhe(j+1), Agrd(j), Bdrf(j),       &
                                updrf(j), Jf(j), dJl, dJr)
      enddo
      ! --- the budget entry of the base boundary condition, both halves at
      ! both of the faces the condition is written between (declaration
      ! above), in g cm^-2 s^-1 and positive outward.  The advective half is
      ! the Riemann face mass flux this evaluation's element transport rode
      ! on (Frho_eval above) times the reconstructed upwinded face helium
      ! mass fraction of the new composition, the face composition the
      ! stages and the relaxation put on the same face.  The face mean of
      ! the cell-centred rho and v is not that flux and is not used: near
      ! the base the cell-centred velocity carries a collocated two-cell
      ! odd-even mode that the Riemann face flux does not carry, and at
      ! face 1 of the certified LHS 1140 b kzz1e9 He/H 0.55 state the face
      ! mean reads -1.07e9 g/s (4 pi r^2 times the flux) against the wind's
      ! 3.93e7 g/s that the Riemann flux carries at every face (MEASURED
      ! 2026-09-24).
      ! species_face_fraction counts the face states it scales back onto
      ! [0,1], and those counters are statements about the trajectory, so a
      ! budget entry evaluated on a state already in hand must not add to
      ! them.
      element_base_flux%evaluation = element_base_flux%evaluation + 1
      element_base_flux%first_solved_face  = jlo - 1
      element_base_flux%advective_measured = have_face_flux
      element_base_flux%advective_base_face   = 0.0d0
      element_base_flux%advective_solved_face = 0.0d0
      if (have_face_flux) then
         call face_mass_flux_cgs(Frho_eval, msum, Frho_face)
         sv_faces_bounded  = n_species_faces_bounded
         sv_face_excursion = species_face_excursion
         call species_face_fraction(Xhe, Frho_face, Yface)
         n_species_faces_bounded = sv_faces_bounded
         species_face_excursion  = sv_face_excursion
         element_base_flux%advective_base_face = Frho_face(0)*Yface(0)
         element_base_flux%advective_solved_face =                        &
                                Frho_face(jlo-1)*Yface(jlo-1)
      endif
      element_base_flux%diffusive_base_face = Jf(0)
      element_base_flux%diffusive_solved_face = Jf(jlo-1)
      if (present(Jface_out)) Jface_out = Jf

      ! --- project the new element totals back into the species vector
      call project_elements(f_sp, Xhe, msum, .true.)

      ! --- trace metals: each element diffuses against the (post-projection)
      ! hydrogen background with its own mass and mean charge.  Carried over
      ! unchanged from the validated trace kernel (it is a trace treatment by
      ! construction: a metal element cannot feed back on the background it
      ! diffuses through); only the background definition is the element count
      ! above rather than HI+HII.  Default OFF.
      if (he_metal_diffusion .and. thereis_metals) then
         trace_ratio_under_zero = -1.0d0
         call element_nucleus_counts(f_sp, nucH, nucHe)
         nH_phys = nucH*rho*n0
         where (nH_phys .lt. 1.0d-30) nH_phys = 1.0d-30
         call carrier_counts(f_sp, carH, carHe, mcarH)
         ntot_phys = (carH + carHe)*rho*n0
         where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
         call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
         ne_phys = max(ne_rel*rho*n0, 1.0d0)
         call carrier_fractions(f_sp, n_hcar, hcar_isp, yH)
         do im = 1, n_melem
            i0m = melem_i0(im)
            top = melem_top(im)
            call trace_element_transport_coefficients(im, f_sp, rho, TK,   &
                     ntot_phys, ne_phys, yH, mcarH, zb1, eEf, dlnpsi,      &
                     nX, nXold, DcoX, GcoX)
            fXbase = nXold(1)/nH_phys(1)                  ! reservoir metal/H
            ! The element's MASS fraction per unit of the mixing ratio the
            ! row is written in, Y_X = A_X n_X/(n0 rho msum) = wYtr fX, and
            ! the factor that returns the code divergence of that mass flux
            ! to a mixing ratio per second: n0 msum v0/R0 makes it a mass
            ! rate, 1/A_X counts nuclei and 1/n_H takes the ratio.
            wYtr  = melem_A(im)*nucH/msum
            cadvX = n0*msum*v0/R0/melem_A(im)/max(nH_phys, 1.0d-30)
            call solve_trace_element_in_hydrogen(nX, nH_phys, DcoX, GcoX,  &
                                                 fXbase, dt_phys, rp,       &
                                                 Frho, wYtr, cadvX, advect)
            ! The metal's diffusive mass flux the step carried: the face
            ! coefficients the solve used and the mixing ratio it returned.
            if (present(JXface_out)) then
               call trace_face_coefficients(nH_phys, DcoX, GcoX, rp,       &
                                            PLX, PRX)
               fXnew = nX/max(nH_phys, 1.0d-30)
               do j = 0, N
                  JXface_out(j,im) = melem_A(im)*mu*                      &
                                     (PLX(j)*fXnew(j) + PRX(j)*fXnew(j+1))
               enddo
            endif
            do j = 1-Ng, N+Ng
               ! Target density from the solved mixing ratio.  There is no
               ! cap at the reservoir ratio: settling piles an element up as
               ! readily as it depletes one, and clipping the pile-up is the
               ! same limiter the memo (section 1) rejects for helium.
               rXsc = nX(j)
               if (rXsc .lt. 0.0d0) rXsc = 0.0d0
               if (nXold(j) .gt. 1.0d-25*nH_phys(j)) then
                  rXsc = rXsc/nXold(j)                    ! scale factor
                  do k = 0, top
                     f_sp(j,mion_fsp(i0m+k)) = f_sp(j,mion_fsp(i0m+k))*rXsc
                  enddo
               else if (rXsc .gt. 1.0d-25*nH_phys(j)) then
                  ! element returned to an exhausted cell: re-seed (neutral)
                  f_sp(j,mion_fsp(i0m)) = rXsc/max(rho(j)*n0, 1.0d-30)
                  do k = 1, top
                     f_sp(j,mion_fsp(i0m+k)) = 0.0d0
                  enddo
               endif
            enddo
         enddo
         ! THE HYDROGEN BACKGROUND RECOILS AGAINST THE METAL FLUXES.  In a
         ! single-fluid mixture the diffusive mass fluxes sum to zero: the
         ! helium flux is balanced by component 1 in the binary solve above,
         ! and a trace element that settles relative to the hydrogen it
         ! diffuses through has to be balanced by that same hydrogen.  The
         ! metal loop wrote densities and no counter-flux, so the mixture
         ! ended carrying a mass its own rho does not have (3.0e-3 of it on
         ! the relaxed synthetic column, MEASURED).  Closing it here puts the
         ! metals' change of mass back on the hydrogen group, at the helium
         ! mass fraction the solve returned, and leaves the metal densities
         ! exactly where their own transport put them.
         call project_elements(f_sp, Xhe, msum, .false.)
      endif

      ! --- elemental census of the metals.  These equations have no sink for a
      ! metal nucleus: an element that is present in the reservoir and absent
      ! from a cell that holds hydrogen cannot have got there by physics.  It is
      ! checked on every step because nothing else in the run can see it --
      ! the steady residual is a residual of the hydro and energy equations,
      ! whose solution with the metals removed is a perfectly good solution of
      ! the equations as posed, and the elemental-flux closure measures a window
      ! that need not contain the cells concerned (section 84).
      call metal_hydrogen_ratio_departure(f_sp, qdep, n_vanished)

      ! --- IS WHAT THE STEP PRODUCED AN ADMISSIBLE COMPOSITION?  The density
      ! is the hydrodynamics' and the step held it fixed, so the species it
      ! hands back have to weigh it: sum_i m_i n_i = rho, which in the mass
      ! fractions of this module reads msum = 1 (mixture_mass_split carries
      ! the mass policy of calc_rho term for term).  Measured on the state
      ! that would be handed back, together with the two tests that finite
      ! and bounded populations make, and the update stands or the entry
      ! composition does.
      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum_out)
      closure = maxval(abs(msum_out(1:N) - 1.0d0))
      element_mass_closure_departure = closure
      outcome = element_step_accepted
      if (.not. every_species_is_finite(f_sp)) then
         outcome = element_step_nonfinite
      else if (max(Xover, Xunder) .gt. element_fraction_bound_tol) then
         outcome = element_step_out_of_bounds
      else if (closure - maxval(abs(msum(1:N) - 1.0d0))                   &
               .gt. element_mass_closure_tol) then
         outcome = element_step_mass_closure_failed
      endif
      if (outcome .ne. element_step_accepted) f_sp = f_entry
      element_step_last_status = outcome
      if (present(status)) status = outcome
      deallocate(f_entry)

      if (diffusion_check_on()) then
         call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum_out)
         call report_step(Xhe, rho_phys,                                  &
              maxval(abs(msum_out(1:N) - msum(1:N))/msum(1:N)),           &
              Xover, Xunder, qdep, n_vanished)
      endif

      end subroutine element_diffusion_step

      ! ------------------------------------------------------------------ !

      subroutine element_base_flux_report(line, evaluation_expected)
      ! The base budget entry written as one line, with the evaluation it
      ! belongs to named in it.  A report is a statement about a state, so
      ! a record no operator has written is said to be absent instead of
      ! being printed as a flux of zero, and a record written by another
      ! evaluation than the one the caller is reporting on is said to be
      ! that one: evaluation_expected (optional) is the count the caller
      ! read before its own operator call, and the line says whether the
      ! record moved past it.
      !
      ! Both halves are helium element mass fluxes in g cm^-2 s^-1, positive
      ! outward.  The advective half is F_rho Y_face on the Riemann face mass
      ! flux the evaluation's transport rode on, and an evaluation that had
      ! none is said to have none; the diffusive half at the base face is
      ! zero by the discretization (the declaration of the record says why),
      ! and the one at the first solved face is the flux the Dirichlet
      ! reservoir composition drives.
      character(len=*), intent(out) :: line
      integer, optional, intent(in) :: evaluation_expected

      if (element_base_flux%evaluation .eq. 0) then
         line = 'base element flux: no element operator evaluation has '// &
                'written a record'
         return
      endif
      if (element_base_flux%advective_measured) then
         write(line,'(A,I0,A,ES15.7,A,ES15.7,A,I0,A,ES15.7,A,ES15.7)')    &
            'base element flux [g/cm2/s, positive outward] of '//         &
            'evaluation ', element_base_flux%evaluation,                  &
            ': base face 0 advective ',                                   &
            element_base_flux%advective_base_face,                        &
            ' diffusive ', element_base_flux%diffusive_base_face,         &
            ' ; first solved face ', element_base_flux%first_solved_face, &
            ' advective ', element_base_flux%advective_solved_face,       &
            ' diffusive ', element_base_flux%diffusive_solved_face
      else
         write(line,'(A,I0,A,ES15.7,A,I0,A,ES15.7,A)')                     &
            'base element flux [g/cm2/s, positive outward] of '//         &
            'evaluation ', element_base_flux%evaluation,                  &
            ': base face 0 diffusive ',                                   &
            element_base_flux%diffusive_base_face,                        &
            ' ; first solved face ', element_base_flux%first_solved_face, &
            ' diffusive ', element_base_flux%diffusive_solved_face,       &
            ' ; advective not measured: the evaluation was handed no '//  &
            'face mass flux and followed no Runge-Kutta stages'
      endif
      if (present(evaluation_expected)) then
         if (evaluation_expected .ne. element_base_flux%evaluation)        &
            line = trim(line)//' (NOT the evaluation asked for)'
      endif
      end subroutine element_base_flux_report

      ! ------------------------------------------------------------------ !

      logical function every_species_is_finite(f_sp) result(ok)
      ! Whether a species vector is a vector of ordinary reals: no NaN and no
      ! infinity anywhere in it.  NaN fails every comparison including with
      ! itself and an infinity exceeds the largest representable finite
      ! value, so the two tests together cover both; written out rather than
      ! taken from ieee_arithmetic so that the generated module dependency
      ! graph stays over the source tree, as finite_real
      ! (ionization_equilibrium.f90) is.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      ok = all(f_sp .eq. f_sp) .and. all(abs(f_sp) .le. huge(1.0d0))
      end function every_species_is_finite

      ! ------------------------------------------------------------------ !

      pure function element_step_outcome_text(outcome) result(text)
      ! One of the five names element_step_last_status (or a status
      ! argument) can carry, in words, for a caller that reports rather
      ! than branches on the integer (the marching controller's log line).
      integer, intent(in) :: outcome
      character(len=32)   :: text
      select case (outcome)
      case (element_step_accepted)
         text = 'accepted'
      case (element_step_solve_failed)
         text = 'nonlinear solve did not converge'
      case (element_step_nonfinite)
         text = 'candidate composition not finite'
      case (element_step_out_of_bounds)
         text = 'helium fraction left [0,1]'
      case (element_step_mass_closure_failed)
         text = 'mass closure worsened'
      case default
         text = 'unrecognized outcome'
      end select
      end function element_step_outcome_text

      ! ------------------------------------------------------------------ !

      subroutine metal_hydrogen_ratio_departure(f_sp, qdep, n_vanished)
      ! Departure of the metal/hydrogen nucleus ratio from the reservoir
      ! abundance melem_ab, over the cells 1..N and the elements the reservoir
      ! states.  qdep is max |(n_X/n_H)/melem_ab - 1| and n_vanished counts the
      ! (cell, element) pairs in which the element has vanished outright while
      ! hydrogen is present.
      !
      ! A nonzero n_vanished is not a tolerance being exceeded, it is an element
      ! that is gone: the transport equation moves metal nuclei with the
      ! hydrogen they are slaved to and destroys none, so the count is zero on
      ! any state these equations can reach.  It is therefore reported as a
      ! warning the first time it happens, whether or not the step diagnostic is
      ! on.  qdep is a magnitude and not an error: he_metal_diffusion moves the
      ! metals relative to hydrogen on purpose, and a restart may carry a
      ! settled column, so only its report is unconditional -- under
      ! EXHALE_DIFFUSION_CHECK -- and not any judgement of it.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8,                                 intent(out) :: qdep
      integer,                                intent(out) :: n_vanished

      real*8, dimension(1-Ng:N+Ng) :: nucH_l, nucHe_l
      real*8  :: nX
      integer :: j, ie, k, jfirst, iefirst
      logical, save :: warned = .false.

      qdep       = 0.0d0
      n_vanished = 0
      if (.not. thereis_metals) return
      if (.not. allocated(melem_ab)) return

      call element_nucleus_counts(f_sp, nucH_l, nucHe_l)
      jfirst  = 0
      iefirst = 0
      do ie = 1, n_melem
         if (melem_ab(ie) .le. 0.0d0) cycle
         do j = 1, N
            if (nucH_l(j) .le. 1.0d-30) cycle
            nX = 0.0d0
            do k = 0, melem_top(ie)
               nX = nX + f_sp(j,mion_fsp(melem_i0(ie)+k))
            enddo
            if (nX .le. 0.0d0) then
               n_vanished = n_vanished + 1
               if (jfirst .eq. 0) then
                  jfirst  = j
                  iefirst = ie
               endif
            endif
            qdep = max(qdep, abs(nX/(nucH_l(j)*melem_ab(ie)) - 1.0d0))
         enddo
      enddo

      if (n_vanished .gt. 0 .and. .not. warned) then
         warned = .true.
         write(*,'(A,I0,A)') ' (element diffusion) WARNING: ',            &
              n_vanished, ' (cell, element) pairs hold no metal nuclei'// &
              ' at all while hydrogen is present.'
         write(*,'(A,I0,A,ES12.5,A,A,A)') '   first at cell ', jfirst,    &
              ', r = ', r(jfirst), ' R_p, element ',                      &
              trim(melem_name(iefirst)),                                  &
              '.  These equations have no sink for a metal nucleus.'
      endif

      end subroutine metal_hydrogen_ratio_departure

      ! ------------------------------------------------------------------ !

      subroutine relax_element_composition(rho, Tcode, f_sp, Frho,         &
                                           omega, map_distance, nstep,     &
                                           status, displacement)
      ! Relax the element composition to its steady state in a FIXED wind.
      !
      ! The step size is the COMPOSITION time scale, not the hydro CFL step.
      ! The implicit operator of element_diffusion_step is unconditionally
      ! stable, so nothing ties its dt to a sound-crossing time, while the
      ! quantities it relaxes evolve on
      !
      !    dr^2/(D_12 + K_zz)   (diffusion)   and   dr/max(|v|, w_s)
      !
      ! with w_s = D_12 |G| the settling drift speed.  At the base of
      ! HD 209458 b the first of these is ~10^6 s against a ~1 s CFL step, so
      ! a relaxation driven at the CFL step moves the composition by ~10^-6 of
      ! the way per step and no practical number of steps reaches the steady
      ! state -- which is what left the outer loop drifting at ~0.75 per pass
      ! and hitting its pass limit.  Here dt starts at the smallest of those
      ! scales and is grown geometrically, so the late steps are direct steady
      ! solves and the coefficients are Picard-updated along the way.
      !
      ! THE ADVECTING FLOW IS THE FACE MASS FLUX OF THIS STATE, GIVEN BY THE
      ! CALLER (Frho), and the advective term is the divergence of the face
      ! element mass fluxes it carries -- the same object, through the same
      ! three routines, that the Runge-Kutta stages advect the elements with
      ! and that the stationary elemental row of element_transport_residual
      ! balances.  The relaxation's fixed point is then the zero of that row:
      ! the elemental Picard alternation and the row are one operator, which
      ! is what makes a converged alternation a state the row reads as
      ! stationary.  The caller reads the flux out of the mass row of the
      ! state it hands over (face_mass_flux_of_state), which refuses a state
      ! whose mass row was never assembled.
      !
      ! It is NOT the non-conservative form rho v dX/dr on a smoothed steady
      ! mass flux mdot/(4 pi r^2).  That form is the conservative equation
      ! only where the CELL-CENTRED rho and v satisfy continuity, and they do
      ! not near the base: below ~1.02 R_p the spread of the cell-centred
      ! r^2 rho v is 10^2-10^4 times its own median (measured on HD 209458 b
      ! and LHS 1140 b), because the cell-centred velocity carries a
      ! collocated two-cell odd-even mode whose sign alternates from cell to
      ! cell and which the Riemann face mass flux does not carry, and the
      ! cell-velocity form relaxed on that field converged to the composition
      ! of a flow that neither conserves mass nor exists.  The conservative
      ! form needs no such repair: div(F_rho X) = X div(F_rho) + F_rho grad X,
      ! so an element
      ! inherits the mass row's own imbalance instead of having it subtracted
      ! out, and a cell can only lose the fraction of its mass the mass row
      ! loses.  The marching path carries the same divergence inside the
      ! stages and passes no flux here.
      !
      ! THE STOPPING MEASURE RUNS OVER EVERY ELEMENT THE PASS RELAXED, not
      ! over helium alone: the largest change of an element's mass fraction
      ! over a step, divided by that element's own reservoir value, maximized
      ! over the elements.  Helium and the trace metals settle against
      ! different masses, charges and diffusion coefficients, so a column
      ! whose metals move on a longer time scale than its helium would be
      ! declared converged while they were still travelling, and the outer
      ! Picard loop would damp and test a distance that is not the distance
      ! to the fixed point.  Each element is measured in its own reservoir
      ! value because the mass fractions differ by four orders of magnitude;
      ! an element the run holds none of at the base carries no equation and
      ! is left out of the maximum.
      !
      ! The measure is ABSOLUTE in that reservoir value.  The relative
      ! measure it replaces, max|dX|/X, is meaningless in a cell the
      ! transport has emptied.
      !
      ! UNDER-RELAXATION.  omega damps the composition update handed back to
      ! the outer Picard loop,
      !
      !    X <- X_old + omega (X_relaxed - X_old),      0 < omega <= 1
      !
      ! because that loop -- solve the wind at fixed composition, relax the
      ! composition at fixed wind -- has no damping of its own and the two can
      ! chase each other instead of converging (measured on the HD 209458 b
      ! Kzz = 0 wind: a limit cycle at 0.15-0.17 over twenty passes).  The
      ! caller owns the schedule.  The returned map_distance is the UNDAMPED
      ! distance max|X_relaxed - X_old|/X_base, not the damped step actually
      ! applied, so omega cannot make it small by moving less; the applied
      ! step is returned separately as displacement.  Neither of the two is
      ! the residual of the returned state: the endpoint of this map is the
      ! fixed point only when the relaxation reached it, and
      ! element_transport_residual_norm measures the elemental transport
      ! balance of the composition that leaves here.
      !
      ! ONLY ADMISSIBLE STEPS ARE KEPT, AND A SMALL MOVEMENT IS NOT A
      ! CONVERGENCE PROOF.  Each transport step is taken on a copy and asked
      ! for its own outcome (element_diffusion_step): a step that did not
      ! solve its row, or whose composition no longer carries the density the
      ! step held fixed, is discarded, the step length is halved and the step
      ! is retried, within the two budgets he_relax_retry_max and
      ! he_relax_grow_floor.  With neither budget left the composition this
      ! pass was given is restored.  The three outcomes the caller is told
      ! apart are the fixed point of the operator
      ! (element_relaxation_converged, the movement measure met), a partial
      ! advance that ran out of steps (element_relaxation_step_budget), and
      ! no admissible advance at all (element_relaxation_failed).  A caller
      ! that asks for none of them still gets the restored composition, and
      ! the refusal is announced once.
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      ! The face mass flux F_rho(j) at r_edg(j) in code units, the one the
      ! Riemann solve of this state returned and the mass row of this state
      ! differences.  Fixed for the whole relaxation, exactly as the wind it
      ! belongs to is.
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: Frho
      real*8,                                 intent(in)    :: omega
      ! THE DISTANCE TO THE ENDPOINT OF THIS INNER MAP, undamped:
      ! max|X_relaxed - X_entry| in each element's own reservoir value,
      ! formed before omega damps the update that is applied.  It is the
      ! distance to the endpoint of a FINITE relaxation, which may exit on
      ! its step budget, so it is not the residual of the state handed back;
      ! element_transport_residual_norm measures that residual.
      real*8,                                 intent(out)   :: map_distance
      integer,                                intent(out)   :: nstep
      integer, optional,                      intent(out)   :: status
      ! WHAT THIS PASS ACTUALLY CHANGED, in the same norm: the distance from
      ! the composition the pass was entered with to the one it hands back,
      ! after the damping and the projection.  It is zero when the entry
      ! composition was restored, and equals omega times the map distance
      ! where the damping is the only thing between the two (helium; a trace
      ! element is not damped).
      real*8, optional,                       intent(out)   :: displacement

      real*8, dimension(:,:), allocatable :: f_pass, f_try
      real*8, dimension(1-Ng:N+Ng) :: dt_code, nucH, nucHe, Xmix
      real*8, dimension(1-Ng:N+Ng) :: msum, mass1, msum_ret
      real*8, dimension(1-Ng:N+Ng) :: TK, Dcl, Gcl, dmcl, zbcl
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, mflx, vadv
      ! The mass fractions of every element this relaxation moves -- helium
      ! in column 1, each trace metal in column 1+im -- at the start of the
      ! pass, at the previous step and now, and the reservoir value each is
      ! measured in.
      real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Ynow, Yprev, Ypass
      real*8, dimension(1+n_melem) :: Yres
      real*8 :: m_1, X_base, tscale, drj, tdiff, tadv, wset, grow, dt_ref
      real*8 :: mflx_median
      integer :: im, j, k, st, outcome, nreject, nkept
      logical, save :: warned_refusal = .false.
      logical, save :: warned_retry_budget = .false.

      map_distance = 0.0d0
      nstep        = 0
      if (present(displacement)) displacement = 0.0d0
      if (present(status)) status = element_relaxation_converged
      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      m_1    = mass_per_H_nucleus_without_He()
      X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)
      tscale = R0/v0

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      ! Only the composition TIME SCALES below are built from this; the
      ! operator recomputes its own D_12 each step from the current state.
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dcl)
      call settling_coefficient(rho, Tcode, f_sp, Gcl, dmcl, zbcl)

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      rho_phys = rho*n0*mu*msum

      ! The material speed the faces of a cell carry, for the advective time
      ! scale below.  It is read from the same flux the operator advects on
      ! and not from the cell velocity, so the step is measured against the
      ! transport that is actually applied.
      do j = 1, N
         vadv(j) = 0.5d0*(abs(Frho(j)) + abs(Frho(j-1)))                   &
                   *n0*mu*msum(j)*v0/max(rho_phys(j), 1.0d-30)
      enddo

      do j = 1, N
         drj   = max((r_edg(j) - r_edg(j-1))*R0, 1.0d0)
         tdiff = drj*drj/max(Dcl(j) + kzz_cell(j), 1.0d-30)
         wset  = Dcl(j)*abs(Gcl(j))
         tadv  = drj/max(vadv(j), wset, 1.0d-30)
         dt_code(j) = min(tdiff, tadv)/tscale
      enddo
      dt_code(1-Ng:0)   = dt_code(1)
      dt_code(N+1:N+Ng) = dt_code(N)
      dt_ref = minval(dt_code(1:N))

      call element_mass_fractions(f_sp, Ypass)
      ! The reservoir each element is measured in.  Helium's is the input
      ! He/H, which is what the operator pins its base to; a trace element's
      ! is its own base cell, the Dirichlet reservoir the trace solve imposes
      ! and the one cell the relaxation never moves.  A zero there is an
      ! element the run holds none of: no equation, and no term in the
      ! maximum.
      Yres    = 0.0d0
      Yres(1) = X_base
      if (he_metal_diffusion .and. thereis_metals) then
         do im = 1, n_melem
            Yres(1+im) = Ypass(1,1+im)
         enddo
      endif
      Yprev   = Ypass
      Ynow    = Ypass
      grow    = 1.0d0
      nreject = 0
      ! A budget exit is the outcome until a step meets the movement measure:
      ! the loop below can only leave by that measure, by a refusal, or by
      ! running out of steps.
      outcome = element_relaxation_step_budget
      allocate(f_pass(1-Ng:N+Ng,n_species), f_try(1-Ng:N+Ng,n_species))
      f_pass = f_sp
      k      = 0
      do while (k .lt. he_relax_maxstep)
         f_try = f_sp
         call element_diffusion_step(rho, Tcode, f_try, dt_code*grow,      &
                                     Frho_in = Frho, status = st)
         if (st .ne. element_step_accepted) then
            ! The operator could not advance the composition over this step.
            ! A shorter step is a better-conditioned one: the time term grows
            ! against the flux divergence it cancels, which is what sets the
            ! conditioning of the tridiagonal solve (solve_mass_fraction).
            if (nreject .ge. he_relax_retry_max .or.                      &
                grow .le. he_relax_grow_floor) then
               ! A RELAXATION THAT HAS KEPT STEPS ENDS AS A PARTIAL ADVANCE.
               ! Each kept step passed element_diffusion_step's own test,
               ! and the outcome table above reserves the refusal for a
               ! relaxation that kept none; ending here as a refusal threw
               ! the kept steps away with the entry composition restored,
               ! and made every later pass a copy of this one.  The retry
               ! budget counts every rejection of the relaxation.
               if (k .gt. 0) then
                  outcome = element_relaxation_step_budget
                  if (.not. warned_retry_budget) then
                     warned_retry_budget = .true.
                     write(*,'(A,I0,A,I0,A)') ' (element diffusion) the'// &
                          ' retry budget ended a relaxation after ', k,    &
                          ' kept step(s) (last step outcome ', st,         &
                          '); the kept steps are handed back as a'//      &
                          ' partial advance'
                  endif
               else
                  outcome = element_relaxation_failed
               endif
               exit
            endif
            nreject = nreject + 1
            grow    = 0.5d0*grow
            cycle
         endif
         f_sp  = f_try
         k     = k + 1
         nstep = k
         call element_mass_fractions(f_sp, Ynow)
         if (element_composition_distance(Ynow, Yprev, Yres)              &
             .lt. he_relax_tol) then
            outcome = element_relaxation_converged
            exit
         endif
         Yprev = Ynow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo
      if (outcome .eq. element_relaxation_failed) then
         ! Nothing this pass did is part of the state handed back, so no step
         ! is reported as taken and the distance to the fixed point is not
         ! measured from a state the caller never receives.
         f_sp         = f_pass
         map_distance = 0.0d0
         nkept        = nstep
         nstep = 0
         deallocate(f_pass, f_try)
         if (present(status)) status = outcome
         if (.not. warned_refusal) then
            warned_refusal = .true.
            write(*,'(A,I0,A,I0,A,ES10.3,A,ES10.3,A)')                    &
                 ' (element diffusion) WARNING: the composition'//        &
                 ' relaxation took ', nkept, ' step(s) and then found'//  &
                 ' no admissible one (outcome ', st,                      &
                 ' at a step length ', grow,                              &
                 ' of the composition time scale, mass-closure'//         &
                 ' departure ', element_mass_closure_departure,           &
                 '); the composition this pass was given is what is'//    &
                 ' handed back.'
         endif
         return
      endif
      deallocate(f_try)
      map_distance = element_composition_distance(Ynow, Ypass, Yres)

      ! Damped update: blend the relaxed composition with the one this pass
      ! started from and project the blend back into the species vector.  The
      ! projection is the same one element_diffusion_step ends on, applied to
      ! the CURRENT f_sp (which carries X_relaxed), so both element totals are
      ! met exactly and no species can go negative.
      !
      ! ONLY HELIUM IS DAMPED.  X is the member of the normalized pair, and
      ! the projection below is the map from it back into the species vector;
      ! a trace element is slaved to the hydrogen background and carries no
      ! such projection, so the trace elements keep the composition the
      ! relaxation left.  The drift reported to the outer loop covers every
      ! element either way, because it measures the distance to the fixed
      ! point and not the step applied.
      if (omega .lt. 1.0d0) then
         Xmix = Ypass(:,1) + omega*(Ynow(:,1) - Ypass(:,1))
         where (Xmix .lt. 0.0d0) Xmix = 0.0d0
         where (Xmix .gt. 1.0d0) Xmix = 1.0d0
         call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
         call project_elements(f_sp, Xmix, msum, .true.)
      endif

      ! The state that leaves this routine stands no farther from the density
      ! it was given than the composition this pass began with, the damped
      ! blend included: the projection above is the inverse of the split it
      ! reads its mass from, so the departure is an invariant of the pass and
      ! the measure is that inversion's round-off.
      call mixture_mass_sum(f_sp, msum_ret)
      element_mass_closure_departure = maxval(abs(msum_ret(1:N) - 1.0d0))
      call mixture_mass_sum(f_pass, msum_ret)
      if (element_mass_closure_departure - maxval(abs(msum_ret(1:N)        &
          - 1.0d0)) .gt. element_mass_closure_tol                          &
          .or. .not. every_species_is_finite(f_sp)) then
         f_sp         = f_pass
         outcome      = element_relaxation_failed
         map_distance = 0.0d0
         nstep        = 0
         if (.not. warned_refusal) then
            warned_refusal = .true.
            write(*,'(A,ES10.3,A)') ' (element diffusion) WARNING: the'//  &
                 ' composition the relaxation assembled does not carry'// &
                 ' the density it was given (departure ',                 &
                 element_mass_closure_departure,                          &
                 '); the composition this pass was given is what is'//    &
                 ' handed back.'
         endif
      endif
      ! WHAT LEAVES THE ROUTINE, MEASURED ON WHAT LEAVES THE ROUTINE.  The
      ! composition here is the one the caller receives: the damped blend
      ! where the damping applied, the relaxed one where it did not, and the
      ! entry composition where the closure above put it back.
      if (present(displacement)) then
         call element_mass_fractions(f_sp, Ynow)
         displacement = element_composition_distance(Ynow, Ypass, Yres)
      endif
      deallocate(f_pass)
      if (present(status)) status = outcome

      ! dt_ref is the shortest composition time scale of the grid; report it
      ! so the log shows what the relaxation was measured against, together
      ! with the mass flux it advected on -- the median of 4 pi r^2 F_rho
      ! over the escape window, which is the wind's own mass-loss rate as the
      ! faces of the solver carry it.
      if (diffusion_check_on()) then
         do j = 1-Ng, N+Ng
            mflx(j) = 4.0d0*pi*(r_edg(j)*R0)**2*Frho(j)*n0*mu*msum(j)*v0
         enddo
         mflx_median = median_of(mflx(j_min:N))
         write(0,'(A,ES10.3,A,ES12.5,A,I0)')                              &
              ' (diffusion) relaxation dt_0 = ', dt_ref*tscale,           &
              ' s, wind 4pi r^2 F_rho = ', mflx_median,                   &
              ' g/s, steps = ', nstep
      endif

      end subroutine relax_element_composition

      ! ------------------------------------------------------------------ !

      real*8 function element_composition_distance(Ya, Yb, Yres)          &
             result(d)
      ! How far apart two element compositions are: the largest change of an
      ! element's mass fraction over the physical cells, in that element's
      ! own reservoir value, maximized over the elements the state carries.
      !
      ! Each element is divided by its OWN reservoir because the transported
      ! mass fractions span four orders of magnitude between helium and a
      ! trace metal, so an absolute maximum over the columns would read the
      ! helium column alone and the metals would never reach the measure.  An
      ! element whose reservoir is zero is one the run holds none of: it has
      ! no equation and is left out.
      real*8, dimension(1-Ng:N+Ng,1+n_melem), intent(in) :: Ya, Yb
      real*8, dimension(1+n_melem),           intent(in) :: Yres
      integer :: ie
      d = 0.0d0
      do ie = 1, 1 + n_melem
         if (Yres(ie) .le. 0.0d0) cycle
         d = max(d, maxval(abs(Ya(1:N,ie) - Yb(1:N,ie)))/Yres(ie))
      enddo
      end function element_composition_distance

      ! ------------------------------------------------------------------ !

      logical function residual_norm_is_admissible(x) result(ok)
      ! A RESIDUAL NORM IS A FINITE NONNEGATIVE NUMBER, and a progress
      ! control may read it only when it is one.
      !
      ! x .eq. x rejects a NaN and nothing else: an infinity compares equal
      ! to itself and would be read as a measurement of a balance, and a
      ! negative value is not a norm of anything.  Both are refusals of the
      ! MEASUREMENT and say nothing about the state, so a caller that meets
      ! one reports the measure unavailable rather than substituting another
      ! quantity for it.  The two tests are written out rather than taken
      ! from ieee_arithmetic, for the reason every_species_is_finite gives.
      real*8, intent(in) :: x
      ok = (x .eq. x) .and. (abs(x) .le. huge(1.0d0)) .and. (x .ge. 0.0d0)
      end function residual_norm_is_admissible

      ! ------------------------------------------------------------------ !

      subroutine element_transport_residual_norm(rho, Tcode, f_sp, Frho,  &
                                 rnorm, jworst, ielem, res_abs, scale_abs,&
                                 measured)
      ! THE ELEMENTAL TRANSPORT BALANCE OF A STATE, AS ONE NUMBER:
      !
      !    rnorm = max_j |res_j| / max(scale_j, floor)
      !
      ! over the helium partition row and over every trace element row the
      ! state carries an equation for, with res and scale exactly those of
      ! element_transport_residual (transport against transport at an
      ! infinite step length, the scale the sum of the row's own terms).
      ! The reduction is the one certification_row_measure forms for the
      ! same rows, so the number read here and the elemental entries of the
      ! certification of the same state are one measure and not two
      ! spellings of it.
      !
      ! WHY A RESIDUAL AND NOT THE INNER MAP'S DISTANCE.  The distance
      ! relax_element_composition returns is the distance to the endpoint of
      ! a finite relaxation, which may exit on its step budget, so it can be
      ! small at a state whose elemental rows are not satisfied.  The
      ! residual is a property of the state handed back and of nothing else.
      !
      ! The measurement is stationary and side-effect free: the operator
      ! below holds its own last-solve diagnostics aside and puts them back,
      ! so no isolated workspace is needed around this call.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      ! The face mass flux of the state the balance rides on, the one its
      ! own mass row returned (face_mass_flux_of_state).
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Frho
      real*8,                                 intent(out) :: rnorm
      ! Where the maximum stands, and in which element's row: 0 is the
      ! helium partition, a positive value the trace element's index.
      integer,                                intent(out) :: jworst, ielem
      ! The two absolute terms behind the ratio at that cell, so a movement
      ! of the scale is visible beside the ratio it divides.
      real*8,                                 intent(out) :: res_abs
      real*8,                                 intent(out) :: scale_abs
      ! False where the state carries no elemental equation at all, or where
      ! the operator could not measure one: the caller then has no element
      ! measure for this state and says so, rather than reading another
      ! quantity in its place.
      logical,                                intent(out) :: measured

      real*8, dimension(1:N)         :: rhe, she
      real*8, dimension(1:N,n_melem) :: rtr, str
      logical, dimension(n_melem)    :: carried
      logical :: ok_he, ok_tr
      real*8  :: q
      integer :: j, im

      rnorm     = 0.0d0
      jworst    = 0
      ielem     = 0
      res_abs   = 0.0d0
      scale_abs = 0.0d0
      measured  = .false.
      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      call element_transport_residual(rho, Tcode, f_sp, Frho, rhe, she,   &
               ok_he, rtr, str, ok_tr, tr_carried = carried)
      if (.not. ok_he) return
      do j = 1, N
         q = abs(rhe(j))/max(she(j), element_row_scale_floor)
         if (q .gt. rnorm) then
            rnorm     = q
            jworst    = j
            ielem     = 0
            res_abs   = abs(rhe(j))
            scale_abs = she(j)
         endif
      enddo
      measured = .true.
      if (.not. ok_tr) return
      do im = 1, n_melem
         if (.not. carried(im)) cycle
         do j = 1, N
            q = abs(rtr(j,im))/max(str(j,im), element_row_scale_floor)
            if (q .gt. rnorm) then
               rnorm     = q
               jworst    = j
               ielem     = im
               res_abs   = abs(rtr(j,im))
               scale_abs = str(j,im)
            endif
         enddo
      enddo
      end subroutine element_transport_residual_norm

      ! ------------------------------------------------------------------ !

      real*8 function median_of(a)
      ! Median of a, by insertion sort of a local copy (the arrays here are
      ! one grid long and this runs once per relaxation).
      real*8, dimension(:), intent(in) :: a
      real*8, dimension(size(a)) :: b
      real*8  :: t
      integer :: i, k, m
      b = a
      do i = 2, size(b)
         t = b(i)
         k = i - 1
         do while (k .ge. 1)
            if (b(k) .le. t) exit
            b(k+1) = b(k)
            k = k - 1
         enddo
         b(k+1) = t
      enddo
      m = size(b)
      if (mod(m,2) .eq. 0) then
         median_of = 0.5d0*(b(m/2) + b(m/2+1))
      else
         median_of = b((m+1)/2)
      endif
      end function median_of

      ! ------------------------------------------------------------------ !

      subroutine element_nucleus_counts(f_sp, nucH, nucHe)
      ! Element (nucleus) counts per unit mass, over every species that
      ! carries them: n_H/rho = f_HI + f_HII + 2 f_H2 + 2 f_H2+ + 3 f_H3+
      ! + f_HeH+, n_He/rho = f_HeI + f_HeII + f_HeIII + f_HeH+.  Columns a
      ! run does not carry are zero, so no flag has to be tested.  He 2^3S
      ! is skipped: it is an excited level of He I and its nucleus is
      ! already inside the HeI column (bsp_is_excited_level).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nucH, nucHe
      integer :: ib

      nucH  = 0.0d0
      nucHe = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nH(ib)  .gt. 0)                                          &
            nucH  = nucH  + dble(bsp_nH(ib)) *f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0)                                          &
            nucHe = nucHe + dble(bsp_nHe(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      where (nucH  .lt. 0.0d0) nucH  = 0.0d0
      where (nucHe .lt. 0.0d0) nucHe = 0.0d0

      end subroutine element_nucleus_counts

      ! ------------------------------------------------------------------ !

      subroutine carrier_counts(f_sp, carH, carHe, mcarH)
      ! Collision-partner (CARRIER) counts per unit mass of the two
      ! components, and the mean hydrogen mass per carrier of component 1:
      !
      !   carH  = sum_c f_c            over hcar_isp  = HI HII H2 H2+ H3+
      !   carHe = sum_s f_s            over hecar_isp = HeI HeII HeIII
      !   mcarH = sum_c m_c f_c / carH   [m_H per carrier]
      !
      ! The friction and the force terms of the diffusion equation count
      ! PARTICLES, not nuclei: an H2 molecule is one collision partner and
      ! carries twice the weight of an H atom (memo section 5).  With the
      ! molecular chemistry off the molecular columns are zero, so carH and
      ! carHe are the nucleus counts and mcarH = 1 exactly, which is why no
      ! flag is tested here and why the atomic region is untouched.
      !
      ! HeH+ carries a nucleus of each element and is in neither list, the
      ! same trace approximation that leaves it out of the mean charges.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: carH, carHe
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: mcarH

      real*8, dimension(1-Ng:N+Ng) :: mHsum
      integer :: ic

      carH  = 0.0d0
      carHe = 0.0d0
      mHsum = 0.0d0
      do ic = 1, n_hcar
         carH  = carH  + max(f_sp(:,hcar_isp(ic)), 0.0d0)
         mHsum = mHsum + hcar_m(ic)*max(f_sp(:,hcar_isp(ic)), 0.0d0)
      enddo
      do ic = 1, n_hecar
         carHe = carHe + max(f_sp(:,hecar_isp(ic)), 0.0d0)
      enddo
      if (present(mcarH)) then
         where (carH .gt. 0.0d0)
            mcarH = mHsum/carH
         elsewhere
            mcarH = m_H_amu
         end where
      endif

      end subroutine carrier_counts

      ! ------------------------------------------------------------------ !
      !  Pair diffusion coefficients (memo 2.6, eqs. 8-10)
      !
      !  All three are the Chapman-Enskog first approximation
      !      D_st = 3 k T / (16 n mu_st Omega_st^(1,1)),
      !      Omega^(1,1) = sqrt(kT/2 pi mu) int_0^inf e^-gam^2 gam^5 Q^(1) dgam,
      !      gam^2 = mu g^2/(2 k T),
      !  evaluated with the momentum-transfer cross section Q^(1) of the
      !  interaction the pair actually has.  n is the total nucleus density in
      !  every case, so the three are directly comparable and the mixture rule
      !  below can add their inverses.
      ! ------------------------------------------------------------------ !

      real*8 function hard_sphere_pair_diffusion(TK, ntot, A_s, A_t)
      ! NEUTRAL-NEUTRAL.  Banks & Kockarts (1973), Aeronomy Part B: the
      ! binary coefficient of two neutrals of atomic weights A_s, A_t is
      !
      !    D = 1.52e18 (1/A_s + 1/A_t)^(1/2) T^(1/2) / n     [cm^2/s]
      !
      ! symmetric in the two species, so one number serves He-in-H and
      ! H-in-He.  Consistency check of the framework: the rigid-sphere
      ! collision integral Omega^(1,1) = pi d^2 (kT/2 pi mu)^(1/2) put into
      ! the Chapman-Enskog D above returns exactly this expression with the
      ! collision diameter d = 2.7 Angstrom, so (8), (9) and (10) are the same
      ! approximation and not three unrelated fits.
      real*8, intent(in) :: TK, ntot, A_s, A_t
      hard_sphere_pair_diffusion = bk_hs_pref*sqrt(1.0d0/A_s + 1.0d0/A_t)  &
                                   *sqrt(TK)/ntot
      end function hard_sphere_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function polarization_pair_diffusion(TK, ntot, alpha_n, mu_g)
      ! ION-NEUTRAL, NON-RESONANT.  The ion polarizes the neutral, giving the
      ! Langevin interaction V = -alpha_n e^2/(2 r^4).  For that potential the
      ! product g Q^(1)(g) is independent of the relative speed,
      !
      !    g Q^(1) = 2.21 pi e (alpha_n/mu)^(1/2),
      !
      ! the constant behind the non-resonant ion-neutral momentum-transfer
      ! collision frequency of Schunk & Nagy (Ionospheres, eq. 4.88),
      ! nu_in = 2.21 pi (n_n m_n/(m_i+m_n)) (gamma_n e^2/mu_in)^(1/2).  The
      ! textbook is not in references/, so the constant below was re-derived
      ! from the Chapman-Enskog integral rather than copied: with
      ! Q^(1) = C/g the integral gives Omega^(1,1) = 3 C/16, and the 3/16 of
      ! the Chapman-Enskog D cancels it exactly, leaving
      !
      !    D = k T / (2.21 pi e n (alpha_n mu)^(1/2))        [cm^2/s]
      !
      ! alpha_n [cm^3] is the neutral partner's static dipole polarizability,
      ! e is in esu and mu_g the reduced mass in grams.  Note D ~ T, not
      ! T^(1/2): the polarization coupling weakens as the pair gets faster.
      ! Valid for NON-RESONANT pairs only; a resonant pair (H+ + H, He+ + He)
      ! is charge-exchange dominated and never appears here, because both of
      ! its partners carry the same element (see the module header).
      !
      ! TEMPERATURE RANGE.  This is the LOW-energy limit of the ion-neutral
      ! interaction: it keeps only the induced-dipole attraction and no
      ! repulsive core, so its friction weakens as T^-1 where a rigid core
      ! would hold it at T^-1/2.  It is therefore NOT called on its own --
      ! ion_neutral_pair_diffusion below adds it to the rigid-core channel,
      ! which is what keeps the friction from vanishing at high temperature.
      real*8, intent(in) :: TK, ntot, alpha_n, mu_g
      polarization_pair_diffusion = kb_erg*TK                              &
           /(2.21d0*pi*e_esu*ntot*sqrt(max(alpha_n*mu_g, 1.0d-60)))
      end function polarization_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function coulomb_pair_diffusion(TK, ntot, ne, Z_s, Z_t, mu_g)
      ! ION-ION.  Screened Coulomb, momentum-transfer cross section
      ! Q^(1) = 4 pi b_90^2 ln(Lambda) with b_90 = Z_s Z_t e^2/(mu g^2) the
      ! impact parameter of a 90 degree deflection.  Carrying it through the
      ! Chapman-Enskog Omega^(1,1) integral (int gam e^-gam^2 dgam = 1/2)
      ! gives
      !
      !    D = 3 (k T)^(5/2)
      !        / [ 4 (2 pi mu)^(1/2) n (Z_s Z_t e^2)^2 ln(Lambda) ]
      !
      ! the standard T^(5/2)/(n Z^2 Z^2 ln Lambda mu^(1/2)) Chapman-Enskog
      ! Coulomb result (Paquette et al. 1986, ApJS 61, 177, their section II;
      ! Schunk & Nagy eq. 4.142 give the equivalent collision frequency).
      ! Neither is in references/, so the numerical constant 3/(4 sqrt(2 pi))
      ! was re-derived here from the collision integral, not copied.
      real*8, intent(in) :: TK, ntot, ne, Z_s, Z_t, mu_g
      real*8 :: zz
      zz = abs(Z_s*Z_t)
      coulomb_pair_diffusion = 3.0d0*(kb_erg*TK)**2.5d0                    &
           /(4.0d0*sqrt(2.0d0*pi*mu_g)*ntot*(zz*e_esu*e_esu)**2            &
             *coulomb_logarithm(TK, ne, zz))
      end function coulomb_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function coulomb_logarithm(TK, ne, zz)
      ! ln(Lambda) = ln(lambda_D/b_90) with the Debye length screened by the
      ! electrons, lambda_D = (k T/(4 pi n_e e^2))^(1/2), and the 90 degree
      ! impact parameter evaluated at the thermal relative energy
      ! mu g^2 = 3 k T:
      !
      !    ln(Lambda) = ln( 3 k T lambda_D / (Z_s Z_t e^2) )
      !
      ! (Spitzer 1962, Physics of Fully Ionized Gases, section 5.2).  Floored
      ! at 1: the expansion the Coulomb cross section comes from needs
      ! Lambda >> 1, and a cell where it is not is one where the ion-ion pair
      ! carries no weight anyway.
      real*8, intent(in) :: TK, ne, zz
      real*8 :: lam_D
      lam_D = sqrt(kb_erg*TK/(4.0d0*pi*max(ne, 1.0d0)*e_esu*e_esu))
      coulomb_logarithm = max(log(3.0d0*kb_erg*TK*lam_D                    &
                                  /(max(zz, 1.0d0)*e_esu*e_esu)), 1.0d0)
      end function coulomb_logarithm

      ! ------------------------------------------------------------------ !

      real*8 function ion_neutral_pair_diffusion(TK, ntot, alpha_n,        &
                                                 A_s, A_t, mu_g)
      ! ION-NEUTRAL, NON-RESONANT: the induced-dipole (polarization) channel
      ! and the rigid-core channel taken TOGETHER,
      !
      !    1/D_in = 1/D_polarization + 1/D_hard sphere .
      !
      ! Reason: the two are momentum-transfer cross sections of the same
      ! encounter -- the long-range attraction and the short-range repulsion
      ! of one interaction potential -- and to first order their Q^(1) add,
      ! so their collision integrals add, so their FRICTIONS add.  Adding
      ! frictions is adding inverse diffusion coefficients (D = 3kT/(16 n mu
      ! Omega^(1,1)) is linear in 1/Omega), which is the same rule
      ! stage_mixture_diffusion uses across stages.  The physical content is
      ! that opening a second channel cannot make a pair MORE mobile: D_in is
      ! bounded above by the weaker-friction channel and can never exceed
      ! either limit.
      !
      ! Both limits are reproduced exactly.  D_pol ~ T and D_hs ~ T^(1/2), so
      ! their ratio crosses unity near 1.5e3 K for an H/He pair: below that
      ! the polarization term dominates the friction and D_in -> D_pol; above
      ! it the rigid core does and D_in -> D_hs.  At the 1e4 K of an ionized
      ! wind D_pol is 2.6 to 4.7 times D_hs, so the combined value sits within
      ! 20-30% of the hard sphere -- taking the polarization value alone there
      ! would understate the friction by a factor of a few and let an element
      ! settle correspondingly faster through the partially ionized layer,
      ! which is where ion-neutral pairs carry the friction at all.
      !
      ! The rigid-core channel uses the same Banks & Kockarts collision
      ! diameter as the neutral-neutral pair: no separate ion-neutral core
      ! radius is available, and the ionic radius differs from the atomic one
      ! by much less than the factor the combination is settling.
      real*8, intent(in) :: TK, ntot, alpha_n, A_s, A_t, mu_g
      real*8 :: Dpol, Dhs

      Dpol = polarization_pair_diffusion(TK, ntot, alpha_n, mu_g)
      Dhs  = hard_sphere_pair_diffusion(TK, ntot, A_s, A_t)
      ion_neutral_pair_diffusion = 1.0d0                                   &
           /(1.0d0/max(Dpol, 1.0d-99) + 1.0d0/max(Dhs, 1.0d-99))

      end function ion_neutral_pair_diffusion

      ! ------------------------------------------------------------------ !

      real*8 function stage_pair_diffusion(TK, ntot, ne, Z_s, m_s, al_s,   &
                                                       Z_t, m_t, al_t)
      ! The pair coefficient of two ionization stages, in whichever of the
      ! three limits above the pair belongs to.  Masses in m_H units,
      ! polarizabilities in cm^3 (used only when that partner is the neutral
      ! one).
      real*8, intent(in) :: TK, ntot, ne, Z_s, m_s, al_s, Z_t, m_t, al_t
      real*8 :: mu_g

      mu_g = m_s*m_t/(m_s + m_t)*mu
      if (Z_s .eq. 0.0d0 .and. Z_t .eq. 0.0d0) then
         stage_pair_diffusion = hard_sphere_pair_diffusion(TK, ntot,       &
                                                           m_s, m_t)
      else if (Z_s .ne. 0.0d0 .and. Z_t .ne. 0.0d0) then
         stage_pair_diffusion = coulomb_pair_diffusion(TK, ntot, ne,       &
                                                       Z_s, Z_t, mu_g)
      else if (Z_s .eq. 0.0d0) then
         stage_pair_diffusion = ion_neutral_pair_diffusion(TK, ntot,       &
                                                 al_s, m_s, m_t, mu_g)
      else
         stage_pair_diffusion = ion_neutral_pair_diffusion(TK, ntot,       &
                                                 al_t, m_s, m_t, mu_g)
      endif

      end function stage_pair_diffusion

      ! ------------------------------------------------------------------ !

      subroutine stage_mixture_diffusion(TK, ntot, ne, nsA, yA, ZA, mA,    &
                                         alA, nsB, yB, ZB, mB, alB,        &
                                         Deff, idom)
      ! Effective binary coefficient between two ELEMENTS whose ionization
      ! stages are resolved (memo eq. 11):
      !
      !    1/D_eff = sum_s sum_t y_s y_t / D_st
      !
      ! The element moves as one body (one velocity per element, not per
      ! stage), so its friction against the other element is the sum of the
      ! frictions its stages feel.  Frictions add, hence the INVERSE
      ! coefficients are averaged, weighted by the carrier fractions y.  This
      ! is Blanc's law, exact for a trace species in a mixture and the
      ! standard approximation otherwise; it is the same rule section 5 of
      ! the memo adopts across molecular carriers, applied across stages too.
      ! With one stage populated in each element it returns that pair exactly.
      !
      ! idom (optional) returns the packed index (s-1)*nsB + t of the pair
      ! carrying the largest share of the friction, for the diagnostic
      ! profile.
      integer,                         intent(in)  :: nsA, nsB
      real*8, dimension(1-Ng:N+Ng),    intent(in)  :: TK, ntot, ne
      real*8, dimension(1-Ng:N+Ng,nsA),intent(in)  :: yA
      real*8, dimension(1-Ng:N+Ng,nsB),intent(in)  :: yB
      real*8, dimension(nsA),          intent(in)  :: ZA, mA, alA
      real*8, dimension(nsB),          intent(in)  :: ZB, mB, alB
      real*8, dimension(1-Ng:N+Ng),    intent(out) :: Deff
      integer, dimension(1-Ng:N+Ng), optional, intent(out) :: idom

      real*8  :: fric, fpair, wgt, Dst, fbest
      integer :: j, is, it, kbest

      do j = 1-Ng, N+Ng
         fric  = 0.0d0
         fbest = 0.0d0
         kbest = 1
         do is = 1, nsA
            if (yA(j,is) .le. 0.0d0) cycle
            do it = 1, nsB
               if (yB(j,it) .le. 0.0d0) cycle
               wgt = yA(j,is)*yB(j,it)
               Dst = stage_pair_diffusion(TK(j), ntot(j), ne(j),           &
                                          ZA(is), mA(is), alA(is),         &
                                          ZB(it), mB(it), alB(it))
               fpair = wgt/max(Dst, 1.0d-99)
               fric  = fric + fpair
               if (fpair .gt. fbest) then
                  fbest = fpair
                  kbest = (is-1)*nsB + it
               endif
            enddo
         enddo
         if (fric .gt. 0.0d0) then
            Deff(j) = 1.0d0/fric
         else
            ! Neither element has a populated stage in this cell: the
            ! settling prefactor vanishes there anyway, and the gradient term
            ! only needs a finite coefficient.  Take the neutral pair.
            Deff(j) = hard_sphere_pair_diffusion(TK(j), ntot(j),           &
                                                 mA(1), mB(1))
         endif
         if (present(idom)) idom(j) = kbest
      enddo

      end subroutine stage_mixture_diffusion

      ! ------------------------------------------------------------------ !

      subroutine carrier_fractions(f_sp, ns, isp_list, yc)
      ! Fraction of an element's collision partners in each of its stages,
      ! sum_s y_s = 1.  f_sp holds particle counts per unit mass, so a ratio
      ! of f_sp entries is a ratio of particles -- which is what the friction
      ! weights, not a ratio of nuclei (an H2 is one collision partner).
      ! A cell with no carrier at all is given to the first (neutral) stage,
      ! so y is always normalized and D_eff always finite.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      integer,                                intent(in)  :: ns
      integer, dimension(ns),                 intent(in)  :: isp_list
      real*8, dimension(1-Ng:N+Ng,ns),        intent(out) :: yc

      real*8, dimension(1-Ng:N+Ng) :: tot
      integer :: is

      tot = 0.0d0
      do is = 1, ns
         yc(:,is) = max(f_sp(:,isp_list(is)), 0.0d0)
         tot      = tot + yc(:,is)
      enddo
      do is = 1, ns
         where (tot .gt. 0.0d0)
            yc(:,is) = yc(:,is)/tot
         elsewhere
            yc(:,is) = 0.0d0
         end where
      enddo
      where (tot .le. 0.0d0) yc(:,1) = 1.0d0

      end subroutine carrier_fractions

      ! ------------------------------------------------------------------ !

      subroutine helium_hydrogen_diffusion(rho, Tcode, f_sp, Deff, Dneut,  &
                                           idom)
      ! The stage-resolved binary coefficient between helium and component 1,
      ! which is what the transport operator carries as D_12.  Dneut is the
      ! all-neutral coefficient of the same cell, returned so the diagnostic
      ! and the reports can quote the suppression D_eff/D_neutral directly.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Deff
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: Dneut
      integer, dimension(1-Ng:N+Ng), optional, intent(out) :: idom

      real*8, dimension(1-Ng:N+Ng)        :: TK, ntot, carH, carHe
      real*8, dimension(1-Ng:N+Ng)        :: zb1, zbHe, ne_rel, ne_phys
      real*8, dimension(1-Ng:N+Ng,n_hecar):: yHe
      real*8, dimension(1-Ng:N+Ng,n_hcar) :: yH
      integer :: j

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      call carrier_counts(f_sp, carH, carHe)
      ntot = (carH + carHe)*rho*n0
      where (ntot .lt. 1.0d0) ntot = 1.0d0
      call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
      ne_phys = max(ne_rel*rho*n0, 1.0d0)

      call carrier_fractions(f_sp, n_hecar, hecar_isp, yHe)
      call carrier_fractions(f_sp, n_hcar,  hcar_isp,  yH)
      call stage_mixture_diffusion(TK, ntot, ne_phys,                      &
                                   n_hecar, yHe, hecar_Z, hecar_m,         &
                                   hecar_alpha,                            &
                                   n_hcar,  yH,  hcar_Z,  hcar_m,          &
                                   hcar_alpha, Deff, idom)

      if (present(Dneut)) then
         do j = 1-Ng, N+Ng
            Dneut(j) = hard_sphere_pair_diffusion(TK(j), ntot(j),          &
                                                  m_He_amu, m_H_amu)
         enddo
      endif

      end subroutine helium_hydrogen_diffusion

      ! ------------------------------------------------------------------ !

      subroutine mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
      ! Mean charge of each component per CARRIER and the free-electron
      ! density per unit mass.  The charge that enters the ambipolar force of
      ! the diffusion equation is the charge of a colliding particle, so the
      ! normalization is the carrier count, not the nucleus count: an H3+
      ! carries one charge for three hydrogen nuclei.  In the atomic region
      ! carriers and nuclei are the same particles and this is the charge per
      ! nucleus, as before.
      !
      ! The component charges use the species that carry ONE element only
      ! (HI/HII/H2/H2+/H3+ against HeI/HeII/HeIII/HeTR); HeH+, which carries
      ! both, is left out of both mean charges as a trace approximation and
      ! the metal charge is left out of component 1 (metal/H ~ 1e-4).  The
      ! electron density, by contrast, counts every charge under the same
      ! policy as utils/calc_ne, because it is the physical n_e the ambipolar
      ! field is built from.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: zb1, zbHe
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: ne_rel
      real*8, dimension(1-Ng:N+Ng) :: carH, carHe
      integer :: ib, im

      zb1    = 0.0d0
      zbHe   = 0.0d0
      ne_rel = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_charge(ib) .eq. 0) cycle
         if (bsp_nH(ib) .gt. 0 .and. bsp_nHe(ib) .eq. 0)                  &
            zb1  = zb1  + dble(bsp_charge(ib))*f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0 .and. bsp_nH(ib) .eq. 0)                  &
            zbHe = zbHe + dble(bsp_charge(ib))*f_sp(:,bsp_fsp(ib))
         ne_rel = ne_rel + dble(bsp_charge(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      if (eos_include_metals .and. thereis_metals) then
         do im = 1, n_mion
            if (mion_stage(im) .gt. 0)                                    &
               ne_rel = ne_rel + dble(mion_stage(im))*f_sp(:,mion_fsp(im))
         enddo
      endif
      call carrier_counts(f_sp, carH, carHe)
      zb1  = zb1 /max(carH,  1.0d-30)
      zbHe = zbHe/max(carHe, 1.0d-30)

      end subroutine mean_charges_and_electrons

      ! ------------------------------------------------------------------ !

      subroutine settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1,   &
                                      eEout, dlnpsi_out)
      ! The settling coefficient G of the header, and the gravity-normalized
      ! relative settling mass it comes from,
      !
      !   dmeff = (m_He - m_c1) - (Zbar_He - Zbar_c1) eE/(m_H g) [m_H units]
      !   G     = dmeff m_H g/(k T) + alpha_T dlnT/dr            [1/cm]
      !
      ! with CARRIER quantities: m_c1 = m_1 mcarH is the mass of component 1
      ! per collision partner (m_1 is its mass per hydrogen nucleus and mcarH
      ! the hydrogen nuclei bound into one carrier), and Zbar_c1, Zbar_He are
      ! the mean charges per carrier.  Gravity and the ambipolar field act on
      ! the particles that collide, so an H2 molecule weighs twice an H atom
      ! here even though it carries the same mass per nucleus (memo section
      ! 5).  In the atomic region mcarH = 1 and this is the nucleus form.
      !
      ! with eE = -k T dln(n_e T)/dr (central difference, one-sided at the
      ! ends).  dmeff is returned so the ambipolar limits can be tested
      ! against (m_He/m_H - 1) neutral, (m_He/m_H - 3/2) in an H+ plasma and
      ! (2 m_He/3 m_H - 1) in a He++ plasma (2.9715, 2.4715, 1.6477 with
      ! m_He/m_H = 3.9715) without a second definition of the field anywhere.  Where the local gravity
      ! vanishes dmeff is reported as zero (it is a ratio to g); G itself is
      ! always built from the forces, never from dmeff.  dmeff is the FORCE
      ! ratio only: the chemistry term -dln(psi)/dr that G also carries in the
      ! molecular region is not a mass and does not belong in it.
      !
      ! psi = n_1c/n_H is the number of collision partners component 1 spreads
      ! one hydrogen nucleus over: 1 in the atomic region, 1/2 where the
      ! hydrogen is H2.  Its gradient is the mole-fraction driver of the
      ! module header, and it is returned (dlnpsi_out) so that the trace-metal
      ! loop uses the one definition of it.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Gco, dmeff, zb1
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: eEout
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: dlnpsi_out

      real*8, dimension(1-Ng:N+Ng) :: carH, carHe, mcarH, zbHe, ne_rel
      real*8, dimension(1-Ng:N+Ng) :: ne_phys, mc1, lnmcar, dlnpsi
      real*8, dimension(1-Ng:N+Ng) :: TK, gphys, eEf, lnpe
      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, mass1, msum
      real*8 :: m_1, dr2
      integer :: j

      m_1 = mass_per_H_nucleus_without_He()
      call carrier_counts(f_sp, carH, carHe, mcarH)
      ! Mass of component 1 per hydrogen nucleus, from the cell's own species
      ! (mixture_mass_split) and not from the reservoir metal/H, so that the
      ! settling force and the mass the projection conserves are one number.
      ! Where a cell holds no hydrogen there is no component-1 carrier to put
      ! a force on and the reservoir value stands in.
      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      do j = 1-Ng, N+Ng
         if (nucH(j) .gt. 1.0d-30) then
            mc1(j) = (mass1(j)/nucH(j))*mcarH(j)
         else
            mc1(j) = m_1*mcarH(j)
         endif
      enddo
      call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      ne_phys = ne_rel*rho*n0
      do j = 1-Ng, N+Ng
         gphys(j) = Dphi(r(j))*v0*v0/R0                    ! [cm/s^2], inward
      enddo

      ! Ambipolar field eE = -k T dln(n_e T)/dr [erg/cm].
      eEf = 0.0d0
      if (he_ambipolar) then
         lnpe = log(max(ne_phys*TK, 1.0d-300))
         do j = 2-Ng, N+Ng-1
            dr2 = (r(j+1) - r(j-1))*R0
            eEf(j) = -kb_erg*TK(j)*(lnpe(j+1) - lnpe(j-1))/dr2
         enddo
         eEf(1-Ng) = eEf(2-Ng)
         eEf(N+Ng) = eEf(N+Ng-1)
      endif

      ! Mole-fraction driver of the molecular region, -dln(psi)/dr with
      ! psi = n_1c/n_H the collision partners per hydrogen nucleus (module
      ! header).  psi = 1 identically wherever the hydrogen is atomic, so the
      ! central difference below is exactly zero there.
      dlnpsi = 0.0d0
      lnmcar = log(max(mcarH, 1.0d-30))      ! ln psi = -ln mcarH
      do j = 2-Ng, N+Ng-1
         dlnpsi(j) = -(lnmcar(j+1) - lnmcar(j-1))                         &
                     /max((r(j+1) - r(j-1))*R0, 1.0d0)
      enddo
      dlnpsi(1-Ng) = dlnpsi(2-Ng)
      dlnpsi(N+Ng) = dlnpsi(N+Ng-1)

      do j = 1-Ng, N+Ng
         Gco(j) = ((m_He_amu - mc1(j))*mu*gphys(j)                   &
                   - (zbHe(j) - zb1(j))*eEf(j))/(kb_erg*TK(j))            &
                  - dlnpsi(j)
         if (abs(gphys(j)) .gt. 0.0d0) then
            dmeff(j) = (m_He_amu - mc1(j))                                &
                       - (zbHe(j) - zb1(j))*eEf(j)/(mu*gphys(j))
         else
            dmeff(j) = 0.0d0
         endif
      enddo

      ! Thermal diffusion: alpha_T dlnT/dr, default he_alphaT = 0 (no-op).
      if (he_alphaT .ne. 0.0d0) then
         do j = 2-Ng, N+Ng-1
            Gco(j) = Gco(j) + he_alphaT*(log(TK(j+1)) - log(TK(j-1)))     &
                     / max((r(j+1)-r(j-1))*R0, 1.0d0)
         enddo
      endif

      ! The field itself, so that every element settling against component 1
      ! -- helium here, the trace metals below -- uses ONE ambipolar field
      ! and there is no second, cruder copy of it anywhere in the module.
      if (present(eEout)) eEout = eEf
      if (present(dlnpsi_out)) dlnpsi_out = dlnpsi

      end subroutine settling_coefficient

      ! ------------------------------------------------------------------ !

      function relative_settling_mass(rho, Tcode, f_sp) result(dmeff)
      ! Diagnostic: the gravity-normalized relative settling mass of helium
      ! against component 1, in m_H units: m_He/m_H - 1 neutral,
      ! m_He/m_H - 3/2 in an H+ plasma, 2 m_He/(3 m_H) - 1 in a fully ionized
      ! He++ plasma (2.9715, 2.4715, 1.6477).  Same definition the operator
      ! uses, there is no second copy of the ambipolar field.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng) :: dmeff
      real*8, dimension(1-Ng:N+Ng) :: Gco, zb1

      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1)

      end function relative_settling_mass

      ! ------------------------------------------------------------------ !

      subroutine drift_and_gradient_face_coefficients(rho_phys, Dco, Gco,  &
                                                      rp, Agrd, Bdrf, updrf)
      ! Face coefficients of the diffusive flux at r_edg(0:N).  The flux
      ! (module header) is
      !
      !     J(f) = -Agrd(f) (X(j+1) - X(j))  -  Bdrf(f) [X(1-X)](f)
      !
      ! with Agrd = rho (D_12 + K_zz)/dr already divided by the face spacing
      ! and Bdrf = rho D_12 G the drift coefficient.  Neither depends on X:
      ! the whole composition dependence of the drift sits in the product
      ! X(1-X), which element_face_flux evaluates at the new level.
      !
      ! updrf selects how that product is taken at the face:
      !    0  central, X_f (1 - X_f) with X_f the face average.  Used where
      !       the drift is resolved, |Bdrf| dr <= 2 (A + E), which is also
      !       where the gradient flux dominates it -- the condition that
      !       keeps the central branch inside the bounds (header).
      !   -1  Bdrf >= 0: helium drifts INWARD, so cell j+1 donates helium and
      !       cell j donates the hydrogen that moves the other way.
      !   +1  Bdrf <  0: helium drifts outward, the roles exchanged.
      !
      ! Faces 0 and N are the two boundaries and carry zero diffusive flux.
      ! These are the DIFFUSIVE face coefficients only; the advection is the
      ! divergence of the hydro's face mass fluxes and is carried in the
      ! Runge-Kutta stages (see the header).  The outer ghost carries the top
      ! cell's own composition, so an inflowing top boundary brings in gas of
      ! the same composition instead of a silent zero flux.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rho_phys, Dco, Gco
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: rp
      real*8, dimension(0:N),       intent(out) :: Agrd, Bdrf
      integer, dimension(0:N),      intent(out) :: updrf

      real*8 :: dr_f, rhof, Df, Gf, Kf
      integer :: j

      Agrd  = 0.0d0
      Bdrf  = 0.0d0
      updrf = 0

      do j = 1, N-1
         dr_f = max(rp(j+1)-rp(j), 1.0d0)
         rhof = 0.5d0*(rho_phys(j)+rho_phys(j+1))
         Df   = 0.5d0*(Dco(j)+Dco(j+1))
         Gf   = 0.5d0*(Gco(j)+Gco(j+1))
         Kf   = 0.5d0*(kzz_cell(j)+kzz_cell(j+1))
         Agrd(j) = rhof*(Df + Kf)/dr_f            ! gradient + eddy, >= 0
         Bdrf(j) = rhof*Df*Gf                     ! drift, either sign
         if (abs(Bdrf(j))*dr_f .le. 2.0d0*rhof*(Df + Kf)) then
            updrf(j) = 0
         else if (Bdrf(j) .ge. 0.0d0) then
            updrf(j) = -1
         else
            updrf(j) = +1
         endif
      enddo

      end subroutine drift_and_gradient_face_coefficients

      ! ------------------------------------------------------------------ !

      subroutine element_face_flux(Xl, Xr, Agr, Bst, upw, Jf, dJl, dJr)
      ! The helium mass flux through one face and its two slopes,
      !   Jf = -Agr (Xr - Xl) - Bst [X(1-X)](f),
      !   dJl = dJf/dXl,  dJr = dJf/dXr,
      ! in the branch drift_and_gradient_face_coefficients selected.  The two
      ! factors of the drift product come from opposite sides of the face:
      ! helium from the cell it leaves, hydrogen from the cell it enters,
      ! which is where the bounds on X come from (module header).
      !
      ! Slopes: dJl >= 0 and dJr <= 0 for any Xl, Xr in [0,1], in all three
      ! branches, which is the M-matrix condition on the Newton Jacobian.
      real*8,  intent(in)  :: Xl, Xr, Agr, Bst
      integer, intent(in)  :: upw
      real*8,  intent(out) :: Jf, dJl, dJr

      real*8 :: Xf

      Jf  = -Agr*(Xr - Xl)
      dJl =  Agr
      dJr = -Agr
      if (upw .eq. 0) then                     ! central, drift resolved
         Xf  = 0.5d0*(Xl + Xr)
         Jf  = Jf  - Bst*Xf*(1.0d0 - Xf)
         dJl = dJl - 0.5d0*Bst*(1.0d0 - 2.0d0*Xf)
         dJr = dJr - 0.5d0*Bst*(1.0d0 - 2.0d0*Xf)
      else if (upw .lt. 0) then                ! inward: He from j+1, H from j
         Jf  = Jf  - Bst*Xr*(1.0d0 - Xl)
         dJl = dJl + Bst*Xr
         dJr = dJr - Bst*(1.0d0 - Xl)
      else                                     ! outward: He from j, H from j+1
         Jf  = Jf  - Bst*Xl*(1.0d0 - Xr)
         dJl = dJl - Bst*(1.0d0 - Xr)
         dJr = dJr + Bst*Xl
      endif

      end subroutine element_face_flux

      ! ------------------------------------------------------------------ !

      subroutine element_advective_divergence(Y, Frho, dvF, dvM, Yf_out,   &
                                              Fs_out)
      ! THE MATERIAL ADVECTION OF ONE ELEMENT, as the divergence of the face
      ! species fluxes the hydrodynamic stages carry:
      !
      !     F_Y(j) = F_rho(j) Y^face(j),
      !     dvF(j) = ( A_+ F_Y(j) - A_- F_Y(j-1) ) / V_j ,
      !
      ! with F_rho the face mass flux of the mass row of this state, Y^face
      ! the reconstructed and upwinded face mass fraction of
      ! species_face_fraction, and A, V the face area and the exact shell
      ! volume of spherical_face_area_and_cell_volume, which the diffusive
      ! half of the same row divides by as well.  dvM is the same expression
      ! with both fluxes in magnitude, for a row scale.  Both are in code
      ! units; the caller multiplies by the factor that returns them to the
      ! units of its own row.
      !
      ! THERE IS ONE OF THESE IN THE ELEMENT OPERATOR, and the backward-Euler
      ! step, the fixed-wind relaxation and the stationary elemental row all
      ! read it here.  Two spellings of one term have two fixed points, so an
      ! alternation of a relaxation with one and a row with the other cannot
      ! converge to a state the row reads as stationary.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Y, Frho
      real*8, dimension(1-Ng:N+Ng), intent(out) :: dvF, dvM
      ! The face mass fraction and the face element mass flux the divergence
      ! was formed from (optional), for a record of the row term by term.
      real*8, dimension(1-Ng:N+Ng), intent(out), optional :: Yf_out, Fs_out

      real*8, dimension(1-Ng:N+Ng) :: Yf, Fs
      real*8, dimension(0:N)       :: fa
      real*8, dimension(1:N)       :: cv
      integer :: j

      call species_face_fraction(Y, Frho, Yf)
      call species_face_flux(Frho, Yf, Fs)
      call spherical_face_area_and_cell_volume(fa, cv)

      dvF = 0.0d0
      dvM = 0.0d0
      do j = 1, N
         dvF(j) = (fa(j)*Fs(j) - fa(j-1)*Fs(j-1))/cv(j)
         dvM(j) = (fa(j)*abs(Fs(j)) + fa(j-1)*abs(Fs(j-1)))/cv(j)
      enddo
      if (present(Yf_out)) Yf_out = Yf
      if (present(Fs_out)) Fs_out = Fs

      end subroutine element_advective_divergence

      ! ------------------------------------------------------------------ !

      subroutine element_advective_face_coefficients(Frho, advj, advm)
      ! The two face coefficients of that divergence, in code units per unit
      ! of the face MASS fraction:
      !
      !     advj(j) =  A_+ F_rho(j)  /V_j ,
      !     advm(j) = -A_- F_rho(j-1)/V_j ,
      !
      ! the same A_+, A_- and V_j the residual divides by.  They are
      ! what the tridiagonal Newton of the implicit step differentiates the
      ! term with, taking the face composition at its donor cell's own value:
      ! the limiter and the second-order part of the reconstruction are left
      ! out of the Jacobian, because they are the strongly nonlinear part of
      ! the operator and reach beyond the tridiagonal band.  The residual is
      ! the full operator either way, so this changes what the Newton
      ! converges AT, not what it converges TO.
      !
      ! The sign pattern is the M-matrix one the diffusive half already has:
      ! an outflowing face puts a nonnegative entry on the diagonal and a
      ! nonpositive one on its donor neighbor, an inflowing face the other
      ! way round, so no off-diagonal entry is positive.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Frho
      real*8, dimension(1:N),       intent(out) :: advj, advm
      real*8, dimension(0:N) :: fa
      real*8, dimension(1:N) :: cv
      integer :: j
      call spherical_face_area_and_cell_volume(fa, cv)
      do j = 1, N
         advj(j) =  fa(j)  *Frho(j)  /cv(j)
         advm(j) = -fa(j-1)*Frho(j-1)/cv(j)
      enddo
      end subroutine element_advective_face_coefficients

      ! ------------------------------------------------------------------ !

      subroutine solve_mass_fraction(Xold, Xhe, rho_phys, dt_phys,          &
                                     Frho, cadvf, Agrd, Bdrf, updrf,        &
                                     X_base, jlo, advect, solved)
      ! One implicit (backward-Euler) step of
      !   rho (X^new - X^old)/dt + div(r^2 J)/r^2 + advect div(F_rho X) = 0
      ! for rows jlo..N.  With advect false the row is the diffusive half
      ! alone, the advection having been carried on the hydro's face mass
      ! fluxes inside the Runge-Kutta stages; with it true the advective term
      ! is the divergence of those same face fluxes, formed here by
      ! element_advective_divergence, for the fixed-wind relaxation, which
      ! has no stages to ride on.  cadvf returns that code divergence to the
      ! [g cm^-3 s^-1] of the row.  The step is NONLINEAR in X^new,
      ! because the drift flux carries the product X(1-X) at the new level,
      ! and it is solved by Newton: each iteration assembles the residual and
      ! its tridiagonal Jacobian and solves for the correction.
      !
      ! Newton, not a Picard sweep on a lagged (1-X): the lag is exactly what
      ! removes the shutoff of the drift at the ends of the composition axis,
      ! and with it the bounds on X (header, and docs/Update_EXHALE_stage1.pdf
      ! section 86).  The Jacobian is an M-matrix for any iterate in [0,1],
      ! and the nonlinearity is quadratic, so the iteration converges in a
      ! few passes; a step that does not reduce the residual is halved.
      !
      ! Xold is the state at the start of the step and Xhe carries the initial
      ! iterate in and the solution out.
      !
      ! WHAT COMES BACK IN Xhe IS THE LAST ITERATE THE ACCEPTANCE RULE BELOW
      ! ACCEPTED, and solved says whether that iterate solves the step.  A
      ! trial the line search could not make better than the iterate it
      ! started from is not a Newton step of this equation: it is a direction
      ! the residual does not decrease along, so adopting it moves the
      ! composition without reference to the row.  It is discarded, the entry
      ! iterate stands, and solved is false.  A NaN or an infinity anywhere in
      ! a trial fails the same comparison and is discarded by the same rule.
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: Xold, rho_phys, dt_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: Frho, cadvf
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: Xhe
      real*8, dimension(0:N),       intent(in)    :: Agrd, Bdrf
      integer, dimension(0:N),      intent(in)    :: updrf
      real*8,                       intent(in)    :: X_base
      integer,                      intent(in)    :: jlo
      logical,                      intent(in)    :: advect
      logical,                      intent(out)   :: solved

      real*8, dimension(1-Ng:N+Ng) :: aa, bb, cc, dd, cpv, dpv, dX, Xtry
      real*8, dimension(0:N)       :: Jf, dJl, dJr
      real*8, dimension(1:N)       :: advj, advm
      ! The one geometry of the grid, in cm: A(f) R0^2 and V(j) R0^3.
      real*8, dimension(0:N)       :: fa
      real*8, dimension(1:N)       :: cv
      real*8  :: R0sq, R0cb
      real*8  :: Kj, mden, sL, sR, cadv, rnorm, rprev, rtry, damp, rstart
      integer :: j, it, ihalf, nit
      logical :: descended
      logical, save :: warned_newton = .false.
      ! Convergence.  The residual is measured RELATIVE to the size of the
      ! terms that make it up, because those terms cancel: at the base the
      ! step is ~1e5 cell diffusion times, so the time term and the flux
      ! divergence are each ~1e5 times their difference and an absolute
      ! residual cannot go below the round-off of the larger one.  The
      ! relative residual reaches ~1e-12 in a few passes; an absolute one in X
      ! stalls at ~1e-11 and burns the iteration limit for nothing.  Below
      ! that the round-off of the tridiagonal solve, amplified by the same
      ! conditioning, is what is left, so an already-small residual that no
      ! longer halves is taken as converged -- Newton is quadratic here, so a
      ! pass that does not gain a factor 2 below newton_floor has reached the
      ! arithmetic floor and further passes only move round-off.  The floor
      ! stall test fires only once the residual is ALREADY small, in one of
      ! two senses: below newton_floor, or newton_drop below the residual the
      ! step started at.  Both are needed because the reachable floor is set
      ! by the conditioning of the tridiagonal solve and varies by orders of
      ! magnitude across the cases here -- 1e-12 in a wind, ~1e-5 in the
      ! K_zz = 2e12 homopause column of test T7, where the time term is
      ! negligible against the eddy term.  Stopping on stagnation ALONE stops
      ! too early: the damped passes of the strong-drift column of T14 stall
      ! twice at a residual of order 1 and then converge, and a step stopped
      ! there loses the bounds on X, which belong to the CONVERGED step (T14
      ! then fails at an excursion of 7e-3).
      integer, parameter :: newton_maxit  = 30
      real*8,  parameter :: newton_tol    = 1.0d-12
      real*8,  parameter :: newton_floor  = 1.0d-8
      ! Set to clear the ~1e-5 floor the paragraph above documents for the
      ! K_zz = 2e12 homopause column of test T7 (MEASURED here at 2e-6 to
      ! 1.1e-5 once the no-descent branch below reads the same floor a
      ! descending pass does): 1e-6 left that column's own floor outside
      ! both floor tests, so a column at its OWN conditioning limit read as
      ! unsolved on every pass and never advanced.
      real*8,  parameter :: newton_drop   = 2.0d-5
      integer, parameter :: newton_halves = 8
      ! A residual this large after every pass is not a conditioning floor,
      ! it is a step that was not solved; that is what gets announced.
      real*8,  parameter :: newton_unsolved = 1.0d-4

      if (jlo .eq. 2) Xhe(1-Ng:1) = X_base
      Xhe(N+1:N+Ng) = Xhe(N)

      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0

      advj = 0.0d0
      advm = 0.0d0
      if (advect) call element_advective_face_coefficients(Frho, advj, advm)

      call composition_residual(Xold, Xhe, rho_phys, dt_phys,             &
                                Frho, cadvf, Agrd, Bdrf, updrf, jlo,      &
                                advect, dd, Jf, dJl, dJr, rnorm)

      rstart   = rnorm
      nit      = 0
      solved   = .false.
      do it = 1, newton_maxit
         if (rnorm .le. newton_tol) then
            solved = .true.               ! converged on its own tolerance
            exit
         endif
         ! --- Jacobian of the residual, tridiagonal by construction
         do j = jlo, N
            Kj    = 1.0d0/(cv(j)*R0cb)
            sL    = fa(j-1)*R0sq
            sR    = fa(j)*R0sq
            aa(j) = -Kj*sL*dJl(j-1)
            bb(j) =  rho_phys(j)/dt_phys(j)                               &
                   + Kj*(sR*dJl(j) - sL*dJr(j-1))
            cc(j) =  Kj*sR*dJr(j)
            ! DONOR-CELL LINEARIZATION OF THE FACE-FLUX DIVERGENCE.  The term
            ! of cell j is
            !
            !   adv(j) = [ A_+ F_rho(j) X(don(j))
            !              - A_- F_rho(j-1) X(don(j-1)) ] / V_j x cadvf(j),
            !
            ! with don(f) = f where F_rho(f) >= 0 and f+1 where it is
            ! negative, so each face contributes one entry, on the diagonal
            ! or on the neighbor the face draws from.  The coefficients do
            ! NOT sum to zero: their sum is div(F_rho) cadvf(j), which is the
            ! mass row's own imbalance, and an element inherits it.
            if (advect) then
               cadv = advj(j)*cadvf(j)
               if (Frho(j) .ge. 0.0d0) then
                  bb(j) = bb(j) + cadv
               else if (j .lt. N) then
                  cc(j) = cc(j) + cadv
               else
                  ! The outer ghost is zero-gradient, X(N+1) = X(N), so the
                  ! entry of an inflowing top face belongs on the diagonal.
                  bb(j) = bb(j) + cadv
               endif
               cadv = advm(j)*cadvf(j)
               if (j .eq. jlo) then
                  ! Cell jlo-1 is data and not an unknown: the Dirichlet
                  ! reservoir, or the ghost of a closed base.  Only an
                  ! inflowing inner face, whose donor is cell j itself,
                  ! leaves an entry.
                  if (Frho(j-1) .lt. 0.0d0) bb(j) = bb(j) + cadv
               else if (Frho(j-1) .ge. 0.0d0) then
                  aa(j) = aa(j) + cadv
               else
                  bb(j) = bb(j) + cadv
               endif
            endif
         enddo
         cc(N) = 0.0d0
         if (jlo .eq. 2) aa(2) = 0.0d0        ! X(1) is the Dirichlet base

         ! --- solve J dX = -residual (dd already carries -residual)
         cpv(jlo) = cc(jlo)/bb(jlo)
         dpv(jlo) = dd(jlo)/bb(jlo)
         do j = jlo+1, N
            mden   = bb(j) - aa(j)*cpv(j-1)
            cpv(j) = cc(j)/mden
            dpv(j) = (dd(j) - aa(j)*dpv(j-1))/mden
         enddo
         dX(N) = dpv(N)
         do j = N-1, jlo, -1
            dX(j) = dpv(j) - cpv(j)*dX(j+1)
         enddo

         ! --- accept the step, halving it while it does not reduce the
         ! residual (the nonlinearity is quadratic, so this is rarely used).
         ! The acceptance rule is descent and nothing else: the shortest
         ! trial is not accepted for being the shortest.
         damp      = 1.0d0
         descended = .false.
         do ihalf = 0, newton_halves
            Xtry = Xhe
            Xtry(jlo:N) = Xhe(jlo:N) + damp*dX(jlo:N)
            if (jlo .eq. 2) Xtry(1-Ng:1) = X_base
            Xtry(N+1:N+Ng) = Xtry(N)
            call composition_residual(Xold, Xtry, rho_phys, dt_phys,      &
                                      Frho, cadvf, Agrd, Bdrf,            &
                                      updrf, jlo, advect, dd,             &
                                      Jf, dJl, dJr, rtry)
            if (rtry .lt. rnorm) then
               descended = .true.
               exit
            endif
            if (ihalf .eq. newton_halves) exit
            damp = 0.5d0*damp
         enddo
         ! Every halving spent without descent.  The direction is not one the
         ! residual of this step decreases along, so no trial along it is an
         ! iterate of this equation: the trial is discarded and Xhe is left
         ! where the last accepted one put it.  A residual that already
         ! stands at the arithmetic floor (newton_floor) or has already
         ! dropped newton_drop below where the step started is the SAME
         ! floor the it >= 3 branch above recognizes on a descending pass;
         ! a step whose last descent reached it and then found no further
         ! trial to improve on has reached that floor too, and is not an
         ! unsolved row for having reached it one pass sooner.  Without this
         ! branch the floor was only ever read on a DESCENDING pass, so a
         ! step whose residual had already fallen to it (rnorm ~ 1e-12 to
         ! 1e-9, MEASURED on the closed-column relaxations of T1a/T7c/T7d/T10)
         ! reported unsolved on the very next pass, whose halving could not
         ! improve on round-off: not a failed step, an already-solved one.
         if (.not. descended) then
            if (rnorm .le. newton_floor .or. rnorm .le. newton_drop*rstart)&
               solved = .true.
            exit
         endif
         Xhe   = Xtry
         rprev = rnorm
         rnorm = rtry
         nit   = it
         if (it .ge. 3 .and. rnorm .gt. 0.5d0*rprev .and.                 &
             (rnorm .le. newton_floor .or.                                &
              rnorm .le. newton_drop*rstart)) then
            solved = .true.        ! the arithmetic floor of the header
            exit
         endif
      enddo

      ! Neither converged nor stalled at the floor: the iteration used up
      ! every pass.  Such a step carries the equation only in so far as its
      ! residual is a conditioning floor rather than an unsolved row, which
      ! is the separation newton_unsolved makes.
      if (.not. solved .and. nit .ge. newton_maxit)                       &
         solved = (rnorm .le. newton_unsolved)

      he_fraction_newton_steps = nit
      he_fraction_newton_resid = rnorm

      ! The bounds on X belong to the CONVERGED step (header), so a step that
      ! used up every pass while still making progress -- neither converged
      ! nor stalled -- is the one case in which the range clip could have work
      ! to do.  It has not been seen in any run measured; it is announced once
      ! if it happens, whether or not the step diagnostic is on.
      if (nit .ge. newton_maxit .and. rnorm .gt. newton_unsolved .and.    &
          .not. warned_newton) then
         warned_newton = .true.
         write(*,'(A,I0,A,ES10.3,A)') ' (element diffusion) WARNING: the'// &
              ' composition step did not converge in ', newton_maxit,     &
              ' Newton passes (relative residual ', rnorm,                &
              '); the bounds on X are not guaranteed for it.'
      endif

      end subroutine solve_mass_fraction

      ! ------------------------------------------------------------------ !

      subroutine composition_residual(Xold, Xhe, rho_phys, dt_phys,       &
                                      Frho, cadvf, Agrd, Bdrf,            &
                                      updrf, jlo, advect,                 &
                                      mres, Jf, dJl, dJr, rnorm,          &
                                      res_out, dsc_out, faces_out,        &
                                      Yface_out)
      ! Residual of the implicit composition step, returned NEGATED (mres is
      ! the right-hand side of the Newton system), together with the face
      ! fluxes, their slopes, and the residual measured RELATIVE to the size
      ! of the terms of its own row (dsc, the same terms with their absolute
      ! values): the terms cancel to many digits at a step long against the
      ! cell diffusion time, so an absolute residual is a measure of their
      ! round-off and not of convergence.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Xold, Xhe, rho_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: dt_phys
      ! The face mass flux the element rides on, and the factor that returns
      ! the code divergence of its face element flux to the [g cm^-3 s^-1]
      ! this row is written in.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Frho, cadvf
      real*8, dimension(0:N),       intent(in)  :: Agrd, Bdrf
      integer, dimension(0:N),      intent(in)  :: updrf
      integer,                      intent(in)  :: jlo
      ! Include the material advection div(F_rho X).  False for the marching
      ! operator, whose advection is that same divergence taken inside the
      ! Runge-Kutta stages; true for the stationary balance and the
      ! fixed-wind relaxation, which have no stages to ride on.
      logical,                      intent(in)  :: advect
      ! res_out and dsc_out (optional) return the SAME residual and the SAME
      ! scale cell by cell, unsigned and unnegated, for a caller that has to
      ! name the cell a balance is out at rather than only the largest value
      ! of the column.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: mres
      real*8, dimension(0:N),       intent(out) :: Jf, dJl, dJr
      real*8,                       intent(out) :: rnorm
      real*8, dimension(1-Ng:N+Ng), intent(out), optional :: res_out, dsc_out
      ! faces_out (optional) returns the four face terms of each row, each
      ! with its face area and divided by the cell volume, signed as they
      ! enter the row and in the row's units: (:,1) the diffusive flux of the
      ! inner face j-1, (:,2) that of the outer face j, (:,3) and (:,4) the
      ! same for the advective face element mass flux (zero with advect
      ! false).  Their sum is the row less its time term, up to the rounding
      ! of the regrouping.  Yface_out (optional) returns the face mass
      ! fraction the advective half carried.  Both are records of the terms
      ! and change nothing the row is formed from.
      real*8, dimension(1-Ng:N+Ng,4), intent(out), optional :: faces_out
      real*8, dimension(1-Ng:N+Ng),   intent(out), optional :: Yface_out

      real*8, dimension(1-Ng:N+Ng) :: dvF, dvM, Yf, Fs
      ! The one geometry of the grid, in cm: A(f) R0^2 and V(j) R0^3.
      real*8, dimension(0:N)       :: fa
      real*8, dimension(1:N)       :: cv
      real*8  :: R0sq, R0cb
      real*8  :: Kj, sL, sR, res, dsc, q
      integer :: j
      logical :: nonfinite_row

      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0
      nonfinite_row = .false.

      do j = 0, N
         call element_face_flux(Xhe(j), Xhe(j+1), Agrd(j), Bdrf(j),       &
                                updrf(j), Jf(j), dJl(j), dJr(j))
      enddo

      dvF = 0.0d0
      dvM = 0.0d0
      Yf  = 0.0d0
      Fs  = 0.0d0
      if (advect) call element_advective_divergence(Xhe, Frho, dvF, dvM,  &
                                                    Yf, Fs)
      if (present(Yface_out)) Yface_out = Yf

      mres  = 0.0d0
      rnorm = 0.0d0
      if (present(res_out)) res_out = 0.0d0
      if (present(dsc_out)) dsc_out = 0.0d0
      if (present(faces_out)) faces_out = 0.0d0
      do j = jlo, N
         Kj  = 1.0d0/(cv(j)*R0cb)
         sL  = fa(j-1)*R0sq
         sR  = fa(j)*R0sq
         res = rho_phys(j)*(Xhe(j) - Xold(j))/dt_phys(j)                  &
             + Kj*(sR*Jf(j) - sL*Jf(j-1))
         dsc = rho_phys(j)*max(abs(Xhe(j)), abs(Xold(j)))/dt_phys(j)      &
             + Kj*(sR*abs(Jf(j)) + sL*abs(Jf(j-1)))
         if (advect) then
            res = res + dvF(j)*cadvf(j)
            dsc = dsc + dvM(j)*cadvf(j)
         endif
         mres(j) = -res
         if (present(res_out)) res_out(j) = res
         if (present(dsc_out)) dsc_out(j) = dsc
         if (present(faces_out)) then
            faces_out(j,1) = -Kj*sL*Jf(j-1)
            faces_out(j,2) =  Kj*sR*Jf(j)
            if (advect) then
               faces_out(j,3) = -fa(j-1)*Fs(j-1)/cv(j)*cadvf(j)
               faces_out(j,4) =  fa(j)*Fs(j)/cv(j)*cadvf(j)
            endif
         endif
         q = abs(res)/max(dsc, 1.0d-300)
         if (q .le. huge(1.0d0)) then
            rnorm = max(rnorm, q)
         else
            nonfinite_row = .true.
         endif
      enddo
      ! A ROW THAT IS NOT AN ORDINARY NUMBER MAKES THE NORM NO MEASUREMENT.
      ! MAX with a NaN argument is processor dependent (the SSE maxsd the
      ! compiler emits returns whichever operand the code generation put
      ! second), so a NaN row could drop out of a max-reduction and the
      ! Newton iteration of solve_mass_fraction read an unsolved step as
      ! converged.  The test is written out instead: a NaN or infinite row
      ! sets the norm to huge, which no trial descends below and no
      ! convergence test accepts.  Finite rows reduce exactly as before.
      if (nonfinite_row) rnorm = huge(1.0d0)

      end subroutine composition_residual

      ! ------------------------------------------------------------------ !

      subroutine project_elements(f_sp, Xhe, msum, scale_metals)
      ! Nonnegative projection of the species vector onto the new element
      ! totals (memo section 3).  Helium-bearing species are scaled by
      ! r_He = n_He^new/n_He^old, every component-1 species by r_1, and a
      ! species carrying a nucleus of each (HeH+) by min(r_1, r_He); the
      ! nuclei this under-counts for the element with the larger factor are
      ! deposited into that element's neutral ground species.  Every
      ! operation is a multiplication by a nonnegative factor or an addition,
      ! so no negative intermediate can arise, and both element totals are
      ! met exactly.  The result is the initial guess for ioniz_eq, which
      ! owns the split WITHIN an element; this step owns the element totals.
      !
      ! THE INVARIANT: THE SPECIES HANDED BACK CARRY THE MASS msum THE
      ! SPECIES HANDED IN CARRIED, to round-off.  The two components are
      ! given mass1 = (1 - X) msum and m_He n_He = X msum, and component 1
      ! is scaled AS A WHOLE, so its mass ends at (1 - X) msum whatever it is
      ! made of: r_1 is the ratio of the two component-1 MASSES and not
      ! n_H^new/n_H^old computed from a reservoir mass per hydrogen nucleus.
      ! Writing it the second way is what let the projection lose mass -- the
      ! reservoir m_1 = mass_per_H_nucleus_without_He() carries the metals at
      ! melem_ab, the cell carries them at whatever the trace transport and
      ! the advection of the element rows left, and the difference between
      ! the two was written into hydrogen at every call (6.7e-3 of the
      ! density at the top of an atomic wind, MEASURED).  The mass a cell's
      ! species carry is read from the cell (mixture_mass_split).
      !
      ! scale_metals says whether the metal ions belong to the scaled group.
      ! They do on the path that transports helium alone and carries the
      ! metals with hydrogen (element_diffusion_step, relax_element_
      ! composition).  They do not where each metal element has its own
      ! advected mass fraction and has already been written at it
      ! (project_element_mass_fractions): there the metals are part of the
      ! mass already spoken for, and what closes the mixture is the hydrogen
      ! group alone.
      !
      ! THE METALS ARE A RATIO, NOT A DENSITY, on the path that scales them.
      ! There component 1 carries the trace metals at the metal/H the cell
      ! holds, so what this projection has to preserve for them is n_X/n_H
      ! and not n_X (r_1 is the component-1 mass ratio, and it is also
      ! n_H^new/n_H^old, so one factor does both).  Multiplying by r_H
      ! does exactly that -- in a cell that HAD hydrogen.  Where the cell had
      ! none the ratio it carries is 0/0, and multiplying by r_H = 0 destroys
      ! it: hydrogen comes back through the shortfall deposit below and the
      ! metals do not, so the cell keeps zero metals for good (every later call
      ! multiplies that zero by something).  A cell emptied of carbon between
      ! two cells at the reservoir C/H is not a state of the atmosphere -- these
      ! equations have no sink for elemental carbon -- and it removes the C I
      ! cooling that sets the temperature there.  The metals therefore return
      ! WITH the hydrogen, at melem_ab, which is the reservoir metal/H that
      ! m_1 = mass_per_H_nucleus_without_He() is built from and that the
      ! hydrogen count of the same cell was taken against just above; it is
      ! also the ratio set_IC and load_IC build an absent element from.  All of it goes into the neutral ground stage, as those two do and
      ! as the trace-metal re-seed of element_diffusion_step does, because a
      ! cell that held no hydrogen held no ionization split either and ioniz_eq
      ! re-solves the split from the element total on the next call.
      ! Section 84 of docs/Update_EXHALE_stage1.pdf.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: Xhe, msum
      logical,                                intent(in)    :: scale_metals

      real*8, dimension(1-Ng:N+Ng) :: nucH_old, nucHe_old, mass1, mnow
      real*8, dimension(1-Ng:N+Ng) :: mmetal
      real*8  :: nucH_new, nucHe_new, rH, rHe, rBoth, gotH, gotHe
      real*8  :: grp_old, grp_new, m1c, m_1
      integer :: j, ib, im, ie, i0e, k

      m_1 = mass_per_H_nucleus_without_He()
      call mixture_mass_split(f_sp, nucH_old, nucHe_old, mass1, mnow,     &
                              mmetal)

      do j = 1-Ng, N+Ng
         ! The mass each component is to end with, and the mass the group
         ! that is about to be scaled starts with.
         nucHe_new = Xhe(j)*msum(j)/m_He_amu
         grp_new   = (1.0d0 - Xhe(j))*msum(j)
         grp_old   = mass1(j)
         if (.not. scale_metals) then
            grp_new = grp_new - mmetal(j)
            grp_old = grp_old - mmetal(j)
         endif
         if (grp_new .lt. 0.0d0) grp_new = 0.0d0
         if (nucH_old(j) .gt. 1.0d-30 .and. grp_old .gt. 1.0d-30) then
            ! The cell's own mass per hydrogen nucleus of the scaled group.
            m1c      = grp_old/nucH_old(j)
            nucH_new = grp_new/m1c
            rH       = nucH_new/nucH_old(j)
         else
            ! No hydrogen to scale: the group comes back at the reservoir
            ! composition, which is the metal/H the re-seed below writes
            ! (and the bare proton mass where the metals are not scaled).
            if (scale_metals) then
               nucH_new = grp_new/m_1
            else
               nucH_new = grp_new/m_H_amu
            endif
            rH = 0.0d0
         endif
         if (nucHe_old(j) .gt. 1.0d-30) then
            rHe = nucHe_new/nucHe_old(j)
         else
            rHe = 0.0d0
         endif
         if (rH .lt. 0.0d0)  rH  = 0.0d0
         if (rHe .lt. 0.0d0) rHe = 0.0d0
         rBoth = min(rH, rHe)

         do ib = 1, n_bsp
            if (bsp_nHe(ib) .gt. 0) then
               if (bsp_nH(ib) .gt. 0) then
                  f_sp(j,bsp_fsp(ib)) = f_sp(j,bsp_fsp(ib))*rBoth
               else
                  f_sp(j,bsp_fsp(ib)) = f_sp(j,bsp_fsp(ib))*rHe
               endif
            else
               ! Every species of component 1, the oxygen carriers included:
               ! CO carries no hydrogen nucleus, but its mass is component
               ! 1's and it is diluted with the rest of it.
               f_sp(j,bsp_fsp(ib)) = f_sp(j,bsp_fsp(ib))*rH
            endif
         enddo
         ! metals carried with hydrogen at the cell's own metal/H, except
         ! where they have already been written at their own mass fraction
         if (scale_metals) then
            if (nucH_old(j) .gt. 1.0d-30) then
               do im = 1, n_mion
                  f_sp(j,mion_fsp(im)) = f_sp(j,mion_fsp(im))*rH
               enddo
            else if (thereis_metals) then
               ! The cell carries no metal/H of its own, so the metals are
               ! set from the reservoir abundance and the NEW hydrogen count:
               ! they come back when hydrogen does and stay at zero while it
               ! has not.  (With the metals off every column here is zero on
               ! both branches.)
               do ie = 1, n_melem
                  i0e = melem_i0(ie)
                  f_sp(j,mion_fsp(i0e)) = melem_ab(ie)*nucH_new
                  do k = 1, melem_top(ie)
                     f_sp(j,mion_fsp(i0e+k)) = 0.0d0
                  enddo
               enddo
            endif
         endif

         ! deposit the shortfall of each element into its neutral ground stage
         gotH  = 0.0d0
         gotHe = 0.0d0
         do ib = 1, n_bsp
            if (bsp_is_excited_level(ib)) cycle
            if (bsp_nH(ib)  .gt. 0)                                       &
               gotH  = gotH  + dble(bsp_nH(ib)) *f_sp(j,bsp_fsp(ib))
            if (bsp_nHe(ib) .gt. 0)                                       &
               gotHe = gotHe + dble(bsp_nHe(ib))*f_sp(j,bsp_fsp(ib))
         enddo
         if (nucH_new  .gt. gotH)                                         &
            f_sp(j,isp_HI)  = f_sp(j,isp_HI)  + (nucH_new  - gotH)
         if (nucHe_new .gt. gotHe)                                        &
            f_sp(j,isp_HeI) = f_sp(j,isp_HeI) + (nucHe_new - gotHe)
      enddo

      end subroutine project_elements

      ! ------------------------------------------------------------------ !

      subroutine element_transport_residual(rho, Tcode, f_sp, Frho,       &
                                 res_he, sc_he, ok_he,                    &
                                 res_tr, sc_tr, ok_tr, tr_carried,        &
                                 he_terms)
      ! THE ELEMENTAL TRANSPORT BALANCES, MEASURED ON A STATE, stationary
      ! and side-effect free.
      !
      ! What the operator solves is a backward-Euler step of
      !   rho (X^new - X^old)/dt + div(r^2 J)/r^2 + div(F_rho X) = 0
      ! and what a converged wind has to satisfy is the same balance with
      ! the time term gone,
      !   div(r^2 J)/r^2 + div(F_rho X) = 0 ,
      ! transport against transport: helium has no source and no sink, so
      ! there is no reaction term on the right.  It is obtained here the way
      ! the carrier balance obtains its stationary form, by evaluating the
      ! step's own residual at a step long enough that the time term is
      ! numerically absent (dt_stationary below), on the composition the
      ! state carries.  Nothing is re-solved and nothing is projected.
      !
      ! res is [g cm^-3 s^-1] for helium (a mass fraction per unit time
      ! times a density) and [s^-1] for a trace element (a mixing ratio per
      ! unit time).  The scale of a row is the sum of the magnitudes of its
      ! own terms, and the floor under it is a rate 1e-20 of the base
      ! composition carried across the domain in one flow time R0/v0 -- far
      ! below any balance either equation can resolve, and in the row's own
      ! units, so a cell where every term vanishes reads zero rather than
      ! dividing by zero.
      !
      ! Cell 1 carries no equation for either element: it is the Dirichlet
      ! reservoir the operator states, so its residual is reported as zero
      ! against the floor.
      !
      ! THE BOUNDARIES ARE THE CALLER'S, at the top as at the base.  This
      ! routine writes no ghost: the transported quantities below are formed
      ! over 1-Ng:N+Ng from the composition the state carries, and the two
      ! outermost rows are closed with whatever stands in the upper ghosts
      ! of f_sp, the WENO3 face of cell N-1 reaching the first of them.  The
      ! boundary the balance is measured under is therefore the caller's
      ! statement and not this operator's.
      !
      ! The condition at the top that the ELEMENT balance is posed under is
      ! zero gradient -- the diffusive face flux vanishes at face N by the
      ! face coefficients (drift_and_gradient_face_coefficients and
      ! trace_face_coefficients run j = 1, N-1), and with the composition of
      ! the upper ghosts equal to the outermost cell's the advective
      ! reconstruction sees no compositional step across the top, so gas
      ! crossing it carries the column's own composition and no element is
      ! created or destroyed there.  Both production callers write that
      ! condition before they call: the fixed-wind relaxation on its own
      ! working arrays (Xhe in element_diffusion_step, fX in
      ! solve_trace_element_in_hydrogen) and the stationary system on the
      ! composition it assembles (write_species_rows_into_composition of
      ! steady_newton.f90, through project_element_mass_fractions).  A
      ! caller that writes neither is measured under a different boundary.
      !
      ! The lower ghosts are the Dirichlet reservoir the state carries: cell
      ! 1 and below are incoming data and carry no equation here.
      !
      ! NO PERSISTENT MUTATION.  The routines called here compute
      ! coefficients and do not write module state; the five diagnostics of
      ! the last solve that the operator does write are held aside and put
      ! back, and the round trip is asserted rather than assumed.
      !
      ! No cell velocity is taken: the material transport of this balance is
      ! the divergence of the FACE mass fluxes of the state, Frho below.

      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      ! The face mass flux this balance rides on, given by the caller: the
      ! one the mass row of the very state being measured returned.  Every
      ! production caller reads it with face_mass_flux_of_state
      ! (steady_residual.f90) immediately before this call, which refuses a
      ! caller that has not assembled that row rather than measuring the
      ! state against another wind.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Frho
      real*8, dimension(1:N),                 intent(out) :: res_he, sc_he
      ! ONE ROW PER TRACE ELEMENT, not the worst element at each cell.  A
      ! stationary system that carries an element's balance as a ROW needs
      ! that element's own residual; a reduction over elements is a
      ! diagnostic and cannot be a row, and it hid an element whose rates
      ! are small beside another's behind the larger one.  A caller that
      ! wants the worst element still takes it, over the second index.
      real*8, dimension(1:N,n_melem),         intent(out) :: res_tr, sc_tr
      logical,                                intent(out) :: ok_he, ok_tr
      ! Which trace elements the state carries a balance FOR: an element the
      ! run holds none of at the base has no equation, and reporting its
      ! empty row as a satisfied one would say the opposite.
      logical, dimension(n_melem), optional,  intent(out) :: tr_carried
      ! THE HELIUM ROW TERM BY TERM (optional), filled from the same
      ! evaluation that forms res_he and sc_he, so that its residual, its
      ! scale and its measure |res|/max(scale, element_row_scale_floor) are
      ! those numbers and not a second computation of them (the record's
      ! declaration lists its columns).  Asking for it changes nothing the
      ! row is formed from.
      type(element_row_terms_record), optional, intent(out) :: he_terms

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, mass1, Xhe
      real*8, dimension(1-Ng:N+Ng) :: carH, carHe, mcarH, dlnpsi
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, TK, Dco, Gco, dt_stat
      real*8, dimension(1-Ng:N+Ng) :: rp, rep, nH_phys, ntot_phys, dmeff
      real*8, dimension(1-Ng:N+Ng) :: nX, nXold, DcoX, GcoX, zb1
      real*8, dimension(1-Ng:N+Ng) :: Dneut, eEf, ne_phys, zbHe, ne_rel
      real*8, dimension(1-Ng:N+Ng) :: mres, resc, dscc, fX, fXold
      ! The face terms of the helium row and its face composition, and the
      ! face mass flux in g cm^-2 s^-1, for the record he_terms.
      real*8, dimension(1-Ng:N+Ng,4) :: row_faces
      real*8, dimension(1-Ng:N+Ng)   :: Yface_he, Fface_cgs
      ! The two factors that put the divergence of a face element mass flux
      ! into the units of each row.
      real*8, dimension(1-Ng:N+Ng) :: cadvf, wYtr, cadvX
      real*8, dimension(1-Ng:N+Ng,n_hcar) :: yH
      integer, dimension(1-Ng:N+Ng) :: idom
      real*8, dimension(0:N) :: Agrd, Bdrf, Jf, dJl, dJr, PL, PR
      integer, dimension(0:N) :: updrf
      ! The five diagnostics of the last solve, held aside (see the header).
      real*8  :: sv_over, sv_under, sv_resid, sv_trace
      integer :: sv_steps
      ! The face counters of species_face_fraction count the run's own
      ! faces; the reconstruction this measurement makes must not add to
      ! them (the same rule as write_element_flux_profile).
      real*8  :: sv_exc
      integer :: sv_bnd
      real*8  :: tscale, X_base, m_1, fXbase, rnorm
      real*8  :: floor_he, floor_tr
      integer :: j, im

      res_he = 0.0d0;  sc_he = 1.0d0;  ok_he = .false.
      res_tr = 0.0d0;  sc_tr = 1.0d0;  ok_tr = .false.
      if (present(tr_carried)) tr_carried = .false.
      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      sv_over  = he_fraction_over_one
      sv_under = he_fraction_under_zero
      sv_steps = he_fraction_newton_steps
      sv_resid = he_fraction_newton_resid
      sv_trace = trace_ratio_under_zero
      sv_bnd   = n_species_faces_bounded
      sv_exc   = species_face_excursion

      m_1 = mass_per_H_nucleus_without_He()
      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      Xhe  = m_He_amu*nucHe/msum

      TK       = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      rp       = r*R0
      rep      = r_edg*R0
      rho_phys = rho*n0*mu*msum
      call carrier_counts(f_sp, carH, carHe, mcarH)
      ntot_phys = (carH + carHe)*rho*n0
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dco, Dneut, idom)
      tscale  = R0/v0
      ! The stationary limit: a step long against every transport time of
      ! the column, so rho (X - X)/dt is zero and the scale of the row is
      ! the transport terms alone.  Same device as the carrier balance.
      dt_stat = 1.0d30
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1, eEf,   &
                                dlnpsi)
      X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)
      call drift_and_gradient_face_coefficients(rho_phys, Dco, Gco, rp,   &
                                                Agrd, Bdrf, updrf)
      ! THE ROW IS THE OPERATOR'S OWN STEP AT AN INFINITE STEP LENGTH.  Both
      ! halves come from composition_residual: the diffusive divergence, and
      ! the material advection as the divergence of F_rho X^face on the same
      ! faces, areas and volume the mass row of this cell uses.  The
      ! fixed-wind relaxation solves that same expression to zero, so the
      ! relaxation's fixed point IS this row's zero and the elemental Picard
      ! alternation is one operator and not two.
      ! Code density times code velocity over code length is a rate in units
      ! of v0/R0, and n0 mu msum is the mass one unit of the code density
      ! carries, so the two factors put the divergence in the row's own
      ! [g cm^-3 s^-1].
      cadvf = n0*mu*msum*v0/R0
      if (present(he_terms)) then
         call composition_residual(Xhe, Xhe, rho_phys, dt_stat,           &
                                   Frho, cadvf, Agrd, Bdrf, updrf, 2,     &
                                   .true., mres, Jf, dJl, dJr, rnorm,     &
                                   resc, dscc, faces_out = row_faces,     &
                                   Yface_out = Yface_he)
      else
         call composition_residual(Xhe, Xhe, rho_phys, dt_stat,           &
                                   Frho, cadvf, Agrd, Bdrf, updrf, 2,     &
                                   .true., mres, Jf, dJl, dJr, rnorm,     &
                                   resc, dscc)
      endif
      do j = 2, N
         floor_he  = 1.0d-20*rho_phys(j)*X_base/tscale
         res_he(j) = resc(j)
         sc_he(j)  = max(dscc(j), floor_he)
      enddo
      sc_he(1) = max(1.0d-20*rho_phys(1)*X_base/tscale, 1.0d-300)
      ok_he    = .true.

      ! The helium row term by term, from the evaluation just made.  Cell 1
      ! is the Dirichlet reservoir and carries no equation, so its terms are
      ! left at zero and its scale is the floor, as res_he and sc_he state.
      if (present(he_terms)) then
         allocate(he_terms%rho_phys(1:N), he_terms%X(1:N),                &
                  he_terms%D12(1:N), he_terms%Kzz(1:N),                   &
                  he_terms%dif_in(1:N), he_terms%dif_out(1:N),            &
                  he_terms%adv_in(1:N), he_terms%adv_out(1:N),            &
                  he_terms%res(1:N), he_terms%scale(1:N),                 &
                  he_terms%floor(1:N), he_terms%F_rho(0:N),               &
                  he_terms%Y_face(0:N), he_terms%J(0:N))
         call face_mass_flux_cgs(Frho, msum, Fface_cgs)
         do j = 1, N
            he_terms%rho_phys(j) = rho_phys(j)
            he_terms%X(j)        = Xhe(j)
            he_terms%D12(j)      = Dco(j)
            he_terms%Kzz(j)      = kzz_cell(j)
            he_terms%dif_in(j)   = row_faces(j,1)
            he_terms%dif_out(j)  = row_faces(j,2)
            he_terms%adv_in(j)   = row_faces(j,3)
            he_terms%adv_out(j)  = row_faces(j,4)
            he_terms%res(j)      = res_he(j)
            he_terms%scale(j)    = sc_he(j)
            he_terms%floor(j)    = 1.0d-20*rho_phys(j)*X_base/tscale
         enddo
         he_terms%floor(1) = sc_he(1)
         do j = 0, N
            he_terms%F_rho(j)  = Fface_cgs(j)
            he_terms%Y_face(j) = Yface_he(j)
            he_terms%J(j)      = Jf(j)
         enddo
      endif

      ! --- the trace elements, each against the hydrogen background ---
      if (.not. (he_metal_diffusion .and. thereis_metals)) then
         call reinstate_diffusion_diagnostics(sv_over, sv_under,          &
                                     sv_steps, sv_resid, sv_trace)
         n_species_faces_bounded = sv_bnd
         species_face_excursion  = sv_exc
         return
      endif
      call element_nucleus_counts(f_sp, nucH, nucHe)
      nH_phys = nucH*rho*n0
      where (nH_phys .lt. 1.0d-30) nH_phys = 1.0d-30
      call carrier_counts(f_sp, carH, carHe, mcarH)
      ntot_phys = (carH + carHe)*rho*n0
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
      ne_phys = max(ne_rel*rho*n0, 1.0d0)
      call carrier_fractions(f_sp, n_hcar, hcar_isp, yH)
      sc_tr = 1.0d-300
      do im = 1, n_melem
         call trace_element_transport_coefficients(im, f_sp, rho, TK,     &
                  ntot_phys, ne_phys, yH, mcarH, zb1, eEf, dlnpsi,        &
                  nX, nXold, DcoX, GcoX)
         fX     = nX/max(nH_phys, 1.0d-30)
         fXbase = nXold(1)/nH_phys(1)
         if (fXbase .le. 0.0d0) cycle
         if (present(tr_carried)) tr_carried(im) = .true.
         floor_tr = 1.0d-20*fXbase/tscale
         ! The element's MASS fraction per unit of the mixing ratio, and the
         ! factor that returns the code divergence of that mass flux to the
         ! mixing ratio per second the row is written in: n0 msum v0/R0 makes
         ! it a mass rate, 1/A_X counts nuclei, 1/n_H takes the ratio.
         wYtr  = melem_A(im)*nucH/msum
         cadvX = n0*msum*v0/R0/melem_A(im)/max(nH_phys, 1.0d-30)
         call trace_face_coefficients(nH_phys, DcoX, GcoX, rp, PL, PR)
         ! The same step at an infinite step length, so the time term is
         ! absent and the row is the transport terms alone -- and it is the
         ! same expression the fixed-wind relaxation of this element solves
         ! to zero.
         fXold = fX
         if (trace_row_terms_diag)                                        &
            write(*,'(A,A)') ' (element row terms) element ',             &
                 trim(melem_name(im))
         call trace_composition_residual(fXold, fX, nH_phys, dt_stat,      &
                                         PL, PR, Frho, wYtr, cadvX,        &
                                         .true., mres, rnorm, resc, dscc)
         do j = 2, N
            res_tr(j,im) = resc(j)
            sc_tr(j,im)  = max(dscc(j), floor_tr)
         enddo
      enddo
      ok_tr = .true.

      call reinstate_diffusion_diagnostics(sv_over, sv_under, sv_steps,   &
                                           sv_resid, sv_trace)
      n_species_faces_bounded = sv_bnd
      species_face_excursion  = sv_exc

      end subroutine element_transport_residual

      ! ------------------------------------------------------------------ !

      subroutine element_nucleus_face_flux(rho, Tcode, f_sp, Frho,        &
                                           Fadv, Jdif, n_el, m_one, n_one)
      ! THE HELIUM NUCLEUS FLUX THROUGH EVERY FACE OF THE COLUMN, the one
      ! object the element row is the divergence of.  For the state
      ! (rho, Tcode, f_sp) and the face mass flux Frho of that state's own
      ! mass row it returns, at the faces f = 0 ... N:
      !
      !   Fadv(f) = F_rho(f) Y_He(f)   the advective face element mass flux,
      !             IN THE UNITS OF THE Frho THE CALLER PASSES.  Y_He(f) is
      !             the reconstructed and upwinded face helium mass fraction
      !             of species_face_fraction, the same one the Runge-Kutta
      !             stages put on the face, and it is dimensionless, so this
      !             routine neither scales nor rescales the flux it is given.
      !             With the code face mass flux of the mass row (rho v in
      !             units of n0 mu and v0, the mixture mass msum(j) still
      !             outside it) Fadv is in those code units and
      !             n0 mu msum v0 Fadv is g cm^-2 s^-1; with a face mass flux
      !             already in g cm^-2 s^-1, which the element flux profile
      !             writer below passes, Fadv is in g cm^-2 s^-1 directly.
      !             Jdif below is in g cm^-2 s^-1 either way, so a caller
      !             that adds the two halves converts the advective one
      !             first.
      !   Jdif(f) = -A_grd (X_r - X_l) - B_drf [X(1-X)](f)   [g cm^-2 s^-1],
      !             the diffusive (gradient, eddy and settling drift) face
      !             flux of element_face_flux.  Faces 0 and N carry none:
      !             their gradient and drift coefficients are zero
      !             (drift_and_gradient_face_coefficients), the base being
      !             the Dirichlet reservoir and the outer face an outflow,
      !             so an element crosses either end only with the gas.
      !             The face is evaluated all the same, so the zero is the
      !             operator's own and is not asserted here.
      !   n_el(f) = the helium NUCLEUS density at the face [cm^-3], the
      !             arithmetic mean of its two cells, which is the face rule
      !             every other coefficient of this operator uses.
      !   m_one(f)= the mass of component 1 per hydrogen nucleus at the face
      !             [g], the same mean.  In the atomic region it is the
      !             reservoir m_1 = mass_per_H_nucleus_without_He(); where
      !             the metals or the molecules have moved it is the cell's
      !             own.
      !   n_one(f)= the HYDROGEN nucleus density at the face [cm^-3],
      !             optional and the same mean, the quantity n_el is for
      !             helium.  A row written per hydrogen nucleus -- the
      !             ionization stages of hydrogen are -- needs it and it is
      !             the same face rule, so it is returned here rather than
      !             rebuilt by the caller.
      !
      ! THE ROW IS WRITTEN FROM THESE ARRAYS.  With A(f) and V(j) the face
      ! area and shell volume of spherical_face_area_and_cell_volume, the
      ! helium row of cell j is, term for term,
      !
      !   res(j) = rho(j) [X(j) - X^old(j)]/dt
      !          + [ A_+ R0^2 Jdif(j) - A_- R0^2 Jdif(j-1) ]/(V_j R0^3)
      !          + cadvf(j) [ A_+ Fadv(j) - A_- Fadv(j-1) ]/V_j ,
      !          cadvf(j) = n0 mu msum(j) v0/R0 ,
      !
      ! which is what composition_residual evaluates; the acceptance row
      ! element_row_is_the_divergence_of_the_exposed_flux holds the two to
      ! the last bit.  There is therefore ONE spelling of the element flux,
      ! and a caller that needs the flux itself (the ionization-stage rows
      ! of L12 stage B, the boundary budget of a diagnostic) reads it here
      ! instead of rebuilding it.
      !
      ! The unit factor cadvf carries the cell's own mixture mass msum(j),
      ! so the column sum telescopes to the two boundary face fluxes up to
      ! the variation of msum between neighbours, which is the state's mass
      ! closure (1e-15 on a closed state) and not a property of the grid.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Frho
      real*8, dimension(0:N),                 intent(out) :: Fadv, Jdif
      real*8, dimension(0:N),                 intent(out) :: n_el, m_one
      real*8, dimension(0:N),       optional, intent(out) :: n_one

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, mass1, Xhe
      real*8, dimension(1-Ng:N+Ng) :: Yf, Fs, nHe_cell, m1_cell, nH_cell
      real*8  :: sv_over, sv_under, sv_resid, sv_trace
      real*8  :: sv_exc
      integer :: sv_steps, sv_bnd
      integer :: j

      Fadv  = 0.0d0
      Jdif  = 0.0d0
      n_el  = 0.0d0
      m_one = 0.0d0
      if (present(n_one)) n_one = 0.0d0
      if (.not. thereis_He) return

      ! THE FACE COUNTERS OF THE RUN COUNT THE RUN'S OWN FACES.  The
      ! advective half below reconstructs a face composition, and
      ! species_face_fraction counts every face state it has to scale back
      ! onto [0,1] and keeps the largest excursion it saw.  Those are
      ! statements about the trajectory the run took, so an evaluation of
      ! the flux of a state that is already in hand -- a diagnostic, a
      ! stage row, an acceptance row -- must not add to them.  The restore
      ! belongs here rather than at each call site, because every caller of
      ! this routine has the same obligation.
      sv_bnd   = n_species_faces_bounded
      sv_exc   = species_face_excursion

      sv_over  = he_fraction_over_one
      sv_under = he_fraction_under_zero
      sv_steps = he_fraction_newton_steps
      sv_resid = he_fraction_newton_resid
      sv_trace = trace_ratio_under_zero

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      Xhe = m_He_amu*nucHe/msum

      ! The diffusive half, face by face, in the branch the Peclet switch
      ! selected -- the same call composition_residual makes.
      call helium_diffusive_face_flux(rho, Tcode, f_sp, Jdif)

      ! The advective half: the face mass flux of the mass row carrying the
      ! reconstructed face mass fraction, the same two routines the stages
      ! and element_advective_divergence call.
      call species_face_fraction(Xhe, Frho, Yf)
      call species_face_flux(Frho, Yf, Fs)
      do j = 0, N
         Fadv(j) = Fs(j)
      enddo

      ! The two face state quantities the stage rows are written against.
      nHe_cell = nucHe*rho*n0
      nH_cell  = nucH *rho*n0
      m1_cell  = mass1*mu/max(nucH, 1.0d-300)
      do j = 0, N
         n_el(j)  = 0.5d0*(nHe_cell(j) + nHe_cell(j+1))
         m_one(j) = 0.5d0*(m1_cell(j)  + m1_cell(j+1))
      enddo
      if (present(n_one)) then
         do j = 0, N
            n_one(j) = 0.5d0*(nH_cell(j) + nH_cell(j+1))
         enddo
      endif

      n_species_faces_bounded = sv_bnd
      species_face_excursion  = sv_exc

      call reinstate_diffusion_diagnostics(sv_over, sv_under, sv_steps,   &
                                           sv_resid, sv_trace)

      end subroutine element_nucleus_face_flux

      ! ------------------------------------------------------------------ !

      subroutine helium_diffusive_face_flux(rho, Tcode, f_sp, Jdif)
      ! The diffusive helium mass flux of the state (rho, Tcode, f_sp) at
      ! the faces f = 0 ... N [g cm^-2 s^-1], positive outward:
      !
      !    Jdif(f) = -A_grd (X_r - X_l) - B_drf [X(1-X)](f) ,
      !
      ! the gradient, eddy and settling-drift flux of element_face_flux in
      ! the branch the Peclet switch selects, with the coefficients of
      ! drift_and_gradient_face_coefficients.  Faces 0 and N carry none.
      ! The one spelling of it for a state in hand: element_nucleus_face_flux
      ! returns it beside the advective half, and the interdiffusion enthalpy
      ! flux of the energy equation is formed from it.  Nothing here writes
      ! module state.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(0:N),                 intent(out) :: Jdif

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, mass1, Xhe
      real*8, dimension(1-Ng:N+Ng) :: rho_phys, Dco, Gco, rp
      real*8, dimension(1-Ng:N+Ng) :: dmeff, zb1, eEf, dlnpsi
      real*8, dimension(0:N)       :: Agrd, Bdrf
      integer, dimension(0:N)      :: updrf
      real*8  :: dJl, dJr
      integer :: j

      Jdif = 0.0d0
      if (.not. thereis_He) return

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      Xhe = m_He_amu*nucHe/msum
      rp       = r*R0
      rho_phys = rho*n0*mu*msum
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dco)
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1, eEf,   &
                                dlnpsi)
      call drift_and_gradient_face_coefficients(rho_phys, Dco, Gco, rp,   &
                                                Agrd, Bdrf, updrf)
      do j = 0, N
         call element_face_flux(Xhe(j), Xhe(j+1), Agrd(j), Bdrf(j),       &
                                updrf(j), Jdif(j), dJl, dJr)
      enddo

      end subroutine helium_diffusive_face_flux

      ! ------------------------------------------------------------------ !

      subroutine trace_element_diffusive_face_flux(rho, Tcode, f_sp, JX)
      ! The diffusive MASS flux of every trace metal element relative to the
      ! hydrogen it diffuses through, at the faces f = 0 ... N
      ! [g cm^-2 s^-1], positive outward, for the state in hand:
      !
      !    JX(f,im) = A_X m_H [ PL(f) fX(f) + PR(f) fX(f+1) ] ,
      !
      ! fX = n_X/n_H and PL, PR of trace_face_coefficients with the
      ! coefficients of trace_element_transport_coefficients, i.e. the flux
      ! the stationary trace row of element_transport_residual is the
      ! divergence of.  Zero unless He_metal_diffusion moves the metals.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(0:N,n_melem),         intent(out) :: JX

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, nH_phys, ntot_phys
      real*8, dimension(1-Ng:N+Ng) :: carH, carHe, mcarH, dlnpsi, TK, rp
      real*8, dimension(1-Ng:N+Ng) :: Gco, dmeff, zb1, eEf, zbHe, ne_rel
      real*8, dimension(1-Ng:N+Ng) :: ne_phys, nX, nXold, DcoX, GcoX, fX
      real*8, dimension(1-Ng:N+Ng,n_hcar) :: yH
      real*8, dimension(0:N) :: PL, PR
      integer :: im, j

      JX = 0.0d0
      if (.not. (he_metal_diffusion .and. thereis_metals)) return

      TK = Tcode*T0
      where (TK .lt. 1.0d0) TK = 1.0d0
      rp = r*R0
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1, eEf,   &
                                dlnpsi)
      call element_nucleus_counts(f_sp, nucH, nucHe)
      nH_phys = nucH*rho*n0
      where (nH_phys .lt. 1.0d-30) nH_phys = 1.0d-30
      call carrier_counts(f_sp, carH, carHe, mcarH)
      ntot_phys = (carH + carHe)*rho*n0
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      call mean_charges_and_electrons(f_sp, zb1, zbHe, ne_rel)
      ne_phys = max(ne_rel*rho*n0, 1.0d0)
      call carrier_fractions(f_sp, n_hcar, hcar_isp, yH)
      do im = 1, n_melem
         call trace_element_transport_coefficients(im, f_sp, rho, TK,     &
                  ntot_phys, ne_phys, yH, mcarH, zb1, eEf, dlnpsi,        &
                  nX, nXold, DcoX, GcoX)
         fX = nX/max(nH_phys, 1.0d-30)
         call trace_face_coefficients(nH_phys, DcoX, GcoX, rp, PL, PR)
         do j = 0, N
            JX(j,im) = melem_A(im)*mu*(PL(j)*fX(j) + PR(j)*fX(j+1))
         enddo
      enddo

      end subroutine trace_element_diffusive_face_flux

      ! ------------------------------------------------------------------ !

      logical function interdiffusion_enthalpy_active()
      ! Whether the energy equation carries the interdiffusion enthalpy
      ! flux: only where the elements move relative to one another, i.e.
      ! with He_diffusion on and helium in the run, and unless the key
      ! "Interdiffusion enthalpy flux: False" removes it to match a
      ! published model that omits it (parameters.f90).
      interdiffusion_enthalpy_active = he_diffusion .and. thereis_He      &
                                       .and. interdiffusion_enthalpy_flux
      end function interdiffusion_enthalpy_active

      ! ------------------------------------------------------------------ !

      subroutine component_specific_enthalpies(Tcode, f_sp, hHe, h1,      &
                                               hHgrp, hX)
      ! THE SENSIBLE SPECIFIC ENTHALPIES [erg/g] OF THE COMPONENTS THE
      ! ELEMENT OPERATOR MOVES, cell by cell (ghosts included):
      !
      !   hHe       helium, every stage, with the free electrons its ions
      !             gave up;
      !   h1        component 1 as the binary operator defines it
      !             (mixture_mass_split): the hydrogen carriers, the oxygen
      !             and carbon carriers, and the metals where their mass is
      !             in the mixture (eos_include_metals), with their
      !             electrons;
      !   hHgrp     component 1 without the metals: the group that recoils
      !             against a metal diffusing through hydrogen;
      !   hX(:,im)  metal element im, every stage, with its electrons.
      !
      ! Per particle, h_s = e_s + p_s/n_s = e_s + kT (Cook 2009, eq. 13),
      ! with e_s the internal energy of the code's caloric equation of state
      ! (caloric_eos): kT/(gamma_ad - 1) for every atom, ion, electron and
      ! for H2+, H3+, HeH+ and the oxygen carriers, plus the rovibrational
      ! energy k u_rv(T) of H2 (the same ladder, same table, and zero under
      ! "Caloric EOS: monatomic").  So h_s = kT gamma_ad/(gamma_ad - 1) for
      ! every particle, H2 adding k u_rv.  Each ion's electrons are counted
      ! with it: the ambipolar field makes the electrons follow the ions with
      ! no current, n_e w_e = sum_s Z_s n_s w_s, so an ion of charge Z moves
      ! Z electrons with its own element.  HeH+ carries a nucleus of each
      ! element; its particle and electron enthalpy is shared between the two
      ! components in the proportion of its nuclei, as its mass is
      ! (mixture_mass_split puts the helium nucleus in the helium component
      ! and the proton in component 1).  The He 2^3S column is a level of
      ! He I and is skipped, as the particle count of the equation of state
      ! skips it (calc_ntot).
      !
      ! FORMATION (IONIZATION, DISSOCIATION) ENERGY IS NOT IN THESE
      ! ENTHALPIES, and must not be: the conserved energy u(3) of this code
      ! is kinetic plus sensible, and the chemical energy of the species is
      ! booked separately through the formation-energy density
      ! (formation_energy_density), whose reference in the coupled source
      ! step is the composition the element transport has just written
      ! (EXHALE_main, u_form_old_csm after element_diffusion_step), so the
      ! formation energy of the nuclei moved rides with them and costs the
      ! thermal energy nothing.  The stationary energy row holds no
      ! formation-energy density at all: chemical energy enters it only as
      ! the reaction heating and cooling rates inside heat - cool, released
      ! where the species react, so nuclei moved by diffusion carry their
      ! formation energy there exactly as advected ones do.  What the
      ! transport of particles carries and no other term books is their
      ! sensible enthalpy.
      !
      ! A component a cell holds none of is given the enthalpy of its
      ! neutral ground species, which is what enters such a cell first; the
      ! flux carrying it is then the gradient flux into an empty cell.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: hHe, h1, hHgrp
      real*8, dimension(1-Ng:N+Ng,n_melem),   intent(out) :: hX

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, mass1, msum, mmetal
      real*8, dimension(n_melem)   :: enth_el, cnt_el
      real*8  :: cp_part, kT_erg, TKc, urv, crv, enth_He, enth_1, enth_met
      real*8  :: part_enth, w_He, dens
      integer :: j, ib, im, iel

      ! gamma_ad/(gamma_ad - 1): the enthalpy per particle in units of kT
      ! of the caloric EOS for every species but H2 (5/2 for 5/3).
      cp_part = gamma_ad/(gamma_ad - 1.0d0)

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum, mmetal)

      do j = 1-Ng, N+Ng
         TKc = Tcode(j)*T0
         if (TKc .lt. 1.0d0) TKc = 1.0d0
         kT_erg = kb_erg*TKc
         urv = 0.0d0
         if (thereis_mol) call h2_rovibrational_energy_and_heat_capacity( &
                                  TKc, urv, crv)
         enth_He  = 0.0d0
         enth_1   = 0.0d0
         enth_met = 0.0d0
         do ib = 1, n_bsp
            if (bsp_is_excited_level(ib)) cycle
            dens = f_sp(j,bsp_fsp(ib))
            if (dens .eq. 0.0d0) cycle
            part_enth = cp_part*kT_erg*(1.0d0 + dble(bsp_charge(ib)))
            if (bsp_fsp(ib) .eq. isp_H2) part_enth = part_enth + kb_erg*urv
            w_He = 0.0d0
            if (bsp_nHe(ib) .gt. 0) w_He = dble(bsp_nHe(ib))              &
                                   /dble(bsp_nHe(ib) + bsp_nH(ib))
            enth_He = enth_He + w_He*dens*part_enth
            enth_1  = enth_1  + (1.0d0 - w_He)*dens*part_enth
         enddo
         enth_el = 0.0d0
         cnt_el  = 0.0d0
         if (eos_include_metals .and. thereis_metals) then
            do im = 1, n_mion
               dens = f_sp(j,mion_fsp(im))
               iel  = mion_elem(im)
               part_enth = cp_part*kT_erg*(1.0d0 + dble(mion_stage(im)))
               enth_el(iel) = enth_el(iel) + dens*part_enth
               cnt_el(iel)  = cnt_el(iel)  + dens
               enth_met     = enth_met + dens*part_enth
            enddo
            enth_1 = enth_1 + enth_met
         endif
         ! Per gram: the masses are those of mixture_mass_split, in units
         ! of the hydrogen atom per unit of the code density, as the
         ! particle counts above are.
         if (m_He_amu*nucHe(j) .gt. 1.0d-30) then
            hHe(j) = enth_He/(m_He_amu*nucHe(j)*mu)
         else
            hHe(j) = cp_part*kT_erg/(m_He_amu*mu)
         endif
         if (mass1(j) .gt. 1.0d-30) then
            h1(j) = enth_1/(mass1(j)*mu)
         else
            h1(j) = cp_part*kT_erg/(m_H_amu*mu)
         endif
         if (mass1(j) - mmetal(j) .gt. 1.0d-30) then
            hHgrp(j) = (enth_1 - enth_met)/((mass1(j) - mmetal(j))*mu)
         else
            hHgrp(j) = cp_part*kT_erg/(m_H_amu*mu)
         endif
         do iel = 1, n_melem
            if (melem_A(iel)*cnt_el(iel) .gt. 1.0d-30) then
               hX(j,iel) = enth_el(iel)/(melem_A(iel)*cnt_el(iel)*mu)
            else
               hX(j,iel) = cp_part*kT_erg/(melem_A(iel)*mu)
            endif
         enddo
      enddo

      end subroutine component_specific_enthalpies

      ! ------------------------------------------------------------------ !

      subroutine interdiffusion_enthalpy_face_flux(Tcode, f_sp, Jhe, qd,  &
                                                   JX)
      ! THE INTERDIFFUSION ENTHALPY FLUX at the faces f = 0 ... N
      ! [erg cm^-2 s^-1], positive outward.
      !
      ! In a mixture whose species move with velocities u + w_s relative to
      ! the mass-weighted mean velocity u of the hydrodynamics, the energy
      ! equation is
      !
      !    dE/dt + div[(E + p) u] = div(tau.u - q_c - q_d) ,
      !    q_d = sum_s h_s J_s ,  J_s = rho_s w_s ,  sum_s J_s = 0 ,
      !
      ! (Cook 2009, Phys. Fluids 21, 055109, eqs. 11-13; the Dufour flux and
      ! the kinetic energy of the relative motion are left out there and
      ! here).  The same term is the difference between the heat flow of a
      ! species measured in the mean frame and in its own drift frame,
      ! q_s* = q_s + (5/2) p_s w_s + ... (Schunk 1977, Rev. Geophys. Space
      ! Phys. 15, 429, eqs. 15b and 16): summing the species energy
      ! equations (Schunk eq. 20c, written about u_s = u + w_s) over s gives
      ! the mixture equation with div sum_s (5/2) p_s w_s = div sum_s h_s J_s
      ! for particles that store (3/2) kT.  Because h_s multiplies a flux
      ! relative to the MASS-weighted velocity, it is the enthalpy and not
      ! the internal energy (Cook 2009, text after eq. 13).
      !
      ! With the binary element model of this module every species of an
      ! element moves with its element, the electrons with their ions, and
      ! J_1 = -J_He, so
      !
      !    q_d = J_He [ h_He - h_1 ]  +  sum_X J_X [ h_X - h_H ] ,
      !
      ! h_c the sensible specific enthalpies of component_specific_
      ! enthalpies.  The metal sum is present only where He_metal_diffusion
      ! moves each metal against hydrogen (J_X of trace_element_diffusive_
      ! face_flux or of the marching step); the hydrogen group then recoils
      ! against it (project_elements), which is why the metal term is taken
      ! against h_H, component 1 without its metals.  Without that option the
      ! metals are part of component 1 and ride in h_1.  The metal term is
      ! also dropped where the metals carry no mass in the mixture
      ! (eos_include_metals off): their transport then moves neither mass
      ! nor energy of the gas the hydrodynamics describes.
      !
      ! J is the operator's whole diffusive flux, gradient, eddy and
      ! settling drift together.  The eddy part is included: the same eddy
      ! mixing that moves the elements carries their enthalpy (Cook 2009,
      ! section II.B: "if J_i is nonzero, then q_d is potentially important,
      ! regardless of whether J_i represents a molecular ..., subgrid-scale
      ! ..., or turbulent diffusion ... flux").
      !
      ! Face values: h_c at a face is the arithmetic mean of its two cells,
      ! the face rule of every coefficient of this operator.  The flux is
      ! zero wherever J is, in particular at faces 0 and N, whose diffusive
      ! coefficients the operator sets to zero.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(0:N),                 intent(in)  :: Jhe
      real*8, dimension(0:N),                 intent(out) :: qd
      real*8, dimension(0:N,n_melem), optional, intent(in) :: JX

      real*8, dimension(1-Ng:N+Ng) :: hHe, h1, hHgrp
      real*8, dimension(1-Ng:N+Ng,n_melem) :: hX
      integer :: j, im

      call component_specific_enthalpies(Tcode, f_sp, hHe, h1, hHgrp, hX)
      do j = 0, N
         qd(j) = Jhe(j)*0.5d0*((hHe(j) - h1(j)) + (hHe(j+1) - h1(j+1)))
      enddo
      if (present(JX) .and. he_metal_diffusion .and. thereis_metals      &
          .and. eos_include_metals) then
         do im = 1, n_melem
            do j = 0, N
               qd(j) = qd(j) + JX(j,im)*0.5d0*                            &
                       ((hX(j,im) - hHgrp(j)) + (hX(j+1,im) - hHgrp(j+1)))
            enddo
         enddo
      endif

      end subroutine interdiffusion_enthalpy_face_flux

      ! ------------------------------------------------------------------ !

      subroutine interdiffusion_enthalpy_divergence(qd, divq)
      ! The divergence of the face flux qd [erg cm^-2 s^-1] in every cell
      ! the element operator carries an equation for, in the code units of
      ! the energy row (q0 = n0 mu v0^3/R0):
      !
      !    divq(j) = [ A_+ qd(j) - A_- qd(j-1) ] / (V_j R0 q0) ,  j = 2 ... N,
      !
      ! with A and V of spherical_face_area_and_cell_volume, the one
      ! geometry of the element rows, so an interior face cancels between
      ! the two cells that share it and sum_j V_j divq(j) R0 q0 is
      ! A_N qd(N) - A_1 qd(1): the energy the term moves is conserved to
      ! round-off within the column.  Cell 1 is the Dirichlet reservoir of
      ! the element operator and carries no element equation (its
      ! composition is prescribed), so it carries no interdiffusion energy
      ! either: the enthalpy that crosses face 1 comes with helium that the
      ! reservoir supplies, exactly as the element flux through face 1 does.
      ! Cell 1 and the ghosts are returned as zero.
      !
      ! THE BUDGET OF THE WHOLE HYDRODYNAMIC COLUMN (cells 1 ... N, lower
      ! face 0). Cell 1's composition is steady, so the element flux the
      ! reservoir puts through face 0 is the flux through face 1, A_0 J_0 =
      ! A_1 J_1, and with it the enthalpy, A_0 q_0 = A_1 qd(1): cell 1 gains
      ! and loses the same energy, divq(1) = 0 is that balance and not an
      ! omission, and
      !    sum_{j=1..N} V_j divq(j) R0 q0 = A_N qd(N) - A_1 qd(1)
      ! is the outflow at face N minus the RESERVOIR'S ENERGY SUPPLY at face
      ! 0, which is A_1 qd(1) and not the qd(0) of the face-flux array (the
      ! element operator carries no flux at face 0). The conservation budget
      ! export states the same (conservation_budget, Sidf).
      real*8, dimension(0:N),       intent(in)  :: qd
      real*8, dimension(1-Ng:N+Ng), intent(out) :: divq

      real*8, dimension(0:N) :: fa
      real*8, dimension(1:N) :: cv
      integer :: j
      divq = 0.0d0
      call spherical_face_area_and_cell_volume(fa, cv)
      do j = 2, N
         divq(j) = interdiffusion_enthalpy_scale()*                     &
                   (fa(j)*qd(j) - fa(j-1)*qd(j-1))/(cv(j)*R0*q0)
      enddo

      end subroutine interdiffusion_enthalpy_divergence

      ! ------------------------------------------------------------------ !

      real*8 function interdiffusion_enthalpy_scale()
      ! A CONTINUATION FACTOR, NOT PHYSICS. Measurement key
      ! EXHALE_INTERDIFF_ENTH_SCALE = s, 0 <= s <= 1 (1 when unset, which
      ! is the equation): interdiffusion_enthalpy_divergence multiplies the
      ! divergence by s. It exists to reach a state of the full equation
      ! (s = 1) from one solved without the term (s = 0) in steps, where the
      ! term moves the solution by O(1) and a Newton solve from the s = 0
      ! state does not converge (the atomic He/H 9.7 cases at reduced XUV,
      ! md/Update_EXHALE_stage3.md section 38). A state solved at 0 < s < 1
      ! is a step of that continuation and not a prediction of the model.
      ! The factor changes the energy equation, so it is part of the
      ! identity of every state written under it: the 'interdiff_enth'
      ! token of the options field carries its value wherever s /= 1
      ! (load_IC, opt_value), and a restart across a change of s must name
      ! that token on its "Restart option change:" line.
      real*8,  save :: enthalpy_scale = 1.0d0
      logical, save :: enthalpy_scale_read = .false.
      character(len=32) :: env
      integer :: st

      if (.not. enthalpy_scale_read) then
         enthalpy_scale_read = .true.
         call get_environment_variable('EXHALE_INTERDIFF_ENTH_SCALE', env, &
                                       status=st)
         if (st .eq. 0 .and. len_trim(env) .gt. 0) then
            read(env,*,iostat=st) enthalpy_scale
            if (st .ne. 0 .or. .not. (enthalpy_scale .ge. 0.0d0 .and.      &
                enthalpy_scale .le. 1.0d0)) then
               write(*,*) '(diffusion) ERROR: EXHALE_INTERDIFF_ENTH_SCALE'//  &
                          ' takes a number from 0 to 1, not "'//            &
                          trim(env)//'".'
               error stop 1
            endif
            write(*,'(A,F8.5,A)') ' (diffusion) interdiffusion enthalpy'//   &
                 ' flux scaled by', enthalpy_scale, ' (continuation factor'//&
                 ' EXHALE_INTERDIFF_ENTH_SCALE; 1 is the equation)'
         endif
      endif
      interdiffusion_enthalpy_scale = enthalpy_scale
      end function interdiffusion_enthalpy_scale

      ! ------------------------------------------------------------------ !

      subroutine interdiffusion_enthalpy_divergence_of_state(rho, Tcode,  &
                                                             f_sp, divq,  &
                                                             qd_out)
      ! The divergence of the interdiffusion enthalpy flux of the state
      ! (rho, Tcode, f_sp) in the code units of the energy row, the term the
      ! stationary energy row adds (steady_residual, assemble_residual):
      ! the element fluxes of the state itself (helium_diffusive_face_flux,
      ! trace_element_diffusive_face_flux), the face flux of
      ! interdiffusion_enthalpy_face_flux and the divergence of
      ! interdiffusion_enthalpy_divergence -- the same three steps the
      ! marching update takes with the fluxes of its own element step.
      ! Zero, with no arithmetic, unless interdiffusion_enthalpy_active().
      ! Writes no module state.  qd_out (optional) returns the face flux.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: divq
      real*8, dimension(0:N), optional,       intent(out) :: qd_out

      real*8, dimension(0:N) :: Jhe, qd
      real*8, dimension(0:N,n_melem) :: JX

      divq = 0.0d0
      if (present(qd_out)) qd_out = 0.0d0
      if (.not. interdiffusion_enthalpy_active()) return
      call helium_diffusive_face_flux(rho, Tcode, f_sp, Jhe)
      call trace_element_diffusive_face_flux(rho, Tcode, f_sp, JX)
      call interdiffusion_enthalpy_face_flux(Tcode, f_sp, Jhe, qd, JX)
      call interdiffusion_enthalpy_divergence(qd, divq)
      if (present(qd_out)) qd_out = qd

      end subroutine interdiffusion_enthalpy_divergence_of_state

      ! ------------------------------------------------------------------ !

      subroutine reinstate_diffusion_diagnostics(over, under, steps,      &
                                                 resid, trace)
      ! Put the five diagnostics of the last solve back as they were, and
      ! assert it: a measurement that changed what the next report reads
      ! would be a measurement with a history.
      real*8,  intent(in) :: over, under, resid, trace
      integer, intent(in) :: steps
      he_fraction_over_one     = over
      he_fraction_under_zero   = under
      he_fraction_newton_steps = steps
      he_fraction_newton_resid = resid
      trace_ratio_under_zero   = trace
      if (he_fraction_over_one .ne. over .or.                             &
          he_fraction_under_zero .ne. under .or.                          &
          he_fraction_newton_steps .ne. steps .or.                        &
          he_fraction_newton_resid .ne. resid .or.                        &
          trace_ratio_under_zero .ne. trace)                              &
         write(*,'(A)') ' (element transport residual) WARNING: the '//   &
              'module diagnostics were NOT reinstated'
      end subroutine reinstate_diffusion_diagnostics

      ! ------------------------------------------------------------------ !

      subroutine trace_element_transport_coefficients(im, f_sp, rho, TK,   &
                     ntot_phys, ne_phys, yH, mcarH, zb1, eEf, dlnpsi,       &
                     nX, nXold, DcoX, GcoX)
      ! The transport coefficients of ONE trace element against the hydrogen
      ! background: its density and mean charge from the species vector, the
      ! stage-resolved binary diffusion coefficient against the hydrogen
      ! carriers, and the settling coefficient G [1/cm] of the mixing ratio.
      ! The one definition of them: the diffusion step solves the element
      ! with these and element_transport_residual measures the same element's
      ! stationary balance with them.
      integer,                                intent(in)  :: im
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho, TK
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: ntot_phys
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: ne_phys
      real*8, dimension(1-Ng:N+Ng,n_hcar),    intent(in)  :: yH
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: mcarH, zb1
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: eEf, dlnpsi
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nX, nXold
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: DcoX, GcoX

      real*8, dimension(1-Ng:N+Ng) :: zbX
      real*8, dimension(1-Ng:N+Ng,n_mstage) :: yX
      real*8, dimension(n_mstage)  :: ZXs, mXs, alXs
      integer, dimension(n_mstage) :: ispX
      real*8  :: mX
      integer :: i0m, top, j, k

      i0m = melem_i0(im)
      top = melem_top(im)
            mX  = melem_A(im)
            nX  = 0.0d0
            zbX = 0.0d0
            do k = 0, top
               nX  = nX  + f_sp(:,mion_fsp(i0m+k))*rho*n0
               zbX = zbX + dble(k)*f_sp(:,mion_fsp(i0m+k))*rho*n0
            enddo
            where (nX .gt. 1.0d-30)
               zbX = zbX/nX
            elsewhere
               zbX = 0.0d0
            end where
            nXold = nX
            ! Stage-resolved friction against the hydrogen carriers, exactly
            ! as for helium: a metal ion in the ionized wind is held to the
            ! protons by the Coulomb coefficient, which is what keeps it from
            ! settling out (Koskinen et al. 2013, section 3.2.2).  Stages
            ! beyond this element's top are given zero weight.
            ispX = mion_fsp(i0m)
            ZXs  = 0.0d0
            mXs  = mX
            alXs = 0.0d0
            alXs(1) = alpha_melem(im)
            do k = 0, top
               ispX(k+1) = mion_fsp(i0m+k)
               ZXs(k+1)  = dble(k)
            enddo
            call carrier_fractions(f_sp, top+1, ispX(1:top+1),             &
                                   yX(:,1:top+1))
            if (top+1 .lt. n_mstage) yX(:,top+2:n_mstage) = 0.0d0
            call stage_mixture_diffusion(TK, ntot_phys, ne_phys,           &
                                         top+1, yX(:,1:top+1),             &
                                         ZXs(1:top+1), mXs(1:top+1),       &
                                         alXs(1:top+1),                    &
                                         n_hcar, yH, hcar_Z, hcar_m,       &
                                         hcar_alpha, DcoX)
            ! Same ambipolar field as helium, computed from the electron
            ! pressure gradient (memo 3a) -- NOT the hydrogen-plasma constant
            ! eE = m_H g/2 this loop used to assume, which is wrong wherever
            ! helium or the metals carry a significant share of the electrons.
            ! The partner mass is the mean HYDROGEN mass per carrier: a metal
            ! atom settles against the particles it collides with, which are H
            ! atoms and ions above the molecular front and H2/H3+ below it.
            ! mcarH is 1 in the atomic region, so nothing changes there.
            ! -dln(psi)/dr enters here for the same reason it enters the
            ! helium equation: the driver is the metal's MOLE fraction among
            ! the hydrogen carriers, n_X/n_1c, while the transported variable
            ! is the mixing ratio n_X/n_H (module header).
            do j = 1-Ng, N+Ng
               GcoX(j) = ((mX - mcarH(j))*mu*(Dphi(r(j))*v0*v0/R0)   &
                          - (zbX(j) - zb1(j))*eEf(j))/(kb_erg*TK(j))      &
                         - dlnpsi(j)
            enddo
            if (he_alphaT .ne. 0.0d0) then
               do j = 2-Ng, N+Ng-1
                  GcoX(j) = GcoX(j) + he_alphaT*(log(TK(j+1))-log(TK(j-1)))&
                            / max((r(j+1)-r(j-1))*R0, 1.0d0)
               enddo
            endif

      end subroutine trace_element_transport_coefficients

      ! ------------------------------------------------------------------ !

      subroutine trace_face_coefficients(nHl, Dco, Gco, rp, PL, PR)
      ! Face coefficients of the trace-element flux in the mixing ratio,
      !   F(face) = PL fX(j) + PR fX(j+1) ,
      ! with the Peclet hybrid of the helium operator: central where the
      ! settling drift is resolved over a cell, upwind toward the settling
      ! direction otherwise.  The outermost and innermost faces carry no
      ! diffusive flux (PL = PR = 0 there), which is the zero-gradient top
      ! and the Dirichlet reservoir at the base.
      ! The one definition of them: the backward-Euler solve and the
      ! stationary residual of the same element both use these.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: nHl, Dco, Gco, rp
      real*8, dimension(0:N),       intent(out) :: PL, PR
      real*8  :: dr_f, nHf, Df, Gf, DK
      integer :: j
      PL = 0.0d0
      PR = 0.0d0
      do j = 1, N-1
         dr_f = max(rp(j+1)-rp(j), 1.0d0)
         nHf  = 0.5d0*(nHl(j)+nHl(j+1))
         Df   = 0.5d0*(Dco(j)+Dco(j+1))
         Gf   = 0.5d0*(Gco(j)+Gco(j+1))
         DK   = Df + 0.5d0*(kzz_cell(j)+kzz_cell(j+1))
         if (abs(Df*Gf)*dr_f .lt. 2.0d0*DK) then
            PL(j) =  nHf*(DK/dr_f - 0.5d0*Df*Gf)
            PR(j) = -nHf*(DK/dr_f + 0.5d0*Df*Gf)
         else if (Gf .ge. 0.0d0) then          ! settles inward: donor is j+1
            PL(j) =  nHf*(DK/dr_f)
            PR(j) = -nHf*(DK/dr_f + Df*Gf)
         else                                  ! rises: donor is j
            PL(j) =  nHf*(DK/dr_f - Df*Gf)
            PR(j) = -nHf*(DK/dr_f)
         endif
      enddo
      end subroutine trace_face_coefficients

      ! ------------------------------------------------------------------ !

      subroutine trace_composition_residual(fXold, fX, nHl, dt_phys,     &
                                            PL, PR, Frho, wY, cadvX,      &
                                            advect, mres, rnorm,          &
                                            res_out, dsc_out)
      ! Residual of the implicit trace-element step in the MIXING RATIO
      ! fX = n_X/n_H, returned NEGATED (mres is the right-hand side of the
      ! Newton system):
      !
      !   res(j) = (fX(j) - fXold(j))/dt
      !          + [ A_+ F(j) - A_- F(j-1) ] / (n_H V_j)
      !          + advect  div(F_rho Y_X)(j) cadvX(j) ,   Y_X = wY fX .
      !
      ! The diffusive face flux is the two-point form of
      ! trace_face_coefficients; the advective term is the divergence of the
      ! face element mass fluxes (element_advective_divergence), the same
      ! object the helium row balances and the Runge-Kutta stages integrate,
      ! with wY the element mass fraction per unit of the mixing ratio and
      ! cadvX the factor that returns the code divergence to the mixing
      ! ratio per second the row is written in.
      !
      ! rnorm is the residual measured RELATIVE to the size of the terms of
      ! its own row (dsc), for the reason composition_residual states: at a
      ! step long against the cell transport time the terms cancel to many
      ! digits and an absolute residual measures their round-off.
      !
      ! There is no equation at cell 1: it is the Dirichlet reservoir.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: fXold, fX, nHl, dt_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Frho, wY, cadvX
      real*8, dimension(0:N),       intent(in)  :: PL, PR
      logical,                      intent(in)  :: advect
      real*8, dimension(1-Ng:N+Ng), intent(out) :: mres
      real*8,                       intent(out) :: rnorm
      real*8, dimension(1-Ng:N+Ng), intent(out), optional :: res_out, dsc_out

      real*8, dimension(1-Ng:N+Ng) :: Y, dvF, dvM
      ! The one geometry of the grid, in cm: A(f) R0^2 and V(j) R0^3.
      real*8, dimension(0:N)       :: fa
      real*8, dimension(1:N)       :: cv
      real*8  :: R0sq, R0cb
      real*8  :: Kj, sL, sR, Fl, Fr, res, dsc, q
      integer :: j
      logical :: nonfinite_row

      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0
      nonfinite_row = .false.

      dvF = 0.0d0
      dvM = 0.0d0
      if (advect) then
         Y = wY*fX
         call element_advective_divergence(Y, Frho, dvF, dvM)
      endif

      mres  = 0.0d0
      rnorm = 0.0d0
      if (present(res_out)) res_out = 0.0d0
      if (present(dsc_out)) dsc_out = 0.0d0
      do j = 2, N
         Kj  = 1.0d0/(max(nHl(j), 1.0d-30)*cv(j)*R0cb)
         sL  = fa(j-1)*R0sq
         sR  = fa(j)*R0sq
         Fr  = PL(j)*fX(j)     + PR(j)*fX(j+1)
         Fl  = PL(j-1)*fX(j-1) + PR(j-1)*fX(j)
         res = (fX(j) - fXold(j))/dt_phys(j) + Kj*(sR*Fr - sL*Fl)
         dsc = max(abs(fX(j)), abs(fXold(j)))/dt_phys(j)                  &
             + Kj*(sR*abs(Fr) + sL*abs(Fl))
         if (advect) then
            res = res + dvF(j)*cadvX(j)
            dsc = dsc + dvM(j)*cadvX(j)
         endif
         mres(j) = -res
         if (present(res_out)) res_out(j) = res
         if (present(dsc_out)) dsc_out(j) = dsc
         q = abs(res)/max(dsc, 1.0d-300)
         if (q .le. huge(1.0d0)) then
            rnorm = max(rnorm, q)
         else
            nonfinite_row = .true.     ! composition_residual says why
         endif
         ! WHAT THIS ROW BALANCES, term by term, where a caller asked for it
         ! (trace_row_terms_diag).  The two fluxes, the divergence each of
         ! them gives, the mixing ratio on both sides of the cell and the
         ! row's own scale: a row far from zero is read from the sizes of
         ! the terms that make it, and at the outermost cell the outer face
         ! carries no diffusive flux (trace_face_coefficients), so Fr is
         ! zero there and the balance is the inward diffusive flux against
         ! the advective divergence.  A print and nothing else.
         if (trace_row_terms_diag .and. j .ge. trace_row_terms_from)       &
            write(*,'(A,I4,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5,' &
                    //'A,ES12.5,A,ES12.5,A,ES12.5)')                      &
                 ' (element row terms) cell ', j,                         &
                 ': fX(j-1) ', fX(j-1), ', fX(j) ', fX(j),                &
                 ', fX(j+1) ', fX(j+1), ', F_diff(out) ', Fr,             &
                 ', F_diff(in) ', Fl, ', diffusive divergence ',          &
                 Kj*(sR*Fr - sL*Fl), ', advective divergence ',           &
                 dvF(j)*cadvX(j), ', row scale ', dsc
      enddo
      ! A row that is not an ordinary number makes the norm no measurement
      ! (composition_residual states why the test is written out).
      if (nonfinite_row) rnorm = huge(1.0d0)

      end subroutine trace_composition_residual

      ! ------------------------------------------------------------------ !

      subroutine solve_trace_element_in_hydrogen(nX, nHl, Dco, Gco, fbase,  &
                                                 dt_phys, rp, Frho,        &
                                                 wY, cadvX, advect)
      ! One backward-Euler step for a TRACE element diffusing relative to a
      ! hydrogen background nHl, with element diffusion coefficient Dco,
      ! settling coefficient Gco and fixed reservoir base ratio
      ! fbase = nX/nHl.  Valid only for an element whose own mass does not
      ! shape the background it moves through (metal/H ~ 1e-4); helium is
      ! NOT solved this way -- that is what the binary operator above exists
      ! for.
      !
      ! The transported variable is the MIXING RATIO fX = nX/n_H,
      !
      !    dfX/dt + advect div(F_rho Y_X)/(n_H A_X/(n0 msum))
      !           = (1/(n_H r^2)) d/dr [ r^2 n_H
      !                        ( (D_1X + K_zz) dfX/dr + D_1X G fX ) ]
      !
      ! With advect false the row is the diffusive half alone: the element's
      ! own mass fraction is advected on the hydro's face mass fluxes inside
      ! the Runge-Kutta stages, where a cell that is outflowing at both faces
      ! loses exactly the fraction of its mass the mass row loses and cannot
      ! be evacuated of the element.  With advect true the fixed-wind
      ! relaxation carries that same divergence here, formed from the same
      ! face mass flux by element_advective_divergence.
      !
      ! DEFERRED CORRECTION, AND WHY THE LOOP.  The diffusive half is linear
      ! in fX and the matrix below carries it exactly, together with the
      ! DONOR-CELL part of the advective term.  What the reconstruction and
      ! the bounding of species_face_fraction add on top of the donor cell is
      ! carried on the right-hand side at the current iterate, so the fixed
      ! point of the iteration solves the FULL operator -- the same one the
      ! stationary elemental row measures -- and not the first-order one.
      ! With advect false there is nothing deferred, the system is linear,
      ! and the loop takes a single pass which is the direct solve.
      !
      ! nX is intent(inout): supply the current density, receive the solved one.
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: nX
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: nHl, Dco, Gco, dt_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: rp
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: Frho, wY, cadvX
      real*8,                       intent(in)    :: fbase
      logical,                      intent(in)    :: advect

      real*8, dimension(0:N)       :: PL, PR
      real*8, dimension(1:N)       :: advj, advm
      real*8, dimension(1-Ng:N+Ng) :: aa, bb, cc, dd, cpv, dpv, fX, fXold
      real*8, dimension(1-Ng:N+Ng) :: mres, dvD, dvFull, Y, dvFa, dvMa
      ! The one geometry of the grid, in cm: A(f) R0^2 and V(j) R0^3.
      real*8, dimension(0:N)       :: fa
      real*8, dimension(1:N)       :: cv
      real*8 :: R0sq, R0cb
      real*8 :: Kj, mden, sL, sR, cadv, rnorm, rprev
      integer :: j, jd, it
      ! The deferred part is a bounded second-order correction to a term the
      ! matrix already carries, so the iteration contracts quickly; the cap
      ! and the stops mirror solve_mass_fraction.
      integer, parameter :: trace_maxit = 20
      real*8,  parameter :: trace_tol   = 1.0d-12

      fXold = nX/max(nHl, 1.0d-30)
      fX    = fXold
      fX(1-Ng:1)   = fbase
      fX(N+1:N+Ng) = fX(N)

      call trace_face_coefficients(nHl, Dco, Gco, rp, PL, PR)
      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0
      advj = 0.0d0
      advm = 0.0d0
      if (advect) call element_advective_face_coefficients(Frho, advj, advm)

      rprev = huge(1.0d0)
      do it = 1, trace_maxit

         ! The advective term at the current iterate, in full and in the
         ! donor-cell part the matrix below carries.  Their difference is
         ! what the right-hand side defers.
         dvD    = 0.0d0
         dvFull = 0.0d0
         if (advect) then
            Y = wY*fX
            call element_advective_divergence(Y, Frho, dvFa, dvMa)
            do j = 2, N
               dvFull(j) = dvFa(j)*cadvX(j)
               jd = j
               if (Frho(j) .lt. 0.0d0) jd = j + 1
               dvD(j) = advj(j)*wY(jd)*fX(min(jd,N))
               jd = j - 1
               if (Frho(j-1) .lt. 0.0d0) jd = j
               dvD(j) = dvD(j) + advm(j)*wY(jd)*fX(jd)
               dvD(j) = dvD(j)*cadvX(j)
            enddo
         endif

         call trace_composition_residual(fXold, fX, nHl, dt_phys,          &
                                         PL, PR, Frho, wY, cadvX, advect,  &
                                         mres, rnorm)

         do j = 2, N
            Kj    = 1.0d0/(max(nHl(j), 1.0d-30)*cv(j)*R0cb)
            sL    = fa(j-1)*R0sq
            sR    = fa(j)*R0sq
            aa(j) = -Kj*sL*PL(j-1)
            bb(j) =  1.0d0/dt_phys(j) + Kj*(sR*PL(j) - sL*PR(j-1))
            cc(j) =  Kj*sR*PR(j)
            dd(j) =  fXold(j)/dt_phys(j)
            if (advect) then
               ! The deferred correction: the full divergence minus the part
               ! the matrix carries, both at the current iterate.
               dd(j) = dd(j) - (dvFull(j) - dvD(j))
               ! donor-cell entries of the two faces
               cadv = advj(j)*cadvX(j)
               jd   = j
               if (Frho(j) .lt. 0.0d0) jd = j + 1
               if (jd .le. j) then
                  bb(j) = bb(j) + cadv*wY(jd)
               else if (j .lt. N) then
                  cc(j) = cc(j) + cadv*wY(jd)
               else
                  ! zero-gradient outer ghost: fX(N+1) = fX(N)
                  bb(j) = bb(j) + cadv*wY(jd)
               endif
               cadv = advm(j)*cadvX(j)
               jd   = j - 1
               if (Frho(j-1) .lt. 0.0d0) jd = j
               if (jd .lt. j) then
                  aa(j) = aa(j) + cadv*wY(jd)
               else
                  bb(j) = bb(j) + cadv*wY(jd)
               endif
            endif
         enddo
         cc(N) = 0.0d0

         dd(2) = dd(2) - aa(2)*fbase                      ! base Dirichlet
         aa(2) = 0.0d0

         cpv(2) = cc(2)/bb(2)
         dpv(2) = dd(2)/bb(2)
         do j = 3, N
            mden   = bb(j) - aa(j)*cpv(j-1)
            cpv(j) = cc(j)/mden
            dpv(j) = (dd(j) - aa(j)*dpv(j-1))/mden
         enddo
         fX(N) = dpv(N)
         do j = N-1, 2, -1
            fX(j) = dpv(j) - cpv(j)*fX(j+1)
         enddo
         fX(1-Ng:1)   = fbase
         fX(N+1:N+Ng) = fX(N)

         if (.not. advect) exit             ! linear: one pass is the solve
         if (rnorm .le. trace_tol) exit
         if (it .ge. 3 .and. rnorm .gt. 0.5d0*rprev) exit
         rprev = rnorm
      enddo

      ! Clip at zero.  Unlike the helium equation this one is linear in fX --
      ! the trace limit drops the (1 - fX) factor -- and the matrix has the
      ! M-matrix sign pattern with a nonnegative right-hand side, which is the
      ! argument for fX >= 0.  It is an argument and not a measurement, so how
      ! far the solve leaves the range is recorded (trace_ratio_under_zero)
      ! and reported with the step diagnostic rather than asserted.
      trace_ratio_under_zero = max(trace_ratio_under_zero,                &
                                   -minval(fX(1:N)))
      where (fX .lt. 0.0d0) fX = 0.0d0
      nX = fX*nHl

      end subroutine solve_trace_element_in_hydrogen

      ! ------------------------------------------------------------------ !

      logical function species_advection_active()
      ! Whether any species mass fraction is transported on the face mass
      ! fluxes in this configuration.  With it false the transported set is
      ! empty, the Runge-Kutta stages carry no species row, and the run is the
      ! one a build without this operator gives.
      species_advection_active = element_rows_advected() .or.             &
                                 carrier_rows_advected()
      end function species_advection_active

      ! ------------------------------------------------------------------ !

      logical function element_rows_advected()
      ! The element mass fractions are transported exactly where the element
      ! diffusion operator transports them.
      element_rows_advected = he_diffusion .and. thereis_He
      end function element_rows_advected

      ! ------------------------------------------------------------------ !

      logical function carrier_rows_advected()
      ! The carriers are transported exactly where the photochemical
      ! transport operator declared them.
      carrier_rows_advected = (n_car .gt. 0)
      end function carrier_rows_advected

      ! ------------------------------------------------------------------ !

      subroutine advected_carrier_reset()
      ! Empty the declared set.  The declaration is made once per run, after
      ! the input keys are parsed, so emptying first makes a second
      ! declaration of the same run's set the same set and not twice it.
      n_car     = 0
      car_ready = .false.
      end subroutine advected_carrier_reset

      ! ------------------------------------------------------------------ !

      subroutine advected_carrier_register(isp, base_imposed)
      ! Declare one carrier species of the transported set.  isp is its
      ! column in f_sp; base_imposed says that a lower-atmosphere handoff
      ! states its composition at the base, so the inflowing face carries the
      ! ghost value rather than the base cell's own.  The mass is read from
      ! the species table, so a carrier weighs here exactly what the equation
      ! of state weighs it.
      integer, intent(in) :: isp
      logical, intent(in) :: base_imposed
      integer :: ib, ipos

      ipos = 0
      do ib = 1, n_bsp
         if (bsp_fsp(ib) .eq. isp) ipos = ib
      enddo
      if (ipos .eq. 0) then
         write(*,*) 'ERROR: advected_carrier_register: species ', isp,    &
                    ' is not in the base species table'
         error stop 1
      endif
      if (n_car .ge. n_car_max) then
         write(*,*) 'ERROR: advected_carrier_register: more than ',       &
                    n_car_max, ' carriers'
         error stop 1
      endif
      n_car = n_car + 1
      car_isp(n_car)          = isp
      car_mass(n_car)         = bsp_mass(ipos)
      car_base_imposed(n_car) = base_imposed

      end subroutine advected_carrier_register

      ! ------------------------------------------------------------------ !

      integer function advected_carrier_count()
      advected_carrier_count = n_car
      end function advected_carrier_count

      ! ------------------------------------------------------------------ !

      logical function advected_carrier_state_ready()
      ! Whether the carriers of the step just taken were advected and are
      ! waiting to be read.  False through an attempt that never reached its
      ! projection, so a reader falls back on the species vector.
      advected_carrier_state_ready = car_ready
      end function advected_carrier_state_ready

      ! ------------------------------------------------------------------ !

      subroutine advected_carrier_fractions(f_sp, fc)
      ! The advected carriers, in the units f_sp counts them in: a number
      ! per unit mass, n_c/(rho n0).  The mass fraction the faces carried is
      ! divided by the mixture mass of the state it is being written onto, so
      ! that a carrier means the same amount of gas before and after any
      ! operator that has moved the mixture in between.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng,n_car),     intent(out) :: fc

      real*8, dimension(1-Ng:N+Ng) :: msum
      integer :: k

      fc = 0.0d0
      if (.not. car_ready) return
      call mixture_mass_sum(f_sp, msum)
      do k = 1, n_car
         fc(:,k) = Ycar(:,k)*msum/car_mass(k)
      enddo

      end subroutine advected_carrier_fractions

      ! ------------------------------------------------------------------ !

      subroutine advected_carrier_state_consumed()
      ! The advected composition has been taken up by the operator that
      ! transports the carriers.  It is a state of ONE step, so it is marked
      ! spent here rather than left to be read a second time by a path that
      ! runs between two steps.
      car_ready = .false.
      end subroutine advected_carrier_state_consumed

      ! ------------------------------------------------------------------ !

      subroutine mixture_mass_sum(f_sp, msum)
      ! Mass of the mixture per unit of f_sp, = rho/(n0 mu).  It is the
      ! denominator every mass fraction in this module is taken against.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: msum

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, mass1

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)

      end subroutine mixture_mass_sum

      ! ------------------------------------------------------------------ !

      subroutine mixture_mass_split(f_sp, nucH, nucHe, mass1, msum, mmet)
      ! THE MASS A CELL'S SPECIES ACTUALLY CARRY, split into the two
      ! components of the module header:
      !
      !   mass1 = m_H n_H + sum_X A_X n_X + A_O n_O(mol) + A_C n_C(mol)
      !   msum  = mass1 + m_He n_He  = rho/(n0 mu)
      !
      ! component 1 being the hydrogen carriers together with the trace
      ! metals and the heavy nuclei bound in the molecular carriers, and
      ! component He the helium of every stage (the helium nucleus of HeH+
      ! included, its proton being component 1's).
      !
      ! IT IS THE STATE'S OWN MASS AND NOT THE MASS A RESERVOIR COMPOSITION
      ! WOULD CARRY.  The metals enter at the abundance the cell holds, not
      ! at melem_ab: the trace transport moves them against hydrogen on
      ! purpose, the element rows advect each of them on its own face flux,
      ! and a restart may carry a settled column, so m_1 n_H + m_He n_He
      ! with the reservoir m_1 is not this cell's mass.  Every mass fraction
      ! of this module is taken against msum and project_elements inverts
      ! THIS split, which is what makes the species handed back carry the
      ! density calc_rho reads out of them.
      !
      ! No species weight is restated here: a base species mass is
      ! n_H m_H + n_He m_He + n_O A_O + n_C A_C exactly by the construction
      ! of bsp_mass (species_table), so the split is that table's own.  The
      ! metal-ion mass follows the eos_metals policy calc_rho follows.  With
      ! the metals off and no oxygen chemistry mass1 is n_H exactly (m_H is
      ! the mass unit), so the atomic and molecular hydrogen columns are
      ! untouched by the split.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nucH, nucHe
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: mass1, msum
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: mmet

      real*8, dimension(1-Ng:N+Ng) :: mmetal, nucOm, nucCm
      integer :: ib, im

      call element_nucleus_counts(f_sp, nucH, nucHe)

      mmetal = 0.0d0
      if (eos_include_metals .and. thereis_metals) then
         do im = 1, n_mion
            mmetal = mmetal + melem_A(mion_elem(im))*f_sp(:,mion_fsp(im))
         enddo
      endif
      if (present(mmet)) mmet = mmetal

      ! The oxygen and carbon bound in OH, H2O and CO: the ionization solve
      ! removes them from the metal-ion totals, so they are counted here and
      ! nowhere else (the same non-double-counting calc_rho makes).
      nucOm = 0.0d0
      nucCm = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nO(ib) .gt. 0)                                           &
            nucOm = nucOm + dble(bsp_nO(ib))*f_sp(:,bsp_fsp(ib))
         if (bsp_nC(ib) .gt. 0)                                           &
            nucCm = nucCm + dble(bsp_nC(ib))*f_sp(:,bsp_fsp(ib))
      enddo

      mass1 = m_H_amu*nucH + mmetal                                       &
              + melem_A(iel_O)*nucOm + melem_A(iel_C)*nucCm
      msum  = mass1 + m_He_amu*nucHe
      where (msum .lt. 1.0d-30) msum = 1.0d-30

      end subroutine mixture_mass_split

      ! ------------------------------------------------------------------ !

      subroutine element_mass_fractions(f_sp, Y)
      ! The transported element mass fractions of a species vector: helium in
      ! column 1, each trace metal in column 1+im when trace-metal diffusion
      ! is on.
      ! The mass sum is the cell's own (mixture_mass_split), the two-component
      ! closure of the module header, so column 1 is exactly the X the
      ! diffusive half solves and the columns are true mass fractions of it.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng,1+n_melem), intent(out) :: Y

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, mass1, nucX
      integer :: j, im, i0m, k

      Y   = 0.0d0
      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      Y(:,1) = m_He_amu*nucHe/msum

      if (.not. (he_metal_diffusion .and. thereis_metals)) return
      do im = 1, n_melem
         i0m  = melem_i0(im)
         nucX = 0.0d0
         do k = 0, melem_top(im)
            do j = 1-Ng, N+Ng
               nucX(j) = nucX(j) + f_sp(j,mion_fsp(i0m+k))
            enddo
         enddo
         Y(:,1+im) = melem_A(im)*nucX/msum
      enddo

      end subroutine element_mass_fractions

      ! ------------------------------------------------------------------ !

      subroutine species_advection_begin_step(f_sp)
      ! Read the composition the step begins with out of the species vector,
      ! for every transported row: the element mass fractions and the
      ! declared carriers.  Called once per attempt, before the first
      ! Runge-Kutta stage, so that a retaken attempt starts from the
      ! composition the checkpoint restored and not from the one a discarded
      ! attempt advected.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp

      car_ready = .false.
      step_flux_ready = .false.

      if (element_rows_advected()) then
         n_etr = 1
         if (he_metal_diffusion .and. thereis_metals) n_etr = 1 + n_melem
         if (.not. allocated(Yetr))                                       &
            allocate(Yetr(1-Ng:N+Ng,1+n_melem), Yetr_n(1-Ng:N+Ng,1+n_melem))
         call element_mass_fractions(f_sp, Yetr)
         Yetr_n = Yetr
      endif

      if (carrier_rows_advected()) then
         if (.not. allocated(Ycar))                                       &
            allocate(Ycar(1-Ng:N+Ng,n_car), Ycar_n(1-Ng:N+Ng,n_car))
         call carrier_mass_fractions(f_sp, Ycar)
         Ycar_n = Ycar
      endif

      end subroutine species_advection_begin_step

      ! ------------------------------------------------------------------ !

      subroutine carrier_mass_fractions(f_sp, Y)
      ! The declared carriers as mass fractions of the mixture, with the
      ! inflow composition of the inner ghosts: the handoff value where one
      ! is stated, the base cell's own partition where none is; and the
      ! outer ghosts continuing cell N's own value, the outflow rule the
      ! stationary carrier operator states in carrier_face_mass_fraction
      ! (diffusive_photochemistry): the outer face is an outflow face of a
      ! transported scalar, its characteristics leave, and the equilibrium
      ! composition the sweep solves in the ghost cell is not what the wind
      ! carried out of cell N. One rule for the stage path and the
      ! stationary path, so the two evaluations of a carrier row agree.
      !
      ! X_ghost = X_N HOLDS IN BOTH DIRECTIONS OF THE FACE MASS FLUX. Where
      ! F_rho(N) reverses there is no reservoir at this end to state an
      ! inflow composition with: the hydrodynamic closure at the same face
      ! continues cell N's own state outward and holds no composition of its
      ! own, so the gas that comes back carries what left. With the ghost
      ! equal to cell N the forward difference of the MC limiter vanishes,
      ! the limited slope of cell N is zero, and the outer face is the donor
      ! cell average itself whichever way the flux points.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng,n_car),     intent(out) :: Y

      real*8, dimension(1-Ng:N+Ng) :: msum
      integer :: k, j

      call mixture_mass_sum(f_sp, msum)
      do k = 1, n_car
         Y(:,k) = car_mass(k)*f_sp(:,car_isp(k))/msum
         if (.not. car_base_imposed(k)) then
            do j = 1-Ng, 0
               Y(j,k) = Y(1,k)
            enddo
         endif
         do j = N+1, N+Ng
            Y(j,k) = Y(N,k)
         enddo
      enddo
      where (Y .lt. 0.0d0) Y = 0.0d0

      end subroutine carrier_mass_fractions

      ! ------------------------------------------------------------------ !

      subroutine species_advection_stage(stage, rho_n, rho_in, rho_new,   &
                                         Frho, dt_loc)
      ! One Runge-Kutta stage of the transported species rows, on the face
      ! mass fluxes, the volumes and the time step of the mass row of the
      ! same stage.
      !
      ! Element rows 1-Ng .. 1 are not advanced: cell 1 is the Dirichlet
      ! reservoir of the element operator, and the inner ghosts carry the same
      ! composition, so the inflowing face at the base of cell 2 reconstructs
      ! the reservoir composition and the base inflow is imposed on the face
      ! flux itself.  The outer ghosts take the zero-gradient copy the outer
      ! boundary uses for every other quantity.
      integer, intent(in) :: stage
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho_n, rho_in, rho_new
      real*8, dimension(1-Ng:N+Ng), intent(in) :: Frho, dt_loc

      real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Ynew
      real*8, dimension(1-Ng:N+Ng,max(n_car,1)) :: Ycnew
      integer :: q

      if (element_rows_advected()) then
         call species_advective_update(stage, n_etr, 1,                   &
                                    Yetr_n(:,1:n_etr), Yetr(:,1:n_etr),   &
                                    Frho, rho_n, rho_in, rho_new, dt_loc, &
                                    2, Ynew(:,1:n_etr))
         Yetr(:,1:n_etr) = Ynew(:,1:n_etr)
         do q = 1, n_etr
            Yetr(N+1:N+Ng,q) = Yetr(N,q)
         enddo
         ! The step's face mass flux, with the SSP-RK3 weights of the stage
         ! (declaration of Frho_step).
         if (allocated(Frho_step)) then
            if (size(Frho_step) .ne. N + 2*Ng) deallocate(Frho_step)
         endif
         if (.not. allocated(Frho_step)) allocate(Frho_step(1-Ng:N+Ng))
         select case (stage)
         case (1)
            Frho_step = Frho/6.0d0
         case (2)
            Frho_step = Frho_step + Frho/6.0d0
         case (3)
            Frho_step = Frho_step + 2.0d0*Frho/3.0d0
            step_flux_ready = .true.
         end select
      endif

      ! THE CARRIERS ARE ADVANCED FROM CELL 1, and the elements are not.
      ! Cell 1 is the Dirichlet reservoir of the element operator, so its
      ! element composition is boundary data; the partition of an element
      ! among its carriers there is not data but a state the base face
      ! carries into, so cell 1 has an advective term like every other cell
      ! and the inner ghosts supply the composition that flows in.
      if (carrier_rows_advected()) then
         call species_advective_update(stage, n_car, 0,                   &
                                    Ycar_n, Ycar,                         &
                                    Frho, rho_n, rho_in, rho_new, dt_loc, &
                                    1, Ycnew(:,1:n_car))
         Ycar = Ycnew(:,1:n_car)
         do q = 1, n_car
            Ycar(N+1:N+Ng,q) = Ycar(N,q)
         enddo
      endif

      end subroutine species_advection_stage

      ! ------------------------------------------------------------------ !

      subroutine species_advection_project(f_sp)
      ! Write the advected element composition into the species vector: the
      ! element totals are set here, the split within an element is left to
      ! the ionization solve, which is the same division of labour
      ! project_elements states.
      !
      ! THE CARRIERS ARE NOT WRITTEN HERE.  A carrier is one species inside
      ! an element, and moving it in the species vector without moving the
      ! stages it competes with would leave the element total wrong for as
      ! long as the two are apart.  Its advected fraction is therefore handed
      ! to the photochemical transport operator as the state that operator's
      ! step starts from, and that operator's write-back restores the element
      ! totals cell by cell, which is where the two are put back together.
      !
      ! Each metal element is put on the amount its own advected mass fraction
      ! asks for, and the hydrogen group closes the mixture around it.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp

      if (carrier_rows_advected()) car_ready = .true.
      if (.not. element_rows_advected()) return
      call project_element_mass_fractions(f_sp, Yetr)

      end subroutine species_advection_project

      ! ------------------------------------------------------------------ !

      subroutine project_element_mass_fractions(f_sp, Y)
      ! Write the element TOTALS of a species vector from the element mass
      ! fractions Y: helium in column 1, each trace metal in column 1+im.
      ! It is the inverse of element_mass_fractions, and it is the one map
      ! from an element mass fraction to a species vector -- the advective
      ! stages and any stationary system that carries an element as an
      ! unknown both go through here, so the two cannot disagree about what
      ! an element fraction means.
      !
      ! The split WITHIN an element is left to the ionization solve, which is
      ! the division of labour project_elements states.  With trace-metal
      ! diffusion on, each metal element is put on the amount its own advected mass
      ! fraction asks for and the hydrogen group takes the mass that is left,
      ! so the mixture's mass is unchanged by the write-back.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng,1+n_melem), intent(in)    :: Y

      real*8, dimension(1-Ng:N+Ng) :: nucH, nucHe, msum, mass1, nucX
      real*8  :: nXnew, rXsc
      integer :: j, im, i0m, k, top

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)

      if (.not. (he_metal_diffusion .and. thereis_metals)) then
         call project_elements(f_sp, Y(:,1), msum, .true.)
         return
      endif

      ! THE METALS ARE WRITTEN FIRST AND THE HYDROGEN GROUP CLOSES THE
      ! MIXTURE.  Each metal element carries its own advected mass fraction
      ! here, so it is not a ratio to hydrogen and does not follow hydrogen's
      ! factor; putting it on its target first and letting project_elements
      ! give the hydrogen group the mass that is left is what keeps
      ! sum_i m_i n_i at msum.  Writing hydrogen first and the metals after
      ! it added the metals' change of mass to the mixture.
      do im = 1, n_melem
         i0m  = melem_i0(im)
         top  = melem_top(im)
         nucX = 0.0d0
         do k = 0, top
            do j = 1-Ng, N+Ng
               nucX(j) = nucX(j) + f_sp(j,mion_fsp(i0m+k))
            enddo
         enddo
         do j = 1-Ng, N+Ng
            nXnew = max(Y(j,1+im), 0.0d0)*msum(j)/melem_A(im)
            if (nucX(j) .gt. 1.0d-25*max(nucH(j), 1.0d-30)) then
               rXsc = nXnew/nucX(j)
               do k = 0, top
                  f_sp(j,mion_fsp(i0m+k)) = f_sp(j,mion_fsp(i0m+k))*rXsc
               enddo
            else if (nXnew .gt. 1.0d-25*max(nucH(j), 1.0d-30)) then
               ! The element returned to a cell it had been emptied of: it
               ! comes back in the neutral ground stage, as every other
               ! re-seed in this module does, and the ionization solve
               ! re-splits it on the next call.
               f_sp(j,mion_fsp(i0m)) = nXnew
               do k = 1, top
                  f_sp(j,mion_fsp(i0m+k)) = 0.0d0
               enddo
            endif
         enddo
      enddo

      call project_elements(f_sp, Y(:,1), msum, .false.)

      end subroutine project_element_mass_fractions

      ! ------------------------------------------------------------------ !

      logical function diffusion_check_on()
      ! EXHALE_DIFFUSION_CHECK=1 turns on the stderr diagnostic written each
      ! step, and the elemental face flux profile of every state the
      ! certification measures (element_flux_profile_on).
      character(len=8) :: envv
      call get_environment_variable('EXHALE_DIFFUSION_CHECK', envv)
      diffusion_check_on = (trim(envv) .eq. '1')
      end function diffusion_check_on

      ! ------------------------------------------------------------------ !

      logical function element_flux_profile_on()
      ! Whether the elemental face flux profile of a measured state is
      ! written: on demand (EXHALE_DIFFUSION_CHECK=1), and unconditionally
      ! whenever a lower-atmosphere profile is in use, because there the
      ! elemental fluxes are not a diagnostic but the quantity the two models
      ! have to agree on, so the run must always leave them behind.
      element_flux_profile_on = diffusion_check_on() .or. lap_in_use
      end function element_flux_profile_on

      ! ------------------------------------------------------------------ !

      subroutine face_mass_flux_cgs(Frho, msum, Fcgs)
      ! A face mass flux in the code units of the mass row, returned in
      ! g cm^-2 s^-1:
      !
      !    F_cgs(f) = F_rho(f) n0 mu v0 [ msum(f) + msum(f+1) ]/2 .
      !
      ! n0 mu msum is the mass one unit of the code density carries (the
      ! factor every row of this operator converts with, cadvf) and v0 the
      ! velocity unit.  The row of cell j converts the flux of both its
      ! faces with its own msum(j), so a face shared by two cells is
      ! converted with two values; the face value here is their mean.  The
      ! two differ by the state's mass closure, msum = 1 to 1e-15 on a
      ! closed state, so the choice moves no digit a record prints.  The
      ! outermost ghost face has no cell beyond it and takes its own cell's
      ! value.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Frho, msum
      real*8, dimension(1-Ng:N+Ng), intent(out) :: Fcgs
      integer :: j, jn
      do j = 1-Ng, N+Ng
         jn = min(j+1, N+Ng)
         Fcgs(j) = Frho(j)*n0*mu*v0*0.5d0*(msum(j) + msum(jn))
      enddo
      end subroutine face_mass_flux_cgs

      ! ------------------------------------------------------------------ !

      subroutine write_element_flux_profile(rho, Tcode, f_sp, Frho)
      ! Radial profile of the elemental face fluxes of ONE STATE, written to
      ! ./output/element_flux_profile.txt and replaced at every call.  The
      ! caller is the certification, which measures every state the run
      ! writes -- the loaded state of an evaluation, the returned state of a
      ! stationary solve and the final state of a march -- so after a run the
      ! file holds the state the run wrote.  The fluxes it writes are:
      !
      !   F_He(r_f) = 4 pi r_f^2 ( F_rho Y_He + J )         [g/s]
      !   F_H (r_f) = 4 pi r_f^2 ( F_rho (1 - Y_He) - J )   [g/s]
      !
      ! A binary mixture has ONE independent diffusive flux, so the hydrogen
      ! element carries -J against the helium element's +J, and the two
      ! advective parts add up to the face mass flux.  Both halves come from
      ! element_nucleus_face_flux, the one public spelling of the element
      ! flux: Y_He(f) is therefore the reconstructed, limited and upwinded
      ! face mass fraction the transport rows put on the face, and not a
      ! one-sided pick of a cell value.  The total mass flux
      ! 4 pi r_f^2 F_rho is written beside them: at a steady state all three
      ! are constant with radius, and the comparison of their radial spreads
      ! is the pass criterion.  dmeff is the face-averaged relative settling
      ! mass of the ambipolar diagnostic (3 neutral, 2.5 H+ plasma,
      ! 5/3 He++ plasma).
      !
      ! THE FACE MASS FLUX THE PROFILE IS WRITTEN ON is Frho, the Riemann
      ! face mass flux of this state that its own mass row differences
      ! (face_mass_flux_of_state), the one the elemental row of the same
      ! state rides on; face_mass_flux_cgs puts it in g cm^-2 s^-1.  The face
      ! mean of the cell-centred rho and v is not that flux: near the base
      ! the cell-centred velocity carries a collocated two-cell odd-even mode
      ! the Riemann flux does not carry, which at faces 1 to 4 of the
      ! certified LHS 1140 b kzz1e9 He/H 0.55 state puts the face mean at
      ! -1.07e9, +1.16e8, -3.8e7 and +5.9e7 g/s against a wind the Riemann
      ! flux carries at 3.93e7 g/s through every face (MEASURED
      ! 2026-09-24).
      !
      ! With a lower-atmosphere profile in use the routine also reduces each
      ! elemental flux to ONE number -- its median and its relative radial
      ! spread -- over each of TWO radial windows, for write_resolved_config
      ! and the flux closure.
      !
      ! (a) The OVERLAP window, the interval both models describe.  Its upper
      !     edge is the radius the profile reaches, r(p_top).  Its lower edge
      !     is MEASURED rather than assumed: the first face at which the
      !     5-face moving spread of the face mass flux 4 pi r^2 F_rho falls
      !     below 10 percent, which leaves out the faces over which the mass
      !     flux of the state is not yet constant -- a base transient of a
      !     state that is still marching.  On a stationary state no face is
      !     left out: the first face the 5-face test can be centred on
      !     qualifies.  When no face qualifies the edge falls back to
      !     1.02 R_p and lap_flux_r_lo_measured says so.  An empty window is
      !     recorded as empty and never replaced by a number from outside it.
      !
      ! (b) The STEADY-FLUX window, r >= r_esc, i.e. the [j_min:N] escape
      !     region the solver uses to declare the wind steady.  At a steady
      !     state the elemental flux 4 pi r^2 (F_rho Y + J) is independent of
      !     radius, so the flux measured there IS the flux through the
      !     matching level.  That identity is why this window is a legitimate
      !     statement about the handoff and not a different quantity, and it
      !     is the only window wide enough on a planet whose lower-atmosphere
      !     column spans ~0.01 R_p.
      !
      ! The last three columns are the stage-resolved friction of memo 2.6:
      ! the effective coefficient D_eff of the cell, the all-neutral
      ! hard-sphere coefficient of the same cell, their ratio (which is the
      ! Coulomb suppression where both elements are ionized), and the name of
      ! the stage pair carrying the largest share of the friction.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      ! The Riemann face mass flux of this state, in the code units of the
      ! mass row.
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: Frho
      real*8, dimension(1-Ng:N+Ng) :: Frho_f, nucH, nucHe, mass1, msum, Xhe
      real*8, dimension(1-Ng:N+Ng) :: Yf, rep, Dco, Dneut, Gco, dmeff
      real*8, dimension(1-Ng:N+Ng) :: zb1, eEf, dlnpsi
      integer, dimension(1-Ng:N+Ng):: idom
      real*8, dimension(0:N)       :: Fadv, Jdif, n_el, m_one
      real*8  :: area, adv, adv_H, dmf, Df, Dnf, sv_exc
      real*8, dimension(1:N) :: FH_win, FHe_win, M_win, r_win
      integer :: j, uu, is, it, nwin, sv_bnd

      if (.not. thereis_He) return
      nwin = 0

      ! THE FACE COUNTERS OF THE RUN COUNT THE RUN'S OWN FACES.  The
      ! evaluation below reconstructs a face composition, and
      ! species_face_fraction counts every face state it has to scale back
      ! onto [0,1] and keeps the largest excursion it saw.  Those are
      ! statements about the trajectory, so a diagnostic that reads the same
      ! state a second time must not add to them.
      sv_bnd = n_species_faces_bounded
      sv_exc = species_face_excursion

      call mixture_mass_split(f_sp, nucH, nucHe, mass1, msum)
      Xhe = m_He_amu*nucHe/msum
      rep = r_edg*R0

      ! The face mass flux of the state in g cm^-2 s^-1.
      call face_mass_flux_cgs(Frho, msum, Frho_f)

      ! Both halves of the element flux, from the one routine that spells
      ! them: the advective half is Frho_f times the face mass fraction, so
      ! it comes back in the units Frho_f was given in.
      call element_nucleus_face_flux(rho, Tcode, f_sp, Frho_f,            &
                                     Fadv, Jdif, n_el, m_one)

      ! The face mass fraction the advective half carries, for the column
      ! that reports it: the same rule, from the same routine the flux above
      ! was built with.
      call species_face_fraction(Xhe, Frho_f, Yf)

      n_species_faces_bounded = sv_bnd
      species_face_excursion  = sv_exc

      ! The friction columns: the coefficients the operator itself forms.
      call helium_hydrogen_diffusion(rho, Tcode, f_sp, Dco, Dneut, idom)
      call settling_coefficient(rho, Tcode, f_sp, Gco, dmeff, zb1, eEf,   &
                                dlnpsi)

      uu = 771
      open(unit=uu, file='./output/element_flux_profile.txt',             &
           status='replace')
      write(uu,'(A)') '# elemental face fluxes of the state the '//       &
         'certification measured'
      write(uu,'(A)') '# (both halves from element_nucleus_face_flux, '// &
         'on the Riemann face mass flux of the state)'
      write(uu,'(A)') '#   F_rho  = Riemann face mass flux            '// &
         '[g/cm2/s]  the one the mass row differences'
      write(uu,'(A)') '#   X_face = face mass fraction of helium '//      &
         '           reconstructed, limited, upwind on F_rho'
      write(uu,'(A)') '#   F_He   = 4 pi r^2 (F_rho X_face + J)      [g/s]'
      write(uu,'(A)') '#   F_H    = 4 pi r^2 (F_rho (1-X_face) - J)  [g/s]'
      write(uu,'(A)') '#   Mdot_f = 4 pi r^2 F_rho                   '//  &
         '[g/s]   (= F_H + F_He identically)'
      ! One whitespace-free token per column, so the schema line can be read
      ! the way every other EXHALE product's is.
      write(uu,'(A)') '# columns: j r_face[R_p] X_face J_diff'//           &
         '[g/cm2/s] F_adv[g/cm2/s] F_He[g/s] F_H[g/s]'//                  &
         ' Mdot_face[g/s] dmeff_face'//                                   &
         ' D_eff[cm2/s] D_neutral[cm2/s] D_eff/D_neutral dominant_pair'
      do j = 1, N-1
         area  = 4.0d0*pi*rep(j)**2
         adv   = Fadv(j)
         ! The hydrogen element carries the rest of the face mass flux, so
         ! that F_H + F_He is the mass flux to the last bit.
         adv_H = Frho_f(j) - adv
         dmf  = 0.5d0*(dmeff(j) + dmeff(j+1))
         Df   = 0.5d0*(Dco(j)   + Dco(j+1))
         Dnf  = 0.5d0*(Dneut(j) + Dneut(j+1))
         is   = (idom(j) - 1)/n_hcar + 1
         it   = idom(j) - (is-1)*n_hcar
         write(uu,'(I6,11ES16.7,2X,A)') j, rep(j)/R0, Yf(j),              &
              Jdif(j), adv, area*(adv + Jdif(j)),                         &
              area*(adv_H - Jdif(j)), area*Frho_f(j), dmf,                &
              Df, Dnf, Df/max(Dnf, 1.0d-99),                              &
              trim(hecar_name(is))//'-'//trim(hcar_name(it))
         if (lap_in_use) then
            nwin = nwin + 1
            r_win(nwin)   = rep(j)
            FHe_win(nwin) = area*(adv + Jdif(j))
            FH_win(nwin)  = area*(adv_H - Jdif(j))
            M_win(nwin)   = area*Frho_f(j)
         endif
      enddo
      close(uu)

      if (lap_in_use) call reduce_element_flux_windows(nwin, r_win, FH_win,&
                                                       FHe_win, M_win)

      end subroutine write_element_flux_profile

      ! ------------------------------------------------------------------ !

      subroutine reduce_element_flux_windows(nf, r_f, FH_f, FHe_f, M_f)
      ! Reduce the face-resolved elemental fluxes to one median and one
      ! relative spread over each of the two windows described above, and
      ! leave them in the lower_atmosphere_profile module for
      ! write_resolved_config.  Nothing here recomputes a flux: it only
      ! selects faces and takes order statistics of the state's own fluxes.
      integer, intent(in) :: nf
      real*8, dimension(nf), intent(in) :: r_f, FH_f, FHe_f, M_f
      real*8, dimension(nf) :: a
      real*8  :: r_lo, r_hi, sprd, med
      integer :: i, i1, i2, n
      ! Half-width of the moving window the flatness test is taken over, and
      ! the spread below which the face mass flux is called constant.
      integer, parameter :: nhalf = 2
      real*8,  parameter :: spread_flat = 0.10d0

      lap_flux_measured = .true.

      ! ---- (a) overlap window ------------------------------------------- !
      ! Lower edge: the first face whose 5-face moving spread of the face
      ! mass flux is below spread_flat.  M_f is 4 pi r^2 F_rho, and a ratio
      ! of a max-minus-min to a median is insensitive to that constant.
      r_lo = -1.0d0
      lap_flux_r_lo_measured = .false.
      do i = 1 + nhalf, nf - nhalf
         med  = median_of(M_f(i-nhalf:i+nhalf))
         sprd = (maxval(M_f(i-nhalf:i+nhalf))                             &
                 - minval(M_f(i-nhalf:i+nhalf)))                          &
                /max(abs(med), 1.0d-99)
         if (sprd .lt. spread_flat) then
            r_lo = r_f(i)
            lap_flux_r_lo_measured = .true.
            exit
         endif
      enddo
      if (.not. lap_flux_r_lo_measured) r_lo = 1.02d0*R0

      r_hi = -1.0d0
      if (lap_r_top_RJ .gt. 0.0d0) r_hi = lap_r_top_RJ*RJ

      lap_flux_r_lo_Rp = r_lo/R0
      lap_flux_r_hi_Rp = r_hi/R0

      i1 = 0; i2 = -1
      do i = 1, nf
         if (r_f(i) .ge. r_lo .and. r_f(i) .le. r_hi) then
            if (i1 .eq. 0) i1 = i
            i2 = i
         endif
      enddo
      n = 0
      if (i1 .gt. 0) n = i2 - i1 + 1
      lap_flux_nface        = n
      lap_flux_window_empty = (n .lt. 1)
      if (n .ge. 1) then
         a(1:n) = FH_f(i1:i2)
         call median_and_spread(a(1:n), lap_FH_median,  lap_FH_spread)
         a(1:n) = FHe_f(i1:i2)
         call median_and_spread(a(1:n), lap_FHe_median, lap_FHe_spread)
         a(1:n) = M_f(i1:i2)
         call median_and_spread(a(1:n), lap_Mdot_median, lap_Mdot_spread)
      endif

      ! ---- (b) steady-flux window, r >= r_esc --------------------------- !
      i1 = 0; i2 = -1
      do i = 1, nf
         if (r_f(i) .ge. r_esc*R0) then
            if (i1 .eq. 0) i1 = i
            i2 = i
         endif
      enddo
      n = 0
      if (i1 .gt. 0) n = i2 - i1 + 1
      lap_steady_nface   = n
      lap_steady_r_lo_Rp = r_esc
      if (n .ge. 1) then
         a(1:n) = FH_f(i1:i2)
         call median_and_spread(a(1:n), lap_steady_FH_median,             &
                                        lap_steady_FH_spread)
         a(1:n) = FHe_f(i1:i2)
         call median_and_spread(a(1:n), lap_steady_FHe_median,            &
                                        lap_steady_FHe_spread)
         a(1:n) = M_f(i1:i2)
         call median_and_spread(a(1:n), lap_steady_Mdot_median,           &
                                        lap_steady_Mdot_spread)
      endif

      end subroutine reduce_element_flux_windows

      ! ------------------------------------------------------------------ !

      logical function element_row_terms_on()
      ! EXHALE_ELEMENT_ROW_TERMS=1 writes the helium element row of every
      ! state the certification measures, term by term, to
      ! output/element_row_terms.txt.  Default off; nothing in the solution
      ! reads the record, so the run it is switched on in is the run it is
      ! switched off in.
      character(len=8) :: envv
      call get_environment_variable('EXHALE_ELEMENT_ROW_TERMS', envv)
      element_row_terms_on = (trim(envv) .eq. '1')
      end function element_row_terms_on

      ! ------------------------------------------------------------------ !

      subroutine element_row_terms_write(fname, terms)
      ! THE HELIUM ELEMENT ROW OF ONE STATE, TERM BY TERM: one line for each
      ! cell, the carrier record's analogue (carrier_row_terms_write).  Every
      ! term is the one element_transport_residual formed the row from, so
      !
      !   diffusive = dif_in + dif_out,   advective = adv_in + adv_out,
      !   residual  = diffusive + advective
      !
      ! up to the rounding of the regrouping (the time term of the row is
      ! identically zero at the stationary step length), and measure =
      ! |residual|/max(scale, element_row_scale_floor) is, bit for bit, the
      ! measure certification_row_measure takes of the same row.  The face
      ! columns are those of the inner face j-1 and the outer face j of the
      ! cell; cell 1 is the Dirichlet reservoir and carries no equation.
      character(len=*),               intent(in) :: fname
      type(element_row_terms_record), intent(in) :: terms
      integer :: u, j
      real*8  :: rin, rout, dif, advt
      if (.not. allocated(terms%res)) return
      open(newunit=u, file=trim(fname), status='replace', action='write')
      write(u,'(A)') '# helium element row terms of the state the '//     &
                     'certification measured (EXHALE_ELEMENT_ROW_TERMS=1)'
      write(u,'(A)') '# every rate is a mass rate density [g cm^-3 s^-1]'// &
                     ', signed as it enters the row; the row is the '//   &
                     'stationary balance transport = 0'
      write(u,'(A)') '# dif_in/dif_out: diffusive (gradient, eddy and '// &
                     'settling drift) flux of the inner/outer face, with'// &
                     ' its area, over the cell volume; adv_in/adv_out: '// &
                     'the same for F_rho X_face'
      write(u,'(A)') '# diffusive = dif_in + dif_out, advective = '//     &
                     'adv_in + adv_out, residual = diffusive + '//        &
                     'advective (to rounding); scale = the sum of the '// &
                     'magnitudes of the row''s terms, floored at floor'
      write(u,'(A)') '# measure = |residual|/max(scale, 1e-300), the '//  &
                     'measure the certification takes of this row; '//    &
                     'F_rho = Riemann face mass flux, J = diffusive '//   &
                     'helium flux [g/cm2/s]'
      write(u,'(A)') '# columns: cell r[R_p] r_in[R_p] r_out[R_p] '//     &
                     'rho[g/cm3] X_He rho_He[g/cm3] D_12[cm2/s] '//       &
                     'K_zz[cm2/s] dif_in dif_out diffusive adv_in '//     &
                     'adv_out advective residual floor scale measure '//  &
                     'F_rho_in[g/cm2/s] F_rho_out[g/cm2/s] X_face_in '//  &
                     'X_face_out J_in[g/cm2/s] J_out[g/cm2/s]'
      do j = 1, N
         rin  = r_edg(j-1)
         rout = r_edg(j)
         dif  = terms%dif_in(j) + terms%dif_out(j)
         advt = terms%adv_in(j) + terms%adv_out(j)
         write(u,'(I6,26(1X,ES23.15E3))') j, r(j), rin, rout,             &
              terms%rho_phys(j), terms%X(j),                              &
              terms%rho_phys(j)*terms%X(j), terms%D12(j), terms%Kzz(j),   &
              terms%dif_in(j), terms%dif_out(j), dif,                     &
              terms%adv_in(j), terms%adv_out(j), advt,                    &
              terms%res(j), terms%floor(j), terms%scale(j),               &
              abs(terms%res(j))/max(terms%scale(j),                       &
                                    element_row_scale_floor),             &
              terms%F_rho(j-1), terms%F_rho(j),                           &
              terms%Y_face(j-1), terms%Y_face(j),                         &
              terms%J(j-1), terms%J(j)
      enddo
      close(u)
      end subroutine element_row_terms_write

      ! ------------------------------------------------------------------ !

      subroutine median_and_spread(a, med, sprd)
      ! Median and relative radial spread (max - min)/|median| of one window.
      real*8, dimension(:), intent(in)  :: a
      real*8,               intent(out) :: med, sprd
      med  = median_of(a)
      sprd = (maxval(a) - minval(a))/max(abs(med), 1.0d-99)
      end subroutine median_and_spread

      ! ------------------------------------------------------------------ !

      subroutine report_step(Xhe, rho_phys, mass_resid,                    &
                             Xover, Xunder, qdep, n_vanished)
      ! Diagnostic written each step (EXHALE_DIFFUSION_CHECK=1): the range of X, the
      ! largest |J_He + J_1| over the faces, the largest relative change of the
      ! mixture's own mass sum across the step (mixture_mass_split, the mass
      ! the operator is not allowed to move), the total helium mass in the
      ! domain, how far the solve left
      ! [0,1] before the clamp, and the metal census.
      !
      ! X over / under are the excursions the clamp removed, signed so that a
      ! negative number means the solve stayed inside the range.  They are the
      ! measurement behind the claim that the clamp is an assertion: the
      ! clamped X cannot distinguish a solve that overshot from one that did
      ! not, because both print 1.0.
      !
      ! metal/H dep is max |(n_X/n_H)/melem_ab - 1| over the elements the
      ! reservoir states, and "gone" counts the (cell, element) pairs in which
      ! an element has vanished while hydrogen is present -- a state these
      ! equations cannot reach, since they have no sink for a metal nucleus.
      !
      ! A binary mixture has ONE independent diffusive flux: the operator
      ! computes J_He and the hydrogen component carries -J_He, so J_He + J_1
      ! is zero by construction and there is nothing there to measure -- it
      ! used to be printed as max|Jf - Jf|, which is an identical zero
      ! wearing the clothes of a measurement, and is gone.  The quantity that
      ! actually guards the discretization is the mass closure: it is what keeps rho -- owned by the hydro -- consistent with
      ! the composition the operator hands back, and it is nonzero the moment
      ! the write-back stops meeting both element totals exactly.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: Xhe, rho_phys
      real*8,                       intent(in) :: mass_resid
      real*8,                       intent(in) :: Xover, Xunder, qdep
      integer,                      intent(in) :: n_vanished
      real*8, dimension(0:N) :: fa
      real*8, dimension(1:N) :: cv
      real*8  :: mHe_tot
      integer :: j

      ! The helium mass of the column, on the same exact shell volumes the
      ! transport rows divide by.
      call spherical_face_area_and_cell_volume(fa, cv)
      mHe_tot = 0.0d0
      do j = 1, N
         mHe_tot = mHe_tot + rho_phys(j)*Xhe(j)*cv(j)*R0**3
      enddo
      write(0,'(A,ES13.6,A,ES13.6,A,ES9.2,A,ES16.9)')                     &
         ' (diffusion) X min ', minval(Xhe(1:N)), ' max ',                &
         maxval(Xhe(1:N)),                                                &
         ' mass closure ', mass_resid, ' He mass ', 4.0d0*pi*mHe_tot
      write(0,'(A,ES10.3,A,ES10.3,A,ES10.3,A,I0)')                        &
         ' (diffusion) X over 1 by ', Xover, ' under 0 by ', Xunder,      &
         ' metal/H dep ', qdep, ' gone ', n_vanished
      write(0,'(A,I0,A,ES10.3,A,ES10.3)')                                 &
         ' (diffusion) Newton steps ', he_fraction_newton_steps,          &
         ' resid ', he_fraction_newton_resid,                             &
         ' trace fX under 0 by ', trace_ratio_under_zero

      end subroutine report_step

      ! End of module
      end module binary_element_diffusion

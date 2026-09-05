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
      ! WHETHER THE BASE FACE IS AN INFLOW is decided by the wind's mass
      ! flux and not by the base cell's velocity, for the reason in sec. 6.
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
      ! ADVECTION IS SECOND ORDER, BY A LIMITED SLOPE ON TOP OF THE ONE-SIDED
      ! UPWIND DIFFERENCE, AND THE UPWIND PART IS ON THE CELL VELOCITY IN
      ! THE INTERIOR, not on face-averaged velocities, and that is a lesson
      ! inherited rather than rediscovered: binary_element_diffusion's header
      ! records that the face-averaged form lets a cell whose two face
      ! velocities straddle zero lose its advective term entirely, which
      ! froze an isolated helium hole at the breathing HD 209458 b base.
      !
      ! THE BASE CELL IS THE ONE EXCEPTION, AND IT IS AN EXCEPTION BECAUSE
      ! THE CELL VELOCITY THERE IS NOT A FLUX THE SCHEME TRANSPORTS.
      ! Measured on the converged hot Uranus (docs/p44_base_sawtooth.md): the
      ! cell-centred rho v r^2 of cell 1 is 160 to 200 times the wind's flux
      ! INWARD, while the face flux reconstructed from the scheme's own mass
      ! residual is 0.36 to 0.98 times the wind's flux OUTWARD, and the
      ! discrepancy decays over four cells. That is the signature of a
      ! collocated odd-even velocity mode, not of infall, and a rule that
      ! reads sign(v(1)) reads the mode. On the oxygen example the same
      ! quantity is 1.5e6 times the wind flux inward.
      !
      ! So cell 1's advective direction, and the base composition condition
      ! of sec. 5, are both decided by the MASS FLUX THE WIND CARRIES --
      ! sum(rho v r^2) over the escape region [j_min:N], the same constant
      ! and the same average "Base velocity: massflux" uses for the ghost
      ! velocity (init.f90, EXHALE_main.f90). In a steady spherical wind that
      ! flux is one number at every face, so using it at the base face is not
      ! an approximation of the steady state; it is a statement of it.
      !
      ! WHY SECOND ORDER, AND WHAT IT COST.  First order donor-cell carries a
      ! numerical diffusivity |v| dr/2.  Measured on the 500-cell grid of the
      ! hot Uranus that is 8.5, 20 and 22 times the molecular D of H2 at 1.15,
      ! 1.30 and 1.50 R_p, and a controlled grid test says it is not a
      ! bookkeeping detail: refining the H2 front by a factor two at UNCHANGED
      ! CFL step (500 -> 1000 cells, the base spacing untouched, so the global
      ! dt moves by 2%) cut the front's advance rate by a factor five and
      ! removed its re-acceleration entirely.  The front this operator
      ! produces was being carried by the scheme as much as by the flow.
      !
      ! The second-order term is a DEFERRED CORRECTION: the residual carries
      !   d f/dr |_j  =  (f_j - f_{j-1})/h  +  (sigma_j - sigma_{j-1})/2
      ! for v >= 0 (and the mirror image for v < 0), with sigma the limited
      ! slope.  On a linear profile the two sigma cancel and on a quadratic
      ! they cancel the donor-cell truncation term exactly, so the scheme is
      ! second order; with sigma = 0 it is EXACTLY the first-order operator
      ! this module had, on the same non-uniform spacing, so the correction is
      ! attributable on its own.
      !
      ! sigma is FROZEN over one Newton solve and refreshed at every transport
      ! step.  That is what keeps the Jacobian block-tridiagonal -- the
      ! correction reaches j-2 and j+2, which block_thomas cannot carry -- and
      ! it is also what keeps the residual smooth, since a limiter is not
      ! differentiable where its branches meet.  At the fixed point of the
      ! marching loop (and of the fixed-wind relaxation) consecutive steps
      ! agree, so the state the code converges to is the second-order one.
      !
      ! THE LIMITER IS VAN LEER, not minmod, for two reasons that are about
      ! this solver rather than about accuracy in general: van Leer's
      ! 2ab/(a+b) is differentiable everywhere except at the sign change,
      ! while minmod has a kink wherever |a| = |b| -- and a kink in the source
      ! of a Newton residual is the defect section 126 spent a campaign on --
      ! and van Leer is the less diffusive of the two, which is the property
      ! the measurement above says is scarce here.
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
      use caloric_eos, only: adiabatic_index_at_T
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
      use steady_residual_mod, only: carrier_row_scale
      use composition, only: base_h2_composition_imposed
      use element_census, only: element_census_state, element_census_take,&
                               element_census_verify

      implicit none
      private

      public :: photochemical_transport_step
      public :: relax_photochemical_composition
      public :: carrier_transport_diagnostics, carrier_co_ceiling_cells
      public :: carrier_diffusion_coefficient
      public :: carrier_set_init
      public :: carrier_slope, n_carrier_max   ! order-of-accuracy unit test
      ! Element-closure unit test (src/tests/element_census_tests.f90): the
      ! write-back's invariant is that the element totals it is handed come
      ! back unchanged, and the test states it on the two routines that make
      ! the pair.
      public :: carrier_state, carrier_write_back
      public :: carrier_steady_residual
      public :: carrier_drift_location
      ! The coupled steady solve (sec. 139): n(H2) is a Newton unknown, so
      ! the solver needs the carrier row, the scale of its unknown, and the
      ! element headroom it must stay inside.
      public :: carrier_column_scale, carrier_row_scale_H2,             &
                carrier_headroom
      public :: n_carrier, carrier_name, ic_H2, ic_OH, ic_H2O, ic_CO
      public :: ic_Hp, carrier_solved

      integer, parameter :: dp = kind(1.0d0)

      ! The transported set, in the order the solve uses.
      !
      ! WHICH SPECIES ARE IN IT IS A PROPERTY OF THE RUN, NOT OF THE MODULE.
      ! H2 is carrier 1 in every configuration -- it is the one carrier
      ! whose Damkohler number is of order unity in the region the option
      ! exists for -- and the oxygen cycle adds OH, H2O and CO behind it.
      ! Measured on the converged He/H = 0.0793 hot Uranus, the H2 chemical
      ! time is 0.02 to 1.3 times the flow time r/|v| across the front,
      ! while the three molecular ions that are left local sit at 1e-9 to
      ! 5e-5 of it and hold at most 2e-8 of the H nuclei, so the split
      ! between carried and local is a measurement and not a convention
      ! (docs/supersonic_molecular_base.md sec. 13.6).
      !
      ! n_carrier_max dimensions every array; n_carrier is how many of them
      ! this run solves, set once by carrier_set_init.  The inactive columns
      ! stay at zero -- the species behind them do not exist without the
      ! oxygen chemistry -- so the whole-array statements below are exact.
      !
      ! THE PROTON IS THE FIFTH, AND IT IS LAST FOR A REASON. Its index has
      ! to be a constant -- the ionization rows, the write-back and the
      ! element budget all name it -- and the four existing indices have to
      ! keep the values they have, so the only place left is behind them.
      ! That makes the solved set non-contiguous in the one configuration
      ! that carries H2 and H+ without the oxygen cycle, which is why
      ! carrier_solved exists below: n_carrier is the LAST index solved and
      ! carrier_solved says which of 1..n_carrier actually are.
      integer, parameter :: n_carrier_max = 5
      integer, save      :: n_carrier = 4   ! set by carrier_set_init
      integer, parameter :: ic_H2 = 1, ic_OH = 2, ic_H2O = 3, ic_CO = 4
      integer, parameter :: ic_Hp = 5
      character(len=3), parameter :: carrier_name(n_carrier_max) =       &
           (/ 'H2 ', 'OH ', 'H2O', 'CO ', 'H+ ' /)
      ! Which of 1..n_carrier this run solves. A carrier that is not solved
      ! keeps the identity row nrho/dt on the diagonal and a zero residual,
      ! so its unknown does not move and the block structure the Thomas
      ! sweep needs is unchanged; it costs one 5x5 block instead of one 2x2
      ! and no chemistry evaluation. Set once by carrier_set_init.
      logical, save :: carrier_solved(n_carrier_max) = .true.

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

      ! THE BASE FACE OF AN INFLOWING CARRIER (design decision D6, amended).
      !
      ! D6 gave the carrier partition a ZERO-FLUX base: the element totals
      ! are boundary data the handoff supplies, but the split of an element
      ! among its carriers is what the oxygen option exists to COMPUTE, so
      ! imposing it at the boundary would have made that option's own A/B
      ! gate circular.  That argument holds exactly when nothing upstream
      ! states the partition.
      !
      ! When a lower-atmosphere handoff DOES state it -- base.inp's
      ! q_H2_base, or the profile's value at the matching level -- the
      ! partition is not a quantity this operator computes but incoming
      ! data, and section 117 already pins it on the ghosts of the
      ! ionization sweep.  Leaving the base cell free while the ghost below
      ! it is pinned puts the same defect section 117 removed one cell
      ! higher up: the equation of state and the species state would again
      ! describe different gas across one face.  So where the handoff speaks,
      ! the advective inflow at the base face carries ITS composition.
      !
      ! The two branches cannot both apply.  Decision D7 refuses q_H2_base
      ! while the oxygen chemistry is on, so a run either has a handoff
      ! partition (Dirichlet here, H2 only -- nothing states OH, H2O or CO)
      ! or computes its own (zero flux, as D6 wrote it).
      !
      ! ONLY THE ADVECTIVE TERM.  Face 0 still carries no diffusive flux.
      ! The handoff supplies the composition of the gas that flows in, not a
      ! mixing rate across a boundary whose gradient is set by the ghost
      ! spacing and by a K_zz that describes an unresolved region.
      logical,  save :: base_dirichlet(n_carrier_max) = .false.
      ! Does the base FACE carry gas into the domain? Set from the wind's
      ! mass flux, never from the base cell's velocity (sec. 6).
      logical,  save :: base_inflow = .true.
      ! Deferred second-order correction to the advective derivative, in the
      ! residual's own units [1/cm], frozen over one Newton solve.
      real(dp), dimension(:,:), allocatable :: adv_corr
      real(dp), save :: fc_base(n_carrier_max) = 0.0d0

      ! Diagnostics of the last step, read by write_output.
      integer,  save :: pct_newton_steps = 0
      real(dp), save :: pct_newton_resid = 0.0d0
      integer,  save :: pct_limited      = 0
      integer,  save :: pct_co_ceiling    = 0
      integer,  save :: pct_worst_j       = 0
      integer,  save :: pct_worst_ic      = 1
      ! Where the last fixed-wind relaxation moved the composition most.
      integer,  save :: pct_drift_j     = 0
      integer,  save :: pct_drift_ic    = 1
      real(dp), save :: pct_worst_limit  = 0.0d0
      ! Molecular diffusion coefficient of each carrier on the grid, kept
      ! from the last step so the run can print the transport time scale
      ! beside the chemical one it already prints.
      real(dp), dimension(:,:), allocatable :: pct_Dco
      ! The sum of each carrier row's own term magnitudes, from the last
      ! residual assembly.  This is the row's scale (sec. 139): the residual
      ! measure, the solver's column scale and the acceptance gate all read
      ! it from here, so there is one definition.
      real(dp), dimension(:,:), allocatable :: row_terms
      ! Scale of the carrier ROW, and separately of its UNKNOWN, from that
      ! same assembly.  The ROW scale is row_terms over one code time unit,
      ! so that the residual divided by it is the dimensionless imbalance
      ! the acceptance test uses.  The COLUMN scale is the unknown's own
      ! magnitude with a floor that does not collapse where it does -- the
      ! same construction as the momentum unknown's |rho v| + rho c_s.
      !
      ! THE TWO MUST DIFFER, and that was measured.  n(H2) alone cannot be
      ! either: it falls four decades across the front and is effectively
      ! zero beyond it, which is the pathology carrier_row_scale was written
      ! to avoid.  Using the ROW scale for both makes the row's merit
      ! contribution the CHEMICAL rate, about 1e3 per code time at the front
      ! against the hydrodynamic rows' 1e-6, so the line search sees nothing
      ! but the carrier and the damped Gauss-Newton escape lands on a state
      ! whose chemistry cannot be certified (merit 89.8, abort at iteration
      ! 2).  Using the COLUMN scale for both makes the Newton propose
      ! carrier steps 84 times the cell's own n(H2) and cuts the line search
      ! to lam = 4.9e-4.
      real(dp), dimension(:), allocatable :: col_scale_H2, row_scale_H2
      ! Element headroom of each cell from that same assembly [cm^-3]: the
      ! hydrogen the carriers may hold, halved because H2 holds two nuclei.
      real(dp), dimension(:), allocatable :: headroom_H2
      logical :: headroom_set = .false.

      contains

      ! Fix the transported set for this run.  Called once from input_read,
      ! after every key is parsed, so that no routine can see a half-decided
      ! carrier set.
      subroutine carrier_set_init()
      carrier_solved = .false.
      carrier_solved(ic_H2) = .true.
      if (thereis_oxychem) then
         carrier_solved(ic_OH)  = .true.
         carrier_solved(ic_H2O) = .true.
         carrier_solved(ic_CO)  = .true.
      endif
      if (ionization_transport) carrier_solved(ic_Hp) = .true.
      ! n_carrier is the last index in the solved set: the loops run
      ! 1..n_carrier and skip the ones carrier_solved leaves out.
      if (ionization_transport) then
         n_carrier = ic_Hp              ! ... H2 [, OH, H2O, CO] and H+
      else if (thereis_oxychem) then
         n_carrier = ic_CO              ! H2, OH, H2O, CO
      else
         n_carrier = ic_H2              ! H2 alone
      endif
      end subroutine carrier_set_init

      ! ------------------------------------------------------------- !

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

      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fc, fc_old
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: Dco
      real(dp), dimension(1-Ng:N+Ng)           :: ntot, TK, gphys, mbar
      real(dp), dimension(1-Ng:N+Ng)           :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng)           :: rp, rep, ntv, dt_phys
      real(dp), dimension(1-Ng:N+Ng)           :: sigrate
      real(dp), dimension(1-Ng:N+Ng)           :: nH_free, nO_free, nC_free
      real(dp), dimension(0:N,n_carrier_max)       :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max)       :: updrf
      ! A1 element-budget assertions.  The transport step holds rho fixed and
      ! its write-back restores the ENTRY element totals cell by cell (that is
      ! what nH_free / nO_free / nC_free are), so the absolute nucleus density
      ! of every element is invariant across both the write-back alone and the
      ! whole step.  Off unless EXHALE_ELEMENT_ASSERT is set.
      type(element_census_state) :: cen_step, cen_wb

      if (.not. thereis_mol)       return
      if (.not. carrier_transport) return
      if (.not. bg_ready)          return

      call element_census_take('photochemical_transport_step',            &
                               rho, f_sp, cen_step)
      call carrier_state(rho, Tcode, f_sp, fc, ntot, nrho, wfac, TK,     &
                         mbar, nH_free, nO_free, nC_free)
      call carrier_base_state(fc, rho, v)
      fc_old = fc

      ! The advective coefficient is built on nrho, not on n_tot: the row
      ! transports the fraction per unit mass (sec. 158).
      call carrier_geometry(v, dt_code, nrho, rp, rep, ntv, dt_phys,     &
                            gphys)
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      if (.not. allocated(pct_Dco)) allocate(pct_Dco(1-Ng:N+Ng,n_carrier_max))
      pct_Dco = Dco
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,     &
                                     Agrd, Bdrf, updrf)
      call carrier_photolysis(rho, TK, f_sp)
      call carrier_advection_correction(fc, rp, ntv)

      call carrier_signal_rate(v, TK, mbar, sigrate)
      call solve_carriers(fc, fc_old, nrho, wfac, dt_phys, rp, rep, ntv, &
                          Agrd, Bdrf, updrf, TK, sigrate,                &
                          nH_free, nO_free, nC_free)

      call element_census_take('carrier_write_back', rho, f_sp, cen_wb)
      call carrier_write_back(rho, f_sp, fc, nrho, nH_free, nO_free,     &
                              nC_free)
      call element_census_verify(cen_wb,   rho, f_sp, rho_is_fixed=.true.)
      call element_census_verify(cen_step, rho, f_sp, rho_is_fixed=.true.)

      if (carrier_debug_on()) then
         write(*,'(a,i4,a,es9.2,a,f6.0,a,i0,a,es9.2,a,f6.0,a,es9.2,'//  &
              'a,es9.2,a,es9.2)')                                        &
            ' (carrier_transport) newton ', pct_newton_steps,            &
            ' resid ', pct_newton_resid, ' at cell '//                  &
            trim(carrier_name(pct_worst_ic))//' j=', dble(pct_worst_j),  &
            ' limited ', pct_limited,                                    &
            ' cells, worst ', pct_worst_limit, ', CO ceiling ',          &
            dble(pct_co_ceiling), '; base x_H2 ',                        &
            2.0d0*fc(1,ic_H2)*nrho(1)/max(nH_free(1), 1.0d-99),          &
            ' <- ', 2.0d0*fc_old(1,ic_H2)*nrho(1)                        &
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
      subroutine carrier_state(rho, Tcode, f_sp, fc, ntot, nrho, wfac,   &
                               TK, mbar, nH_free, nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rho, Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(out) :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: ntot, TK
      ! nrho = rho*n0, the density the carrier fraction is a fraction OF --
      ! the one the hydrodynamic mass row conserves exactly.  wfac converts
      ! that fraction into the PARTICLE mixing ratio the diffusive flux is
      ! driven by: X = wfac*y, wfac = rho*n0/n_tot = mbar/m_amu.
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nrho, wfac
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

      ! THE CARRIER UNKNOWN IS THE FRACTION PER UNIT MASS, not per particle.
      ! Section 158: the marching loop's hydrodynamic stage moves rho at
      ! frozen f_sp, so the carrier density it hands on changes by
      ! n_c d(ln rho)/dt.  A carrier advected as n_c/n_tot instead composes
      ! with that stage into the conservation law PLUS a spurious
      ! n_c v d(ln mbar)/dr, which lives exactly across the dissociation and
      ! ionization fronts where mbar varies -- measured at 13 times the row's
      ! own steady residual in the wind.  n_c/(rho n0) is the variable the
      ! mass row conserves, and it is the variable f_sp already stores.
      nrho = nd
      wfac = nd/max(ntot, 1.0d-99)
      fc(:,ic_H2)  = f_sp(:,isp_H2)
      fc(:,ic_OH)  = f_sp(:,isp_OH)
      fc(:,ic_H2O) = f_sp(:,isp_H2O)
      fc(:,ic_CO)  = f_sp(:,isp_CO)
      ! The proton column is filled whether or not it is solved -- it costs
      ! one assignment and it keeps fc a complete picture of the state -- but
      ! only a run with ionization_transport lets the solve move it.
      fc(:,ic_Hp)  = f_sp(:,isp_HII)
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
      ! THAT ARGUMENT DOES NOT REACH AS FAR AS IT WAS TAKEN TO, and hydrogen
      ! is where it breaks.  "The next sweep re-solves the stages" is true of
      ! the sweep; it is not true of the step in between.  This step's
      ! write-back does not re-solve H I, H II, H2+ and H3+ -- it rescales
      ! them in proportion to what they already held -- so nuclei the
      ! carriers take are nuclei those stages simply lose, and a carrier
      ! holding the whole element leaves the sweep a cell with no electrons.
      ! So the ELEMENT total is the right headroom for the LIMITER's outer
      ! question (which element, and how much of it exists), and
      ! hydrogen_available_to_carriers is the right one for what the
      ! carriers may actually take.  nH_free below is the first; the limiter
      ! and the chemistry rows both go through the second.
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

      ! Cell and carrier that the last fixed-wind relaxation moved most.
      subroutine carrier_drift_location(j, ic)
      integer, intent(out) :: j, ic
      j  = pct_drift_j
      ic = pct_drift_ic
      end subroutine carrier_drift_location

      ! ------------------------------------------------------------- !

      ! Is the state (rho, v, T, f_sp) a steady state OF THE CARRIER
      ! EQUATION, and where is it least so?
      !
      ! The steady solver never sees this equation: it holds the carrier
      ! partition fixed and solves the wind, and an outer Picard iteration
      ! relaxes the carriers at that wind (EXHALE_main). Their common fixed
      ! point is a steady state of both, but ONLY the hydrodynamic half of it
      ! is measured by ||R||. This is the other half.
      !
      ! It is the residual the transport step already assembles, evaluated
      ! with the time term switched off (dt -> infinity with fc_old = fc
      ! makes it identically zero), divided by carrier_row_scale -- the same
      ! place, and the same principle, as every hydrodynamic row scale
      ! (steady_residual.f90). A no-op, reported as zero, when the carriers
      ! are not transported.
      subroutine carrier_steady_residual(rho, v, Tcode, f_sp, rcmax,     &
                                         jworst, icworst, rvol, rlegacy, &
                                         res_out, terms_out)
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho, v
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real(dp), intent(out) :: rcmax
      integer,  intent(out) :: jworst, icworst
      ! The volume-weighted measure, the one a gate should read: the larger
      ! of the layer's and the wind's ratio of sums, exactly as
      ! relnorm_over_cells forms the hydrodynamic rows.
      real(dp), intent(out), optional :: rvol
      ! The same residual on the retired n_H(|v|+c_s)/dr scale, so that a run
      ! can report both while the two are being compared.
      real(dp), intent(out), optional :: rlegacy
      real(dp), dimension(1:N,n_carrier_max), intent(out), optional ::   &
                                                 res_out, terms_out

      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fc, Dco
      real(dp), dimension(1-Ng:N+Ng) :: ntot, TK, mbar, gphys
      real(dp), dimension(1-Ng:N+Ng) :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng) :: rp, rep, ntv, dt_big, sigrate
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp), dimension(0:N,n_carrier_max) :: Agrd, Bdrf, Jf, dJl, dJr
      integer,  dimension(0:N,n_carrier_max) :: updrf
      real(dp), dimension(1:N,n_carrier_max) :: res
      real(dp) :: rnorm_unused, cs, sc, rel, w
      real(dp) :: num_l, den_l, num_w, den_w
      integer  :: j, ic

      rcmax = 0.0d0;  jworst = 0;  icworst = 1
      if (present(rvol))      rvol      = 0.0d0
      if (present(rlegacy))   rlegacy   = 0.0d0
      if (present(res_out))   res_out   = 0.0d0
      if (present(terms_out)) terms_out = 0.0d0
      if (.not. thereis_mol)       return
      if (.not. carrier_transport) return
      if (.not. bg_ready)          return

      if (.not. allocated(row_terms)) allocate(row_terms(1:N,n_carrier_max))
      row_terms = 0.0d0
      call carrier_state(rho, Tcode, f_sp, fc, ntot, nrho, wfac, TK,     &
                         mbar, nH_free, nO_free, nC_free)
      call carrier_base_state(fc, rho, v)
      call carrier_geometry(v, spread(1.0d0,1,N+2*Ng), nrho, rp, rep,    &
                            ntv, dt_big, gphys)
      dt_big = 1.0d30                 ! kills the time term exactly
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      call carrier_face_coefficients(ntot, TK, mbar, gphys, Dco, rp,     &
                                     Agrd, Bdrf, updrf)
      call carrier_photolysis(rho, TK, f_sp)
      call carrier_advection_correction(fc, rp, ntv)
      call carrier_signal_rate(v, TK, mbar, sigrate)
      call carrier_residual(fc, fc, nrho, wfac, dt_big, rp, rep, ntv,    &
                            Agrd, Bdrf, updrf, nH_free, nO_free,         &
                            nC_free, sigrate, res, Jf, dJl, dJr,         &
                            rnorm_unused)

      ! THE SCALE IS THE ROW'S OWN LARGEST TERMS, and that is a change.
      ! It used to be carrier_row_scale, n_H(|v|+c_s)/dr -- the flux
      ! divergence the cell would carry if its WHOLE hydrogen inventory moved
      ! at the signal speed. That is a bound on the row, but it is not the
      ! row's own term, and measured on the He/H = 0.0793 hot Uranus it
      ! exceeds the sum of the terms the row actually balances by 10^2.7 to
      ! 10^6.4 (docs/p50_carrier_wind_alternation.md section 2.4). A state
      ! the steady solver accepted at info = 0 then reported a carrier
      ! residual of 1.6e-4 while carrying 0.60 of its own terms at the H2
      ! front, and the outer loop had no way to see that the state it was
      ! handing over was not a carrier steady state at all.
      !
      ! Section 133 already made "a bound on that row's own largest term" the
      ! principle for every hydrodynamic row. This is the carrier row on the
      ! same principle, and the quantity it needs is the one
      ! carrier_residual has always computed for its own Newton: row_terms.
      if (.not. allocated(col_scale_H2))                                 &
         allocate(col_scale_H2(1:N), row_scale_H2(1:N), headroom_H2(1:N))
      num_l = 0.0d0;  den_l = 0.0d0;  num_w = 0.0d0;  den_w = 0.0d0
      do j = 1, N
         cs = sqrt(adiabatic_index_at_T(j, TK(j)/T0)*kb_erg                &
                  *max(TK(j), 1.0d0)/max(mbar(j), 1.0d-40))
         sc = carrier_row_scale(j, nH_free(j), v(j)*v0, cs)
         w  = r(j)*r(j)*dr_j(j)
         row_scale_H2(j) = row_terms(j,ic_H2)*R0/v0
         col_scale_H2(j) = fc(j,ic_H2)*nrho(j)                           &
                         + row_terms(j,ic_H2)/max(sigrate(j), 1.0d-300)
         ! THE ELEMENT HEADROOM IS FROZEN AT THE OUTER ITERATE, and the
         ! reference matters. It is a constraint the solver tests trial
         ! states against, and it depends on the very unknown it constrains
         ! -- n_H_free counts H2's own two nuclei -- so evaluating it at the
         ! trial makes it move with the step. Refreshing it at every state
         ! that is not a Jacobian probe was tried and is worse than that:
         ! the last LINE-SEARCH trial then leaves its own headroom behind,
         ! and the Gauss-Newton escape that follows finds the ITERATE
         ! outside it, so every backtracked step is refused down to
         ! lam = 5e-4 and the solve aborts holding a descent direction
         ! (measured). weno_mode = 1 is the one evaluation per outer
         ! iteration that is the iterate itself.
         if (weno_mode .eq. 1 .or. .not. headroom_set)                   &
            headroom_H2(j) =                                             &
               0.5d0*hydrogen_available_to_carriers(j, nH_free(j))
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            rel = abs(res(j,ic))/max(row_terms(j,ic), 1.0d-300)
            if (rel .gt. rcmax) then
               rcmax = rel;  jworst = j;  icworst = ic
            endif
            if (present(rlegacy)) rlegacy = max(rlegacy,                 &
                                     abs(res(j,ic))/max(sc, 1.0d-300))
            if (j .lt. j_min) then
               num_l = num_l + abs(res(j,ic))*w
               den_l = den_l + row_terms(j,ic)*w
            else
               num_w = num_w + abs(res(j,ic))*w
               den_w = den_w + row_terms(j,ic)*w
            endif
         enddo
      enddo
      ! The two regions are normed separately and combined by the larger,
      ! for the reason residual_norms gives: the wind carries almost all of
      ! sum r^2 dr, so any single sum lets it average the inner column away.
      if (present(rvol)) rvol = max(num_l/max(den_l, 1.0d-300),          &
                                    num_w/max(den_w, 1.0d-300))
      headroom_set = .true.
      if (present(res_out))   res_out   = res(1:N,:)
      if (present(terms_out)) terms_out = row_terms(1:N,:)
      end subroutine carrier_steady_residual

      ! ------------------------------------------------------------- !

      ! COLUMN scale of the carrier unknown in cell j [cm^-3]: n(H2) itself
      ! plus the amount of it a signal crossing would carry at the size of
      ! the row's own terms.  Exactly the construction of the momentum
      ! unknown's |rho v| + rho c_s -- the unknown where it means something,
      ! a floor that does not collapse where it does not.  It sets the
      ! Newton step and the finite-difference step, so it has to be the size
      ! of the variable; using the ROW scale for it instead was measured to
      ! propose carrier steps 84 times the cell's own n(H2) and to send the
      ! line search to lam = 4.9e-4.  Zero before the first assembly, which
      ! the caller reads as "no scale yet".
      double precision function carrier_column_scale(j) result(D)
      integer, intent(in) :: j
      D = 0.0d0
      if (allocated(col_scale_H2)) D = col_scale_H2(j)
      end function carrier_column_scale

      ! ROW scale of the carrier row in cell j [cm^-3]: the row's own
      ! largest terms accumulated over one code time unit, so that the
      ! residual divided by it is the dimensionless imbalance the acceptance
      ! test uses.  It is NOT the column scale, and the two must differ:
      ! the row is stiff -- its natural rate is the chemical one, 1e3 per
      ! code time at the front -- while its unknown is small, so one vector
      ! cannot scale both without either drowning the merit or overshooting
      ! the step.
      double precision function carrier_row_scale_H2(j) result(D)
      integer, intent(in) :: j
      D = 0.0d0
      if (allocated(row_scale_H2)) D = row_scale_H2(j)
      end function carrier_row_scale_H2

      ! The largest n(H2) cell j may hold [cm^-3]: half the hydrogen the
      ! carriers are allowed, because H2 holds two nuclei.  A trial state
      ! above it is not a state the element budget admits, and the coupled
      ! solve refuses it rather than clamping it -- clamping would put a
      ! kink in the residual, which is the defect section 138 removed, and
      ! it is what made the fixed-wind relaxation return a state pinned on
      ! the ceiling in a quarter of the grid.
      double precision function carrier_headroom(j) result(nmax)
      integer, intent(in) :: j
      nmax = huge(1.0d0)
      if (allocated(headroom_H2)) nmax = headroom_H2(j)
      end function carrier_headroom

      ! ------------------------------------------------------------- !

      ! The cell's signal-crossing rate (|v| + c_s)/dr [1/s], from the same
      ! background the carrier rows are assembled on.
      subroutine carrier_signal_rate(v, TK, mbar, sigrate)
      real(dp), dimension(1-Ng:N+Ng), intent(in)  :: v, TK, mbar
      real(dp), dimension(1-Ng:N+Ng), intent(out) :: sigrate
      real(dp) :: cs
      integer  :: j
      do j = 1-Ng, N+Ng
         cs = sqrt(adiabatic_index_at_T(j, TK(j)/T0)*kb_erg                &
                  *max(TK(j), 1.0d0)/max(mbar(j), 1.0d-40))
         sigrate(j) = (abs(v(j))*v0 + cs)                                &
                    / max((r_edg(min(max(j,1),N)) -                      &
                           r_edg(min(max(j,1),N)-1))*R0, 1.0d0)
      enddo
      end subroutine carrier_signal_rate

      ! ------------------------------------------------------------- !

      ! Van Leer limited slope of carrier ic in cell j, per unit length.
      ! a and b are the two one-sided differences; the harmonic mean 2ab/(a+b)
      ! is zero when they disagree in sign, which is what makes the scheme
      ! monotone across the H2 front.
      double precision function carrier_slope(fc, ic, j, rp) result(sig)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fc
      integer,                                      intent(in) :: ic, j
      real(dp), dimension(1-Ng:N+Ng),               intent(in) :: rp
      real(dp) :: a, b
      sig = 0.0d0
      if (j .lt. 2 .or. j .gt. N-1) return
      a = (fc(j,ic)   - fc(j-1,ic))/max(rp(j)   - rp(j-1), 1.0d0)
      b = (fc(j+1,ic) - fc(j,ic)  )/max(rp(j+1) - rp(j),   1.0d0)
      if (a*b .le. 0.0d0) return
      sig = 2.0d0*a*b/(a + b)
      end function carrier_slope

      ! ------------------------------------------------------------- !

      ! The deferred second-order correction of the advective derivative,
      ! built once per transport step from the state the step was handed and
      ! held fixed through the Newton (see sec. 6).  Cell 1 is left at zero:
      ! its direction comes from the base face and it has no upstream slope
      ! inside the domain.
      subroutine carrier_advection_correction(fc, rp, ntv)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fc
      real(dp), dimension(1-Ng:N+Ng),               intent(in) :: rp, ntv
      integer :: j, ic
      if (.not. allocated(adv_corr)) allocate(adv_corr(1:N,n_carrier_max))
      ! FROZEN INSIDE THE COUPLED NEWTON'S MODEL, on the same switch and for
      ! the same reason as the WENO3 weights (sec. 121): the limited slope
      ! reaches j-2 and j+2, which the banded Jacobian's stencil does not
      ! carry, and it is the strongly nonlinear part of the operator. Mode 2
      ! is the inner model -- Jacobian columns and Krylov products -- and it
      ! reuses the slopes stored at the outer iterate. Modes 0 and 1 rebuild
      ! them, so the marching path, the relaxation and the line search all
      ! see the slopes of the state they are evaluating.
      if (weno_mode .eq. 2) return
      adv_corr = 0.0d0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 2, N-1
            if (ntv(j) .ge. 0.0d0) then
               adv_corr(j,ic) = 0.5d0*(carrier_slope(fc, ic, j,   rp)     &
                                     - carrier_slope(fc, ic, j-1, rp))
            else
               adv_corr(j,ic) = -0.5d0*(carrier_slope(fc, ic, j+1, rp)    &
                                      - carrier_slope(fc, ic, j,   rp))
            endif
         enddo
      enddo
      end subroutine carrier_advection_correction

      ! ------------------------------------------------------------- !

      ! The composition the base face hands the carriers, captured before
      ! the Newton starts.  It is read from the lower ghost, which the
      ! ionization sweep has just pinned to the handoff value (sec. 117);
      ! capturing it here rather than reading fc(0,:) inside the residual is
      ! deliberate, because the Newton's own trial states overwrite the
      ! ghosts with the base cell on every line-search try.
      subroutine carrier_base_state(fc, rho, v)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in) :: fc
      real(dp), dimension(1-Ng:N+Ng),               intent(in) :: rho, v
      base_dirichlet = .false.
      fc_base        = 0.0d0
      base_inflow    = (wind_mass_flux(rho, v) .gt. 0.0d0)
      if (.not. base_inflow) return
      if (.not. base_h2_composition_imposed()) return
      base_dirichlet(ic_H2) = .true.
      fc_base(ic_H2)        = max(fc(0,ic_H2), 0.0d0)
      ! THE PROTON GETS NO DIRICHLET VALUE HERE, and the reason is the one
      ! D6 gave: a partition is imposed at the base face only where
      ! something upstream STATES it. The handoff states the molecular
      ! partition -- base.inp's q_H2_base, or the profile's value at the
      ! matching level -- because the lower atmosphere is where the
      ! dissociating field and the mixing that sets it act. It states no
      ! ionization fraction, and there is no key for one: the gas below the
      ! base is shielded and its H+ is whatever the base cell's own balance
      ! makes of it. So the proton keeps the zero-flux partition boundary of
      ! D6, i.e. the inflow carries the base cell's own composition, and the
      ! ionization sweep continues to own the ghosts. If a handoff ever
      ! states an ionization fraction this is the one place that changes.
      end subroutine carrier_base_state

      ! ------------------------------------------------------------- !

      ! The mass flux the wind carries, sum(rho v r^2) averaged over the
      ! escape region [j_min:N], in code units.  This is the SAME quantity
      ! and the SAME average that "Base velocity: massflux" uses for the
      ! ghost velocity -- seeded in init.f90 and updated each step in
      ! EXHALE_main.f90 -- so the two places that need to know what the base
      ! face carries agree by construction rather than by coincidence.  It
      ! is evaluated here from the state handed in, so it is available
      ! whether or not that option is on.
      double precision function wind_mass_flux(rho, v) result(F)
      real(dp), dimension(1-Ng:N+Ng), intent(in) :: rho, v
      if (j_min .gt. N) then
         F = 0.0d0
      else
         F = sum(rho(j_min:N)*v(j_min:N)*r(j_min:N)*r(j_min:N))          &
             /dble(N - j_min + 1)
      endif
      end function wind_mass_flux

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
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(out) :: Dco
      real(dp) :: y(n_bkg), ysum, fric, Dpair, nd, nc
      integer  :: j, ic, is
      integer  :: bkg_isp(n_bkg)
      bkg_isp = (/ isp_HI, isp_HII, isp_H2, isp_H2p, isp_H3p,            &
                   isp_HeI, isp_HeII, isp_HeIII, isp_HeHp /)
      Dco = 0.0d0
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
            if (.not. carrier_solved(ic)) cycle
            ! THE PROTON DIFFUSES AS AN ION, AND THAT IS NOT THIS OPTION.
            ! The hard-sphere pair coefficient below is the mutual diffusion
            ! of a NEUTRAL through neutrals. A proton in a mostly neutral
            ! gas is not that: it is coupled to the electrons by the
            ! ambipolar electric field, and to the other charges by Coulomb
            ! collisions whose cross section is orders of magnitude larger
            ! than the gas-kinetic one and is a function of the ionization
            ! fraction. Using the neutral coefficient for it would put a
            ! number in the equation that is not the physics of the term.
            !
            ! So the molecular diffusion of the proton is set to zero and it
            ! is carried by the bulk flow alone -- which is the term this
            ! option exists for: what Model A's Figure 15 shows carrying 43
            ! per cent of the proton production at 2.0 r_base and 56 per
            ! cent at 2.4 is the ADVECTION, not a diffusive flux. The eddy
            ! coefficient K_zz is a separate matter and is still applied to
            ! this carrier (carrier_face_coefficients): it is a bulk mixing
            ! coefficient of the gas, blind to the charge of what it mixes.
            ! Adding the ambipolar coefficient is a physics decision of its
            ! own and is not taken here.
            if (ic .eq. ic_Hp) then
               Dco(j,ic) = 0.0d0
               cycle
            endif
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
      case (ic_Hp)
         m = bsp_mass(2)
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
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: Dco
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rp
      real(dp), dimension(0:N,n_carrier_max),       intent(out) :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max),       intent(out) :: updrf
      real(dp), dimension(1-Ng:N+Ng) :: Gco
      real(dp) :: dr_f, ntf, Df, Gf, Kf
      integer  :: j, ic

      Agrd  = 0.0d0
      Bdrf  = 0.0d0
      updrf = 0

      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
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
      real(dp), dimension(1-Ng:N+Ng) :: nH2, nH2O, nOH, nH_nuc
      real(dp), dimension(1-Ng:N+Ng) :: c1, c2, c3, fsh, trl, klw
      real(dp), dimension(1-Ng:N+Ng) :: pls, pla
      real(dp), dimension(1-Ng:N+Ng,n_fuv_band) :: tau_b, jh2o, joh
      real(dp) :: amax, ovl
      integer  :: j
      nH2  = rho*n0*f_sp(:,isp_H2)
      nH2O = rho*n0*f_sp(:,isp_H2O)
      nOH  = rho*n0*f_sp(:,isp_OH)
      ! Total hydrogen NUCLEUS density, the density axis of the
      ! self-shielding table: every carrier weighted by the H nuclei it
      ! holds, the same sum ionization_equilibrium forms for nh.
      nH_nuc = rho*n0*(f_sp(:,1) + f_sp(:,2)                             &
                     + 2.0_dp*(f_sp(:,isp_H2) + f_sp(:,isp_H2p))         &
                     + 3.0_dp*f_sp(:,isp_H3p) + f_sp(:,isp_HeHp)         &
                     + f_sp(:,isp_OH) + 2.0_dp*f_sp(:,isp_H2O))
      call fuv_lw_photon_field(nH2, nH2O, nOH, TK, nH_nuc, c1, c2, c3,   &
                               fsh, trl, tau_b, klw, pls, pla,           &
                               jh2o, joh, amax, ovl)
      do j = 1-Ng, N+Ng
         cph_klw(j)      = klw(j)
         cph_jh2o(j,:)   = jh2o(j,:)
         cph_joh(j,:)    = joh(j,:)
      enddo
      end subroutine carrier_photolysis

      ! ------------------------------------------------------------- !

      ! The H nuclei cell j's transported carriers may hold [cm^-3].
      !
      ! It is NOT the hydrogen element total.  The H sitting in H II, H2+ and
      ! H3+ is frozen for the duration of this step (sec. 4), and the
      ! write-back does not re-solve those stages -- it rescales them in
      ! proportion to what they already held.  So a carrier state that used
      ! the whole element would leave the closure species nothing: the
      ! remainder handed to H I, H II, H2+ and H3+ goes to zero, the free
      ! electron density with it, and the very next ionization sweep is given
      ! a cell it cannot find a root in.  That is not hypothetical -- it is
      ! what the transported hot Uranus did in 80 of 500 cells between
      ! 1.032 and 1.156 R_p, ending in "persistent non-root chemical
      ! equilibrium" at r = 1.112 with n_e = 3.7 cm^-3.
      !
      ! This is the same budget carrier_source uses to close atomic H, and
      ! the two must not be able to drift apart, so both call this.
      double precision function hydrogen_available_to_carriers(j, nH_free) &
                                result(nH)
      ! The hydrogen nuclei the carriers of this cell may hold: the element
      ! total less the nuclei sitting in the stages this step holds FROZEN
      ! and rescales rather than re-solves.
      !
      ! WHICH STAGES THOSE ARE DEPENDS ON THE RUN, and that is the whole
      ! change the proton option makes here. Without it H+ is one of the
      ! frozen stages and its nuclei are unavailable; with it H+ is a
      ! carrier, solved by the same Newton as H2, so its nuclei are on the
      ! carrier side of the budget and are subtracted by the caller through
      ! nc(ic_Hp) instead of here. Leaving cbg_nhii in the subtraction with
      ! the proton carried would charge the proton twice -- once as frozen
      ! background and once as a carrier -- and cap the carriers below the
      ! hydrogen that exists.
      integer,  intent(in) :: j
      real(dp), intent(in) :: nH_free
      if (ionization_transport) then
         nH = max(nH_free - 2.0d0*cbg_nh2p(j) - 3.0d0*cbg_nh3p(j), 0.0d0)
      else
         nH = max(nH_free - cbg_nhii(j) - 2.0d0*cbg_nh2p(j)              &
                          - 3.0d0*cbg_nh3p(j), 0.0d0)
      endif
      end function hydrogen_available_to_carriers

      ! ------------------------------------------------------------- !

      ! Chemical production minus loss of the four carriers in cell j, at the
      ! trial densities nc(1..4) [cm^-3].  The rows are the SAME rows the
      ! local equilibrium solve uses (sec. 4); nothing is rewritten here.
      subroutine carrier_source(j, nc, nH_free, nO_free, src)
      integer,  intent(in)  :: j
      real(dp), intent(in)  :: nc(n_carrier_max), nH_free, nO_free
      real(dp), intent(out) :: src(n_carrier_max)
      real(dp) :: fv(10)
      real(dp) :: n_hi, n_hii, n_h2p, n_h3p, n_hehp
      real(dp) :: n_hei, n_heii, n_heiii, n_heiTR, n_heiSI, n_e, n_o0

      ! THE PROTON IS EITHER FROZEN BACKGROUND OR THE TRIAL UNKNOWN.
      ! With the ionization state carried, n(H+) is what this Newton is
      ! solving for in this cell, so the row has to be evaluated at the
      ! trial value and not at the value the last sweep left behind --
      ! otherwise the recombination sink and the H closure below both refer
      ! to a different proton than the one the row is moving, and the
      ! Jacobian column for it is zero.
      if (ionization_transport) then
         n_hii = nc(ic_Hp)
      else
         n_hii = cbg_nhii(j)
      endif
      n_h2p   = cbg_nh2p(j)
      n_h3p   = cbg_nh3p(j)
      n_hehp  = cbg_nhehp(j)
      n_heii  = cbg_nheii(j)
      n_heiii = cbg_nheiii(j)
      n_heiTR = cbg_nheiTR(j)
      n_hei   = cbg_nhei(j)
      n_heiSI = max(n_hei - n_heiTR, 0.0d0)
      ! THE ELECTRONS THAT RECOMBINE WITH THE TRIAL PROTONS ARE THE TRIAL
      ! PROTONS' OWN. The rest of the electron budget -- He+, He++, the
      ! molecular ions, the metals -- is background and stays frozen with
      ! it, but holding the H+ electrons frozen while H+ itself moves would
      ! break charge neutrality inside the solve and, where H+ carries the
      ! electrons (the wind of this configuration), would evaluate the
      ! recombination sink a_hii n_e n_hii at an n_e belonging to a
      ! different ionization state. n_e is a sum over charges, so replacing
      ! one term by its trial value is exact.
      n_e     = cbg_ne(j)
      if (ionization_transport)                                              &
         n_e = max(cbg_ne(j) - cbg_nhii(j) + n_hii, 0.0d0)

      ! Atomic H closes the hydrogen budget: what the element has left after
      ! the ions, the molecular ions and the trial carriers.  This is the
      ! coupling that makes the H2 source a function of the unknown.
      ! nH_free is the element headroom, so the frozen ions come off here.
      ! With the proton carried it is one of the carriers this closure
      ! subtracts (one H nucleus each), and hydrogen_available_to_carriers
      ! has stopped subtracting it as frozen background.
      if (ionization_transport) then
         n_hi = max(hydrogen_available_to_carriers(j, nH_free)           &
                    - 2.0d0*nc(ic_H2) - nc(ic_OH) - 2.0d0*nc(ic_H2O)     &
                    - nc(ic_Hp), 0.0d0)
      else
      n_hi = max(hydrogen_available_to_carriers(j, nH_free)              &
                 - 2.0d0*nc(ic_H2) - nc(ic_OH) - 2.0d0*nc(ic_H2O), 0.0d0)
      endif
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
                        ieq_cell%P_HeITR, ieq_cell%P_H2,                 &
                        ieq_cell%P_H2_di, ieq_cell%P_H2_dd,              &
                        ieq_cell%P_H2_nd, ieq_cell%k_LW,                 &
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
      ! Row (1) of mol_heh_rows is the complete proton balance and is taken
      ! whole: photoionization of H including the secondaries, the
      ! dissociative and double photoionization branches of H2, the H2+ + H
      ! and He+ + H2 channels and the He(2^3S) Penning term, less radiative
      ! recombination and less the H+ + H2 reactions R10 and R13. There is
      ! no proton chemistry written here that the local solve does not
      ! already have; the difference is only where the row is evaluated.
      src(ic_Hp)  = 0.0d0
      if (ionization_transport) src(ic_Hp) = fv(1)

      end subroutine carrier_source

      ! ------------------------------------------------------------- !

      ! Face flux of one carrier and its two derivatives.  The drift term is
      ! linear in f (not the X(1-X) form of the element operator, whose two
      ! factors exist because that variable is one half of a binary pair).
      subroutine carrier_face_flux(fl, fr, wcl, wcr, Agr, Bst, upw,      &
                                   Jf, dJl, dJr)
      ! fl, fr are the unknown -- the fraction per unit mass -- and wcl, wcr
      ! turn each into the PARTICLE mixing ratio X = wc*y whose gradient
      ! drives molecular diffusion and whose value the settling drift acts
      ! on.  The two are not the same variable wherever mbar varies.
      real(dp), intent(in)  :: fl, fr, wcl, wcr, Agr, Bst
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
      Jf  = -Agr*(wcr*fr - wcl*fl) - Bst*(wl*wcl*fl + wr*wcr*fr)
      dJl =  (Agr - Bst*wl)*wcl
      dJr = (-Agr - Bst*wr)*wcr
      end subroutine carrier_face_flux

      ! ------------------------------------------------------------- !

      ! Residual of the backward-Euler transport-chemistry system, and the
      ! face-flux derivatives the Jacobian needs.  rnorm is the largest
      ! RELATIVE row imbalance, built the way composition_residual builds
      ! its own: the row divided by the sum of the magnitudes of its terms.
      subroutine carrier_residual(fc, fc_old, nrho, wfac, dt_phys, rp,   &
                                  rep, ntv, Agrd, Bdrf, updrf,           &
                                  nH_free, nO_free, nC_free, sigrate,    &
                                  res, Jf, dJl, dJr, rnorm)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)  :: fc, fc_old
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: dt_phys
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: rp, rep, ntv
      real(dp), dimension(0:N,n_carrier_max),       intent(in)  :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max),       intent(in)  :: updrf
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nC_free
      ! Signal-crossing rate (|v| + c_s)/dr of each cell [1/s].  It carries
      ! the absolute floor of the row scale below, which used to be carried
      ! by 1/dt -- see there.
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: sigrate
      real(dp), dimension(1:N,n_carrier_max),       intent(out) :: res
      real(dp), dimension(0:N,n_carrier_max),       intent(out) :: Jf, dJl, dJr
      real(dp),                                 intent(out) :: rnorm
      ! Cell and species carrying the worst relative row imbalance, for the
      ! debug line: a Newton that stops short says WHERE it stopped.
      real(dp) :: nc(n_carrier_max), src(n_carrier_max)
      real(dp) :: Kj, sL, sR, cadv, dsc, rr, nfl
      integer  :: j, ic
      ! Row scale of every row of this assembly, kept so that the worst
      ! relative imbalance can be located AFTER the cell loop rather than
      ! inside it -- see there.
      real(dp) :: rsc(1:N,n_carrier_max)

      res = 0.0d0
      Jf  = 0.0d0
      dJl = 0.0d0
      dJr = 0.0d0
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N-1
            call carrier_face_flux(fc(j,ic), fc(j+1,ic), wfac(j),        &
                                   wfac(j+1), Agrd(j,ic),                &
                                   Bdrf(j,ic), updrf(j,ic),              &
                                   Jf(j,ic), dJl(j,ic), dJr(j,ic))
         enddo
      enddo

      ! CELL-LOCAL, AND THE ONE THING THAT WAS NOT.  Every term of a carrier
      ! row is built from this cell's own state and its two stored face
      ! fluxes, and carrier_source evaluates the chemistry at this cell
      ! through module state that is threadprivate (ieq_cell, the mk* rate
      ! coefficients of System_HeH_mol and the oj* photolysis channels), so
      ! each thread fills its own copy per cell.  The exception was the
      ! running maximum: rnorm together with the cell and carrier that
      ! attain it, which is a max and an argmax and cannot be taken inside a
      ! parallel loop without making the diagnostic depend on which thread
      ! got there first.
      !
      ! So the loop takes the maximum only -- a max is exact and independent
      ! of the order it is taken in, so rnorm is bitwise the sequential value
      ! at any number of threads -- and the cell that attains it is found by
      ! a scan afterwards, which is what makes pct_worst_j / pct_worst_ic
      ! deterministic.  The sequential loop updated on a STRICT improvement,
      ! so its winner was the FIRST row in (j, then ic) order attaining the
      ! maximum, and the scan below stops at exactly that row.  With no row
      ! above the initial zero it leaves both untouched, as the sequential
      ! loop did.
      rnorm = 0.0d0
      rsc   = 0.0d0
      !$omp parallel do default(shared) schedule(static)                 &
      !$omp   private(j,ic,nc,src,Kj,sL,sR,cadv,dsc,rr,nfl)              &
      !$omp   reduction(max:rnorm)
      do j = 1, N
         nc = fc(j,:)*nrho(j)
         call carrier_source(j, nc, nH_free(j), nO_free(j), src)
         Kj = 1.0d0/(rp(j)**2*max(rep(j) - rep(j-1), 1.0d0))
         sL = rep(j-1)**2
         sR = rep(j)**2
         do ic = 1, n_carrier
            ! A carrier this run does not solve keeps a zero row and an
            ! identity block, so its unknown never moves and it contributes
            ! nothing to the residual norm.
            if (.not. carrier_solved(ic)) then
               res(j,ic) = 0.0d0
               rsc(j,ic) = 0.0d0
               if (allocated(row_terms)) row_terms(j,ic) = 0.0d0
               cycle
            endif
            rr  = nrho(j)*(fc(j,ic) - fc_old(j,ic))/dt_phys(j)
            dsc = abs(nrho(j)*fc(j,ic)/dt_phys(j))                       &
                + abs(nrho(j)*fc_old(j,ic)/dt_phys(j))
            rr  = rr + Kj*(sR*Jf(j,ic) - sL*Jf(j-1,ic))
            dsc = dsc + Kj*(sR*abs(Jf(j,ic)) + sL*abs(Jf(j-1,ic)))
            ! One-sided upwind on the CELL velocity; the base face carries
            ! no advective term, which is the zero-flux partition boundary.
            if (j .eq. 1) then
               ! THE BASE CELL: direction from the base FACE, magnitude from
               ! the cell (sec. 6).  With an inflowing face the upstream is
               ! the ghost, and the partition it carries is either the
               ! handoff's or -- with no handoff -- the base cell's own,
               ! which is D6's zero-gradient statement and contributes
               ! nothing.
               if (base_inflow) then
                  if (base_dirichlet(ic)) then
                     cadv = abs(ntv(1))/max(rp(1) - rp(0), 1.0d0)
                     rr   = rr + cadv*(fc(1,ic) - fc_base(ic))
                     dsc  = dsc + abs(cadv)*max(abs(fc(1,ic)),           &
                                                abs(fc_base(ic)))
                  endif
               else
                  cadv = -abs(ntv(1))/max(rp(2) - rp(1), 1.0d0)
                  rr   = rr + cadv*(fc(2,ic) - fc(1,ic))
                  dsc  = dsc + abs(cadv)*max(abs(fc(1,ic)),              &
                                             abs(fc(2,ic)))
               endif
            else if (ntv(j) .ge. 0.0d0) then
               cadv = ntv(j)/max(rp(j) - rp(j-1), 1.0d0)
               rr   = rr + cadv*(fc(j,ic) - fc(j-1,ic))
               dsc  = dsc + abs(cadv)*max(abs(fc(j,ic)),                 &
                                          abs(fc(j-1,ic)))
            else
               if (j .lt. N) then
                  cadv = ntv(j)/max(rp(j+1) - rp(j), 1.0d0)
                  rr   = rr + cadv*(fc(j+1,ic) - fc(j,ic))
                  dsc  = dsc + abs(cadv)*max(abs(fc(j,ic)),              &
                                             abs(fc(j+1,ic)))
               endif
            endif
            ! The second-order half, deferred: a frozen source, so the
            ! Jacobian below stays the first-order block-tridiagonal one.
            if (allocated(adv_corr)) then
               rr  = rr + ntv(j)*adv_corr(j,ic)
               dsc = dsc + abs(ntv(j)*adv_corr(j,ic))
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
            ! volumetric-rate units: a density that small cannot change any
            ! observable, so a row below it IS converged.
            !
            ! IT IS CONVERTED BY THE LARGER OF TWO RATES, and the second one
            ! is why this line changed. Dividing by the step alone made the
            ! floor vanish in the STEADY evaluation, which substitutes
            ! dt = 1e30 to kill the time term (carrier_steady_residual): the
            ! floor went with it, and the far wind was then measured with no
            ! floor at all. The cell's signal-crossing rate is the same
            ! statement written without a step, so it survives dt -> infinity
            ! and agrees with the step form whenever the step is the cell's
            ! own crossing time, which is what the relaxation uses.
            if (ic .eq. ic_H2) then
               nfl = nH_free(j)
            else if (ic .eq. ic_CO) then
               nfl = min(nO_free(j), nC_free(j))
            else
               nfl = nO_free(j)
            endif
            dsc = dsc + 1.0d-20*nfl*max(1.0d0/dt_phys(j), sigrate(j))
            res(j,ic)  = rr
            rsc(j,ic)  = dsc
            if (allocated(row_terms)) row_terms(j,ic) = dsc
            if (abs(rr)/max(dsc, 1.0d-300) .gt. rnorm)                   &
               rnorm = abs(rr)/max(dsc, 1.0d-300)
         enddo
      enddo
      !$omp end parallel do

      if (rnorm .gt. 0.0d0) then
         outer: do j = 1, N
            do ic = 1, n_carrier
               if (.not. carrier_solved(ic)) cycle
               if (abs(res(j,ic))/max(rsc(j,ic), 1.0d-300)               &
                   .eq. rnorm) then
                  pct_worst_j  = j
                  pct_worst_ic = ic
                  exit outer
               endif
            enddo
         enddo outer
      endif
      end subroutine carrier_residual

      ! ------------------------------------------------------------- !

      ! Newton solve of the coupled system, block-tridiagonal in space with
      ! 4x4 blocks.  The chemistry Jacobian is a forward difference of the
      ! same rows the residual calls, so a change to the network reaches the
      ! Jacobian without a second edit.
      subroutine solve_carriers(fc, fc_old, nrho, wfac, dt_phys, rp,     &
                                rep, ntv, Agrd, Bdrf, updrf, TK, sigrate,&
                                nH_free, nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout) :: fc
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)    :: fc_old
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nrho, wfac
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: dt_phys
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rp, rep
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: ntv, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: sigrate
      real(dp), dimension(0:N,n_carrier_max),       intent(in)    :: Agrd, Bdrf
      integer,  dimension(0:N,n_carrier_max),       intent(in)    :: updrf
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free

      real(dp), dimension(1:N,n_carrier_max)  :: res, rhs
      real(dp), dimension(0:N,n_carrier_max)  :: Jf, dJl, dJr
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max) :: aa, bb, cc
      real(dp), dimension(1:N,n_carrier_max)  :: dfc
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: ftry
      real(dp) :: nc(n_carrier_max), s0(n_carrier_max), s1(n_carrier_max)
      real(dp) :: Kj, sL, sR, cadv, dn, nref
      real(dp) :: rnorm, rprev, rstart, rtry, damp
      integer  :: j, ic, kc, it, ihalf

      call carrier_residual(fc, fc_old, nrho, wfac, dt_phys, rp, rep,    &
                            ntv, Agrd, Bdrf, updrf, nH_free, nO_free,    &
                            nC_free, sigrate, res, Jf, dJl, dJr, rnorm)
      rstart = rnorm
      rprev  = rnorm
      pct_newton_steps = 0

      do it = 1, newton_maxit
         if (rnorm .le. newton_tol) exit

         ! ---- assemble the block tridiagonal
         aa = 0.0d0
         bb = 0.0d0
         cc = 0.0d0
         ! CELL-LOCAL ASSEMBLY.  Row j of the block-tridiagonal system is
         ! built from cell j's own geometry, its two stored face-flux
         ! derivatives and n_carrier + 1 evaluations of carrier_source at
         ! cell j, and it writes aa/bb/cc at index j alone.  carrier_source
         ! reaches the chemistry through threadprivate module state, so each
         ! thread carries its own rate coefficients.  There is no reduction
         ! here, so the assembled blocks are bitwise the sequential ones at
         ! any number of threads.  The Thomas sweep that consumes them is
         ! sequential in space and stays serial.
         !$omp parallel do default(shared) schedule(static)              &
         !$omp   private(j,ic,kc,Kj,sL,sR,cadv,nc,s0,s1,dn,nref)
         do j = 1, N
            Kj = 1.0d0/(rp(j)**2*max(rep(j) - rep(j-1), 1.0d0))
            sL = rep(j-1)**2
            sR = rep(j)**2
            do ic = 1, n_carrier
               bb(j,ic,ic) = nrho(j)/dt_phys(j)
               ! The identity row of a carrier this run does not solve: the
               ! diagonal above with a zero right-hand side gives it a zero
               ! update, which is what keeps it out of the system without
               ! taking it out of the block.
               if (.not. carrier_solved(ic)) cycle
               if (j .lt. N) then
                  bb(j,ic,ic) = bb(j,ic,ic) + Kj*sR*dJl(j,ic)
                  cc(j,ic,ic) = cc(j,ic,ic) + Kj*sR*dJr(j,ic)
               endif
               if (j .gt. 1) then
                  bb(j,ic,ic) = bb(j,ic,ic) - Kj*sL*dJr(j-1,ic)
                  aa(j,ic,ic) = aa(j,ic,ic) - Kj*sL*dJl(j-1,ic)
               endif
               if (j .eq. 1) then
                  if (base_inflow) then
                     if (base_dirichlet(ic)) then
                        ! The ghost is data, not an unknown: the term is on
                        ! the diagonal alone.
                        cadv = abs(ntv(1))/max(rp(1) - rp(0), 1.0d0)
                        bb(1,ic,ic) = bb(1,ic,ic) + cadv
                     endif
                  else
                     cadv = -abs(ntv(1))/max(rp(2) - rp(1), 1.0d0)
                     bb(1,ic,ic) = bb(1,ic,ic) - cadv
                     cc(1,ic,ic) = cc(1,ic,ic) + cadv
                  endif
               else if (ntv(j) .ge. 0.0d0) then
                  cadv = ntv(j)/max(rp(j) - rp(j-1), 1.0d0)
                  bb(j,ic,ic) = bb(j,ic,ic) + cadv
                  aa(j,ic,ic) = aa(j,ic,ic) - cadv
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
            nc = fc(j,:)*nrho(j)
            call carrier_source(j, nc, nH_free(j), nO_free(j), s0)
            do kc = 1, n_carrier
               if (.not. carrier_solved(kc)) cycle
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
                  if (.not. carrier_solved(ic)) cycle
                  bb(j,ic,kc) = bb(j,ic,kc)                              &
                              - (s1(ic) - s0(ic))/dn*nrho(j)
               enddo
            enddo
         enddo
         !$omp end parallel do

         rhs = -res
         call block_thomas(aa, bb, cc, rhs, dfc)

         ! ---- damped update, with f >= 0 enforced on every try
         damp = 1.0d0
         do ihalf = 0, newton_halves
            ftry = fc
            do ic = 1, n_carrier
               if (.not. carrier_solved(ic)) cycle
               do j = 1, N
                  ftry(j,ic) = max(fc(j,ic) + damp*dfc(j,ic), 0.0d0)
               enddo
               ! Only the OUTER ghosts mirror the interior.  The lower
               ! ghost is the inflow reservoir, not a copy of the base
               ! cell: it carries the composition the handoff states, and
               ! carrier_base_state read it before the Newton started.
               ! Writing the base cell into it made the base Dirichlet
               ! degenerate into a zero-gradient condition from the second
               ! relaxation pass on -- the marching loop hid that because
               ! its ionization sweep re-pins the ghost every step, and
               ! relax_photochemical_composition, which takes many
               ! transport steps between two sweeps, did not (sec. 158).
               ftry(N+1:N+Ng,ic) = ftry(N,ic)
            enddo
            call carrier_residual(ftry, fc_old, nrho, wfac, dt_phys, rp, &
                                  rep, ntv, Agrd, Bdrf, updrf, nH_free,  &
                                  nO_free, nC_free, sigrate, res, Jf,    &
                                  dJl, dJr, rtry)
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

      call limit_to_element_budget(fc, nrho, TK, nH_free, nO_free,        &
                                   nC_free)

      end subroutine solve_carriers

      ! ------------------------------------------------------------- !

      ! Block-tridiagonal Thomas sweep with dense n_carrier x n_carrier
      ! blocks.  Written here because the element operator's scalar Thomas
      ! sweep does not generalize: the species system is diagonal in space
      ! but dense in species through the chemistry.
      subroutine block_thomas(aa, bb, cc, dd, xx)
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max), intent(in)  :: aa,bb,cc
      real(dp), dimension(1:N,n_carrier_max),           intent(in)  :: dd
      real(dp), dimension(1:N,n_carrier_max),           intent(out) :: xx
      real(dp), dimension(1:N,n_carrier_max,n_carrier_max) :: cp
      real(dp), dimension(1:N,n_carrier_max)           :: dp_
      real(dp) :: mm(n_carrier_max,n_carrier_max)
      real(dp) :: rhs(n_carrier_max,n_carrier_max+1)
      real(dp) :: piv, fac
      integer  :: j, i, k, l, ip

      cp  = 0.0d0
      dp_ = 0.0d0
      xx  = 0.0d0
      do j = 1, N
         ! mm = bb(j) - aa(j) cp(j-1)
         mm(1:n_carrier,1:n_carrier) = bb(j,1:n_carrier,1:n_carrier)
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
         rhs(1:n_carrier,1:n_carrier) = cc(j,1:n_carrier,1:n_carrier)
         rhs(1:n_carrier,n_carrier+1) = dd(j,1:n_carrier)
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
         cp(j,1:n_carrier,1:n_carrier) = rhs(1:n_carrier,1:n_carrier)
         dp_(j,1:n_carrier)  = rhs(1:n_carrier,n_carrier+1)
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
      subroutine limit_to_element_budget(fc, nrho, TK, nH_free, nO_free,  &
                                         nC_free)
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(inout) :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nrho, TK
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      real(dp) :: nc(n_carrier_max), got, sc, over, nco_eq, nH_car
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
      ! The OUTER ghosts take the partition of the cell they mirror -- the
      ! zero-gradient outflow boundary -- and then go through the same
      ! limiter as an interior cell.  They HAVE to: a ghost carries its own
      ! density, so the same fraction is a different number of nuclei there
      ! and can ask for more than the ghost's element holds.  Leaving them
      ! out of the limiter was measured, and it broke the oxygen and carbon
      ! budgets by 2% in a 400-step run while every interior cell closed at
      ! round-off.  The LOWER ghosts are not mirrored and not limited: this
      ! operator does not own them (see solve_carriers).
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         fc(N+1:N+Ng,ic) = fc(N,ic)
      enddo
      do j = 1, N+Ng
         nc  = max(fc(j,:)*nrho(j), 0.0d0)
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
         ! hydrogen: scale the three H-bearing carriers together, against
         ! the SAME budget the chemistry rows close atomic H with
         ! (hydrogen_available_to_carriers) -- not against the element
         ! total, which would let the carriers take the nuclei the frozen
         ! ion stages are holding and leave the write-back nothing to give
         ! back to them.
         ! The proton holds ONE H nucleus and, where it is carried, it is
         ! in this sum with the others. It is bounded the same way -- the
         ! whole H-bearing carrier set scaled by one factor when it asks for
         ! more hydrogen than the cell has -- rather than clamped on its
         ! own: a clamp on one carrier is a kink in the residual at the
         ! surface where it binds, and section 162 measured what a kink at
         ! a front does to a Newton. Scaling the set keeps the state inside
         ! the simplex along a straight line through the interior, which is
         ! the same construction the positivity limiter of the
         ! reconstruction uses.
         !
         ! WHERE IT BINDS IS THE OTHER END OF THE WIND FROM THE H2 FRONT.
         ! Below ~1.5 r_base the hydrogen is nearly all neutral and the
         ! proton asks for a small share; the budget is tight where
         ! x(H+) -> 1, high in the wind, and there H2 has long since gone,
         ! so the two carriers are rarely tight at the same place.
         nH_car = hydrogen_available_to_carriers(j, nH_free(j))
         got = 2.0d0*nc(ic_H2) + nc(ic_OH) + 2.0d0*nc(ic_H2O)
         if (carrier_solved(ic_Hp)) got = got + nc(ic_Hp)
         if (got .gt. nH_car) then
            over = got/max(nH_car, 1.0d-300) - 1.0d0
            pct_worst_limit = max(pct_worst_limit, over)
            sc = nH_car/max(got, 1.0d-300)
            nc(ic_H2)  = nc(ic_H2) *sc
            nc(ic_OH)  = nc(ic_OH) *sc
            nc(ic_H2O) = nc(ic_H2O)*sc
            if (carrier_solved(ic_Hp)) nc(ic_Hp) = nc(ic_Hp)*sc
            if (over .gt. limit_report) hit = .true.
         endif
         if (hit) pct_limited = pct_limited + 1
         fc(j,:) = nc/max(nrho(j), 1.0d-99)
      enddo
      end subroutine limit_to_element_budget

      ! ------------------------------------------------------------- !

      ! Write the solved carriers back into f_sp, moving the nucleus
      ! difference into the closure species of each element so that no
      ! element total and no mass density is changed by the step.  bsp_mass
      ! of a carrier is exactly the sum of its nuclei's masses, so the mass
      ! moved out of H I, O I and C I is the mass moved into the carriers.
      subroutine carrier_write_back(rho, f_sp, fc, nrho, nH_free,        &
                                    nO_free, nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max), intent(in)    :: fc
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nrho
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: nC_free
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp) :: nd, nH2n, nOHn, nH2On, nCOn, nHpn, rest, held, sc
      integer  :: j, k, i0

      ! THE LOWER GHOSTS ARE NOT WRITTEN.  They hold the composition of the
      ! gas flowing in, which the ionization sweep pins from the handoff and
      ! which this operator reads as its base Dirichlet value; rewriting
      ! them with the base cell's own partition would erase the boundary
      ! condition the step is being driven by.
      do j = 1, N+Ng
         nd = rho(j)*n0
         if (nd .le. 0.0d0) cycle
         nH2n  = fc(j,ic_H2) *nrho(j)
         nOHn  = fc(j,ic_OH) *nrho(j)
         nH2On = fc(j,ic_H2O)*nrho(j)
         nCOn  = fc(j,ic_CO) *nrho(j)
         nHpn  = 0.0d0
         if (carrier_solved(ic_Hp)) nHpn = fc(j,ic_Hp)*nrho(j)
         f_sp(j,isp_H2)  = nH2n /nd
         f_sp(j,isp_OH)  = nOHn /nd
         f_sp(j,isp_H2O) = nH2On/nd
         f_sp(j,isp_CO)  = nCOn /nd
         if (carrier_solved(ic_Hp)) f_sp(j,isp_HII) = nHpn/nd

         ! The nuclei the carriers did NOT take are shared out over the
         ! element's other species IN PROPORTION TO WHAT THEY ALREADY HELD.
         ! Two reasons for the proportional rule rather than putting the
         ! whole remainder in the neutral stage: it cannot drive a stage
         ! negative, and it leaves the ionization fraction of the cell where
         ! it was, which is the right starting point for the sweep that
         ! re-solves it a moment later.  HeH+ is excluded from the hydrogen
         ! share because rescaling it would move a helium nucleus too.
         ! The nuclei the carriers did not take, shared over the stages
         ! that are still local. WHEN THE PROTON IS CARRIED IT LEAVES THAT
         ! SHARE: its density has just been written from the transported
         ! value above, and rescaling it here in proportion to what it held
         ! would immediately undo the transport. It moves to the taken side
         ! of the ledger instead, one nucleus each, exactly as H2 takes two.
         if (carrier_solved(ic_Hp)) then
            rest = nH_free(j) - 2.0d0*nH2n - nOHn - 2.0d0*nH2On - nHpn
            held = f_sp(j,isp_HI)*nd                                     &
                 + 2.0d0*f_sp(j,isp_H2p)*nd + 3.0d0*f_sp(j,isp_H3p)*nd
         else
         rest = nH_free(j) - 2.0d0*nH2n - nOHn - 2.0d0*nH2On
         held = f_sp(j,isp_HI)*nd + f_sp(j,isp_HII)*nd                   &
              + 2.0d0*f_sp(j,isp_H2p)*nd + 3.0d0*f_sp(j,isp_H3p)*nd
         endif
         if (held .gt. 0.0d0 .and. rest .gt. 0.0d0) then
            sc = rest/held
            f_sp(j,isp_HI)  = f_sp(j,isp_HI) *sc
            if (.not. carrier_solved(ic_Hp))                             &
               f_sp(j,isp_HII) = f_sp(j,isp_HII)*sc
            f_sp(j,isp_H2p) = f_sp(j,isp_H2p)*sc
            f_sp(j,isp_H3p) = f_sp(j,isp_H3p)*sc
         else if (held .le. 0.0d0) then
            f_sp(j,isp_HI) = max(rest, 0.0d0)/nd
         else
            f_sp(j,isp_HI)  = 0.0d0
            if (.not. carrier_solved(ic_Hp)) f_sp(j,isp_HII) = 0.0d0
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
                                                 trust, drift, nstep)
      ! Advance the carriers at the fixed wind by AT MOST `trust` of the
      ! largest H2 mixing ratio on the grid, and hand that state over whole.
      !
      ! WHAT THIS REPLACES, AND WHY (docs/p50_carrier_wind_alternation.md).
      ! It used to relax all the way to the fixed-wind steady state and then
      ! blend the result into the entry state at an under-relaxation factor.
      ! Both halves of that were measured on the He/H = 0.0793 hot Uranus,
      ! from a wind the steady solver had just accepted at info = 0:
      !
      !  * The fixed-wind steady state is not near the entry state and is not
      !    a solution of anything. It drives x2 = 2n(H2)/n_H to the hydrogen
      !    ELEMENT CEILING, exactly 1.0, in 123 of 500 cells between 1.047
      !    and 1.227 R_p -- a state pinned by limit_to_element_budget over a
      !    quarter of the grid. Handing it over threw ||R|| from 1.5e-4 to
      !    1.9e-2 and the next solve aborted.
      !  * The blend is worse than either state it interpolates. Carrier
      !    steady residual: entry 1.6e-4, relaxed exit 2.2e-3, the omega=0.5
      !    blend 8.6e-3, and it threw ||R|| to 6.7e-2 rather than 1.9e-2.
      !    A blend of two states solves neither equation. Reducing omega did
      !    not help and is measured to hurt (omega = 1.8e-3 gave 2.0e-1).
      !  * A bounded advance, handed over whole, does work: at 1e-2 the
      !    steady solve returned info = 0 on 20 consecutive passes with
      !    ||R|| between 2.4e-4 and 3.1e-4, where the blend aborted on the
      !    first. The bound is where the margin is: the solve still accepts
      !    at 0.2 and stops accepting at 0.4, so 1e-2 stands a factor 20
      !    inside the boundary, and it is met by the FIRST step of the
      !    operator (one cell crossing time moves the composition by 1.1e-2).
      !
      ! What a bounded pass does NOT do is reach a fixed point, and nothing
      ! here claims it: the H2 front advances about one cell in each pass, at
      ! 83 per cent of the local gas speed. That is the front of sections
      ! 128.6 and 134.4 still propagating, and it is why the carrier
      ! transport stays default off.
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: rho, v
      real(dp), dimension(1-Ng:N+Ng),           intent(in)    :: Tcode
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real(dp),                                 intent(in)    :: trust
      real(dp),                                 intent(out)   :: drift
      integer,                                  intent(out)   :: nstep

      real(dp), dimension(1-Ng:N+Ng) :: dt_code
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: fprev, fnow, fentry
      real(dp), dimension(1-Ng:N+Ng,n_carrier_max) :: Dco
      real(dp), dimension(1-Ng:N+Ng) :: ntot, TK, mbar, nrho, wfac
      real(dp), dimension(1-Ng:N+Ng) :: nH_free, nO_free, nC_free
      real(dp) :: drj, tdiff, tadv, grow, tscale, dmax, x_ref
      integer  :: j, k, ic
      type(element_census_state) :: cen_relax

      drift = 0.0d0
      nstep = 0
      if (.not. thereis_mol)      return
      if (.not. carrier_transport) return
      if (.not. bg_ready)         return

      call element_census_take('relax_photochemical_composition',         &
                               rho, f_sp, cen_relax)
      tscale = R0/v0
      call carrier_state(rho, Tcode, f_sp, fprev, ntot, nrho, wfac, TK, &
                         mbar, nH_free, nO_free, nC_free)
      call carrier_base_state(fprev, rho, v)
      fentry = fprev
      call carrier_diffusivities(f_sp, rho, TK, ntot, Dco)
      do j = 1, N
         drj   = max((r_edg(j) - r_edg(j-1))*R0, 1.0d0)
         tdiff = drj*drj/max(maxval(Dco(j,1:n_carrier))               &
                            + kzz_cell(j), 1.0d-30)
         tadv  = drj/max(abs(v(j))*v0, 1.0d-30)
         dt_code(j) = min(tdiff, tadv)/tscale
      enddo
      dt_code(1-Ng:0)   = dt_code(1)
      dt_code(N+1:N+Ng) = dt_code(N)

      ! THE BOUND IS ON THE COMPOSITION, NOT ON THE TIME, and the difference
      ! is what makes it a bound at all.  Section 128 tried the time form --
      ! cap how many flow times one pass may advance -- and measured that it
      ! controls nothing: at one tenth of a cell's flow time the pass still
      ! crosses about sixty cells at the H2 front and the move is of order
      ! unity either way (drifts 0.536 and 0.533 at caps of 1.0 and 0.1).
      ! The front is advected, and advection is not slow on any time scale
      ! such a bound can be written in.  It IS slow on a composition scale:
      ! the first step, at the cell's own crossing time, moves the largest
      ! mixing ratio by 1.1e-2, so a bound of 1e-2 is met immediately and a
      ! bound of 0.4 is not met until the pass is most of the way to the
      ! clamped ceiling state.
      grow  = 1.0d0
      do k = 1, relax_maxstep
         call photochemical_transport_step(rho, v, Tcode, f_sp,          &
                                           dt_code*grow)
         call carrier_state(rho, Tcode, f_sp, fnow, ntot, nrho, wfac,    &
                            TK, mbar, nH_free, nO_free, nC_free)
         nstep = k
         dmax  = 0.0d0
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            do j = 1, N
               dmax = max(dmax, abs(fnow(j,ic) - fprev(j,ic)))
            enddo
         enddo
         x_ref = max(maxval(fnow(1:N,ic_H2)), 1.0d-30)
         drift = 0.0d0
         do ic = 1, n_carrier
            if (.not. carrier_solved(ic)) cycle
            do j = 1, N
               drift = max(drift, abs(fnow(j,ic) - fentry(j,ic)))
            enddo
         enddo
         if (drift/x_ref .gt. trust) exit
         if (dmax/x_ref .lt. relax_tol) exit
         fprev = fnow
         if (grow .lt. 1.0d12) grow = grow*1.5d0
      enddo

      ! How far the pass moved, on the same measure the element relaxation
      ! reports (its own drift is max|X_relaxed - X_old|/X_base).  The
      ! normalization is the largest H2 mixing ratio anywhere on the grid,
      ! one number and not the cell's own value, so this is an ABSOLUTE
      ! measure: a cell holding 1e-5 of the gas cannot report a large change
      ! by being small.  It is reported, and it is not a convergence gate --
      ! see the caller.
      drift = 0.0d0
      pct_drift_j  = 0
      pct_drift_ic = 1
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         do j = 1, N
            if (abs(fnow(j,ic) - fentry(j,ic)) .gt. drift) then
               drift = abs(fnow(j,ic) - fentry(j,ic))
               pct_drift_j  = j
               pct_drift_ic = ic
            endif
         enddo
      enddo
      drift = drift/max(maxval(fentry(1:N,ic_H2)),                       &
                        maxval(fnow(1:N,ic_H2)), 1.0d-30)

      call element_census_verify(cen_relax, rho, f_sp, rho_is_fixed=.true.)

      end subroutine relax_photochemical_composition

      ! End of module

      end module diffusive_photochemistry

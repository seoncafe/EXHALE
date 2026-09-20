   module ionization_equilibrium
	! Evaluate the ionization structure and the heating and cooling functions for a given temperature

	use mol_rates, only: h2_thermochemistry_init, h2_thermochemistry_ready
	use global_parameters
   use ion_cell_state, only: ieq_cell, ion_rates
   use species_table, only: n_mion, mion_fsp, n_melem, melem_i0,        &
                            melem_top, mion_stage,                       &
                            isp_HII, isp_HeII, isp_HeIII,                &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp,          &
                            isp_OH, isp_H2O, isp_CO, iel_O, iel_C
   use utils
   use utils_ion_eq
   use Cooling_Coefficients      ! eval_cool, recombination/ionization rates
   use System_HeH                ! Equilibrium equations
	use System_HeH_TR
	use System_HeH_mol            ! molecular network
	use System_HeH_mol_metals     ! merged molecular network + metals
	use lower_column, only: q_h2_equilibrium
	! The lower boundary's own module, for the two quantities this sweep
	! MEASURES in the ghost it solves: the ghost's particle counts and the
	! closure of its imposed molecular partition. Both are statements about
	! the state; the prescribed reservoir is not written from here.
	use base_boundary, only: set_base_ghost_counts, set_base_ghost_closure, &
	                         set_base_ghost_state_record,                 &
	                         set_base_ghost_fixed_point,                  &
	                         set_previous_sweep_ghost_rows,               &
	                         ghost_composition_seed_armed,                &
	                         lower_ghost_seed_rows, write_ghost_seed_rows
	use composition, only: base_h2_nuclei_fraction,                       &
	                       base_h2_composition_imposed
	use element_census, only: element_census_state, element_census_take,  &
	                          element_census_verify
	use h2_vibrational_relaxation, only: h2_vibrational_heat_fraction,     &
	                                     h2_energy_per_bound_fluorescence_eV
	! Oxygen chemistry (the A2 option, docs/a2_oxygen_option_design.md):
	! the FUV photolysis bands and the CO reservoir.
	use water_photolysis, only: n_fuv_band
	use oxygen_rates, only: co_equilibrium_density,                       &
	                        oxygen_chemical_equilibrium_fractions
   use System_HeH_metals
   use System_HeH_TR_metals      ! merged He-triplet + metals system
	! Constrained element-conserving chemical equilibrium in positive
	! species densities, tracked from the radiation-free molecular limit to
	! the cell's field (docs/supersonic_molecular_base.md section 11.5-B).
	use constrained_chemical_equilibrium, only:                           &
	                             equilibrium_from_molecular_limit
   use charge_exchange, only: cx_set_cell, cx_metal_base,  &   ! Huang Table 4 charge exchange
                              he_h_cx_rates,               &   ! He <-> H pair (group B)
                              cx_add_to_turnover               ! cx bound of a row's turnover
   use System_H
   use newton_solver, only: solve_ieq   ! Task 2: analytic-Jacobian Newton (+ hybrd1 fallback)
   use opacity_models            ! opacity_pT_factor for the 'P' model

   implicit none

	! WIDTH, IN CELLS, OF THE BLOCK THE XUV FIELD IS HELD FIXED OVER while
	! the equilibrium sweep solves it (ioniz_eq, xuv_block_H / xuv_block_HHe).
	! Zero, the default, means ONE block: the whole grid, i.e. the field of
	! the composition the sweep was handed. Any positive width makes the
	! sweep carry the columns of the cells it has already solved from one
	! block to the next, which is a substitution in the star-ward direction.
	!
	! IT CHANGES NO FIXED POINT. At the fixed point a cell's entry and
	! returned composition are one state, so every width builds the same
	! columns and the same field. What the width trades is the pass count of
	! the outer coupled loop against the width of the parallel front: inside
	! a block the cells do not depend on each other and are solved together,
	! across blocks they are solved in order.
	!
	! MEASURED 2026-09-07, wasp_full at 300 steps, one private binary per
	! width, wall time at one and at sixteen threads:
	!
	!     width        passes    1 thread   16 threads
	!     whole grid    10.10     136.3 s      17.0 s
	!        64          9.51     133.4 s      18.6 s
	!        16          8.05     109.4 s      28.0 s
	!         8          7.13     103.4 s      42.5 s
	!         4          7.14      99.0 s      47.9 s
	!         2          7.19     106.3 s      62.2 s
	!         1          8.10     112.3 s      99.3 s
	!
	! THE SUBSTITUTION DOES NOT PAY, and the reason is in the table rather
	! than in the physics. It removes the part of the lag that walks the
	! ionization front inward one cell per pass, and the pass count falls
	! from 10.10 to 7.13 at its best width; but a pass costs the same field
	! work either way, and the field of a whole grid is a single loop over
	! independent cells that a machine solves all at once. A narrow block
	! gives that width up, and sixteen threads recover more than the pass
	! count saves. The width that keeps the parallel front (64) has almost
	! no pass saving left. Hence zero.
	!
	! The width is a compile-time constant and NOT a function of the thread
	! count, deliberately: the run's answer must not depend on how many
	! threads it was given, and the widths above reach measurably different
	! states at the 1e-8 fixed-point tolerance while the thread count moves
	! the present code by 7e-14 (MEASURED on the same case).
	integer, parameter :: xuv_field_block_cells = 0

	! HOW MANY TIMES A CELL'S OWN ATTENUATION AND ITS OWN COMPOSITION ARE
	! SOLVED AGAINST EACH OTHER before the sweep leaves the cell.
	!
	! A cell's field is exp(-tau_out) times the cell mean of its OWN
	! attenuation, (1 - exp(-dtau))/dtau, and dtau = sum_abs sigma_nu n_abs
	! dr is built from the very densities the cell is being solved for
	! (photoionization_field_at_cell_HHe).  One pass therefore solves the
	! chemistry in the field of the composition the cell was HANDED, and the
	! composition it returns implies a slightly different field.  A second
	! pass re-forms the field from the returned densities and solves again;
	! COST3 MEASURED that mode contracting by 0.17 a pass on wasp_full, so a
	! handful of passes exhausts it.
	!
	! ONE, the default, is one field build and one solve, which is the
	! scheme this code has always run.  IT CHANGES NO FIXED POINT: at the
	! fixed point the composition a cell enters with and the one it returns
	! are one state, so the field built from either is the same field, and
	! the outer coupled loop of EXHALE_main converges to the same pair.
	! What a larger count trades is the pass count of that outer loop
	! against the field and solve work inside the cell.
	!
	! THAT ARGUMENT IS ABOUT A FIXED POINT AND NOT ABOUT THEIR NUMBER, and
	! the sweep of a shielded molecular column has more than one. MEASURED
	! on the 1 microbar hot Uranus reloaded (mol_diffusion, He diffusion on):
	! the sweep repeated at one state from two admissible seeds reaches two
	! compositions that are each exact roots of every cell's network
	! (normalized reaction residual 1e-25, acceptance class 1, no amnesty)
	! and each self-consistent with the field it builds, differing by a
	! factor 8.9 in x(H2) at the H2 front and 19 percent in x(H I) across
	! the ionization front. Raising this count to 20 leaves that difference
	! unchanged to every digit, so it is not the cell's own field lag.
	! Around ONE such composition the sweep is also marginal in the total
	! molecular hydrogen content of the shielded layer: the fast chemistry
	! there cycles H2 -> H2+ -> H3+ and back without changing the H2 nuclei
	! total, so the local equilibrium fixes only the partition and 0.77 of
	! any seed perturbation of the total survives every pass (MEASURED
	! identically at seed perturbations of 1e-6 and 1e-2). The content is
	! upstream and transport data, which is what `Molecular carrier
	! transport` makes it: with H2 a Newton unknown the same measurement
	! falls from 3.2e-03 to 9.7e-10 in the acceptance measure.
	!
	! MEASURED 2026-09-07, one private binary per count, wasp_full at 300
	! steps, wall time at one and at sixteen threads:
	!
	!     passes   coupled passes   1 thread   16 threads
	!        1         10.09          139.3 s     17.7 s
	!        2         10.01          229.9 s     25.2 s
	!        8         10.01          259.9 s     26.8 s
	!
	! THE SELF-CONSISTENT FIELD DOES NOT PAY, and the reason is that the
	! outer pass count is not set by this mode.  Closing the cell's own lag
	! moves that count by 0.8 percent and costs 1.65x; on mol_base_handoff
	! (1500 steps) it moves it the wrong way, 7.19 -> 7.24, and costs 2.0x
	! (273.5 s -> 548.6 s).
	!
	! AND IT IS NOT THE MODE THE CARRIED-COLUMN SCHEME LEAVES EITHER, which
	! is what this measurement was expected to show.  COST3 cut the pass
	! count 10.10 -> 7.13 by carrying the star-ward columns forward
	! (xuv_field_block_cells = 8) and read what remained as the cell's own
	! optical depth.  MEASURED with one binary carrying BOTH levers
	! (block 8, four self-field passes): 7.20 passes, 190.1 s, against
	! COST3's 7.13 and 103.4 s for the columns alone.  With both of the
	! stellar field's lags closed, seven passes remain, so what sets them is
	! not the field: it is the temperature/composition alternation of the
	! coupled step and the two couplings documented as lagged one sweep --
	! the H(n=2) Balmer continuum and the locally absorbed He recombination
	! radiation.  Which of those dominates is not measured here.  Hence one.
	!
	! The loop does NOT narrow the parallel front, unlike the other lever:
	! a cell reads only the star-ward columns of the cells outside it, which
	! are fixed for the whole block, and its own iterate, so the sweep stays
	! one parallel region over the grid at any pass count.  MEASURED on
	! wasp_full 300 at two passes, one thread against sixteen: byte-identical
	! output, as the default scheme is on the same case.
	integer, save :: xuv_self_field_passes = 1
	! One cell's accepted equilibrium state, reported by every sweep that
	! touches it: the acceptance class, the normalized reaction residual the
	! state was accepted on, the solver exit code and the hydrogen partition.
	! A caller asking whether the elimination pins a composition, or only
	! brackets it within the acceptance tolerance, needs the residual of the
	! accepted state and not the sweep's maximum. Cell index from
	! EXHALE_IEQ_REPORT_CELL. The default has to lie outside 1-Ng..N+Ng --
	! the ghosts carry negative indices and 0 is one of them -- so it is the
	! most negative integer and no cell can equal it.
	integer, parameter :: ieq_report_no_cell = -huge(1)
	integer, save :: ieq_report_cell = ieq_report_no_cell
	! THE CELLS THE CHEMICAL-DECAY DIAGNOSTIC MEASURES AT, and what it
	! measured there. The reaction Jacobian of a cell can only be formed
	! where that cell's own rate coefficients are live, which is inside the
	! sweep's loop, so the cells are named before the sweep
	! (set_molecular_decay_cells) and the rates are read after it
	! (molecular_decay_rate_of_the_cell). Untouched by a run that names
	! none, which is every run by default.
	integer, parameter :: ieq_decay_ncell_max = 8
	! The four molecular rows of the network, in the layout order of
	! System_HeH_mol: H2, H2+, H3+, HeH+.
	integer, parameter :: n_mol_decay_row = 4
	integer, save :: ieq_decay_cells(ieq_decay_ncell_max) = 0
	integer, save :: ieq_decay_ncell = 0
	real*8,  save :: ieq_decay_row(ieq_decay_ncell_max,n_mol_decay_row) = 0.0d0
	real*8,  save :: ieq_decay_slow(ieq_decay_ncell_max) = 0.0d0
	real*8,  save :: ieq_decay_fast(ieq_decay_ncell_max) = 0.0d0
	logical, save :: ieq_decay_have(ieq_decay_ncell_max) = .false.
	! A cell leaves the loop above early once its accepted stage fractions
	! stop moving by more than this.  A decade tighter than the coupled
	! loop's own composition tolerance (csm_comp_tol = 1e-8, EXHALE_main),
	! so the inner loop is never what stops the outer one short.
	real*8,  parameter :: xuv_self_field_tol    = 1.0d-9

	! ---- THE BASE HANDOFF'S OWN CLOSURE, IN A LOWER GHOST CELL ----
	!
	! The handoff states the molecular partition of the NON-IONIZED
	! hydrogen, so what it imposes on a ghost is x_H2 = x2 (1 - x_ion), and
	! x_ion is the ionization the wind's field produces in that same ghost:
	! the two are one system and neither is data for the other. Closed here
	! by the cell's own passes -- impose, solve the ghost's ionization
	! balance against it, impose again at the ionization the solve returned
	! -- and NOT by alternation across sweeps, which is what left the ghost
	! carrying an iteration history of the composition the state was entered
	! with (9.06 rounding floors of the base continuity row on the hot-Uranus
	! molecular state, docs/lhs1140b_stationary_D5a_20260918.md section 8).
	!
	! THE TOLERANCE IS ON THE STATED EQUATION and is absolute in x_H2, a
	! fraction of the cell's hydrogen nuclei: the pair is closed when
	! |x_H2 - x2 (1 - x_ion)| at the state the pass returned is below this.
	! The passes reach the root by a secant step on that scalar equation
	! (the ghost reservoir block of the sweep). What the solve underneath
	! resolves is MEASURED on the LHS 1140 b molecular seed ghosts: the
	! secant takes the residual to 2e-16 to 6e-16 in 10 or 11 passes, so the
	! tolerance asks nothing the ionization solve cannot deliver.
	!
	! A FAILURE IS A REFUSED BOUNDARY. Exhausting the passes without reaching
	! the tolerance means the pair has no fixed point the sweep can reach,
	! and the ghost composition would then be an arbitrary iterate; the run
	! stops and says so rather than handing the boundary a state nothing
	! stands behind.
	real*8,  parameter :: base_ghost_closure_tol    = 1.0d-10
	integer, parameter :: base_ghost_closure_passes = 30
	! What the closure reached in this sweep: the largest residual of
	! x_H2 - x2 (1 - x_ion) left in a lower ghost, and the most passes any
	! of them took. Reduced over the cell loop and handed to base_boundary,
	! which reports them with the reservoir.
	real*8,  save :: ghost_closure_res_sweep  = 0.0d0
	integer, save :: ghost_closure_pass_sweep = 0

	! ---- THE GHOST COMPOSITION IS A FIXED POINT OF THE SOLVE THAT RETURNS IT ----
	!
	! THE INVARIANT. The composition the two lower ghost cells leave a sweep
	! with is the composition that same sweep returns when it is entered at
	! it. The lower boundary is built on the ghost, and a ghost that still
	! remembers what it was seeded with is a boundary that is a function of
	! the entry state as well as of the physical column and the declared
	! inputs.
	!
	! WHY IT NEEDS ITS OWN STATEMENT. Every measure the sweep already applies
	! is satisfied at ghosts 5 per cent apart in their trace ions: on one
	! frozen physical column, reservoir, spectrum and executable, two entry
	! compositions return ghosts 4.7e-3 (H II, He II) to 5.5e-2 (H2+) apart,
	! both accepted at reaction residuals of about 3e-26, both with the
	! imposed H2 partition closed to 0.0, both with hybrd1 answering info = 1,
	! and the difference does not close over six decades of the inner stopping
	! accuracy (READ, docs/lhs1140b_p1_step1b_20260919.md sections 1 and 3).
	! The acceptance tests are blind in the direction the base mass row reads.
	!
	! WHAT IS COMPARED, AND IN WHICH UNITS. Each species of the two lower
	! ghost cells, in the f_sp layout, which is the species' number density
	! per unit mass of the cell: that is the weight with which the species
	! enters the cell's particle count, and, times its own mass, the cell's
	! mass and the base face flux built on it. The move is therefore taken in
	! the composition's own units and NOT as each species against its own
	! value: no quantity the lower boundary builds -- the mass flux, the
	! charge, the particle count, the opacity, the caloric energy -- reads a
	! stage at 1e-10 of the cell with the weight it reads one carrying a
	! third of it, and the composition solve underneath pins the ghost
	! ABSOLUTELY and not stage by stage (below). The ghost's heavy-particle
	! plus electron count is the second leg, taken against its own value
	! because it is an O(1) quantity: that sum is what turns the ghost's
	! temperature into its pressure and its caloric energy, so the contract
	! is not closed until the rates the ghost was solved with belong to the
	! ghost's own count.
	!
	! THE ACCURACY, AND WHY THIS NUMBER. It is bracketed from above by what
	! the base continuity row can resolve and from below by what the
	! composition solve delivers, and one decade satisfies both.
	!   Above: that row is a near-cancellation of two face fluxes, so it reads
	! the ghost with a large gain. MEASURED on the LHS 1140 b hot-Uranus
	! molecular state, a ghost that moves by 4.1e-9 in its hydrogen partition
	! (x_HII 8.848e-7 against 8.888e-7, x_H2 by 3.9e-9) takes that row from
	! 3.90 to 111.3 of its own rounding floors, i.e. 107.4 floors per 4.1e-9
	! of composition, so ONE rounding floor of the row is bought by 3.8e-11
	! (READ, docs/lhs1140b_p1_step1_20260919.md sections 5.1 and 5.2; the
	! f_sp layout and that partition differ by the hydrogen nuclei per unit
	! mass, 0.75 there, which is an O(1) factor). A move of 1e-11 is a quarter
	! of one rounding floor of the row that reads it.
	!   Below: the composition solve stops at sqrt(dpmpar(1)) = 1.49e-8 and
	! pins the returned ghost to a plateau in the composition's own units.
	! MEASURED on examples/15_molecular with the carriers transported, where
	! the map's gain in the H3+ direction is close to one, the sweep settles
	! into an exact two-cycle whose largest move is 6.5e-12 (H I) and whose
	! H3+ leg is 2.2e-12; tightening the inner stop to 1e-12 takes the whole
	! sweep under this accuracy in two applications
	! (docs/lhs1140b_p1_step2b_20260920.md).
	!   A relative statement per species cannot be made at all there: 2.2e-12
	! of composition is 5.1e-3 OF H3+'s own value at 4.2e-10 of the cell, and
	! 1.5e-2 of H2+'s at 4.7e-18, while the same move is 0.06 of one rounding
	! floor of the base row.
	!
	! THE PASS BOUND. The furthest seed the boundary can present is the
	! reservoir row, 9.7e-1 from the ghost in H II, and the map's MEASURED
	! gain from seed to returned ghost is 4.9e-3 (the same memo, section 6.1),
	! so four applications carry that seed to 5e-10. Twenty passes is five
	! times what the measured contraction needs from the worst seed.
	!
	! A BOUNDARY THAT CANNOT REACH IT FAILS LOUDLY, with the numbers, exactly
	! as the ghost's H2 partition closure does: a ghost that is not a fixed
	! point of its own solve is an iterate, and nothing stands behind a
	! boundary built on one.
	!
	! VALIDITY: the lower ghost cells of a run whose network is molecular and
	! carries helium, where the boundary's composition holds the trace ions
	! the base row reads and where every measurement above was made. The
	! atomic ghost is not under the contract and its seed dependence is not
	! measured.
	real*8,  parameter :: ghost_composition_fixed_point_move   = 1.0d-11
	real*8,  parameter :: ghost_count_fixed_point_move         = 1.0d-6
	integer, parameter :: ghost_composition_fixed_point_passes = 20
	! What the last sweep reached: the largest move of a ghost species
	! density over the last application of the map, in the composition's own
	! units, the relative move of the ghost's particle plus electron count
	! over the same application, and how many applications the sweep took.
	! Handed to base_boundary for the boundary report.
	real*8,  save :: ghost_fixed_point_move_sweep    = 0.0d0
	real*8,  save :: ghost_thermal_move_sweep        = 0.0d0
	integer, save :: ghost_fixed_point_passes_sweep  = 0
	! Whether this application of the composition map is the one that
	! advances the non-root streak, which counts consecutive SWEEPS and not
	! applications. True on the first application of every sweep.
	logical, save :: sweep_advances_the_streak = .true.

	! molecular species densities (cols 1 H2, 2 H2+, 3 H3+, 4 HeH+;
	! zero unless thereis_mol).  Module state: written by the equilibrium
	! solve, read by write_output for the extra output columns.
	real*8, dimension(:,:), allocatable :: nmol_eq

	! Lyman-Werner photodissociation diagnostics, filled by the equilibrium
	! solve when a band flux is supplied and read by write_lyman_werner:
	! star-ward H2 column [cm^-2], the self-shielding factor the rate
	! carries (Richings, Schaye & Oppenheimer 2014; lyman_werner.f90
	! sec. 2), and the resulting dissociation rate [s^-1]. Untouched
	! (shielding 1, rate 0) when the run supplies no band flux.
	real*8, dimension(:), allocatable :: NH2_col_lw
	real*8, dimension(:), allocatable :: f_shield_lw
	real*8, dimension(:), allocatable :: k_lw_diss
	! The two branching ratios of one Lyman-Werner absorption at each cell,
	! from the same level-resolved table as k_lw_diss (utils_ion_eq
	! fuv_lw_photon_field says what each one counts).  p_lw_single is the
	! fluorescence heat return, p_lw_absorbed the band photon ledger.  Both
	! zero when the run supplies no band flux.
	real*8, dimension(:), allocatable :: p_lw_single, p_lw_absorbed
	! Transmission of the Lyman-Werner LINES down to each cell's inner face,
	! 1 - A: the fraction of the 912-1110 A band that H2 has NOT already
	! taken, which the H2O and OH continua of the same interval multiply
	! their own attenuation by (water_photolysis.f90 sec. 3). Identically 1
	! without a Lyman-Werner band flux. a_lines_lw is the largest A the
	! column reached, i.e. the deepest cell's removed fraction; the run
	! reports once if the H2 column leaves the checked range of either of
	! the two self-shielding fits (lyman_werner.f90 sec. 2).
	real*8, dimension(:), allocatable :: tr_lines_lw
	! H2 photoionization rate of the last sweep [s^-1], kept for the H2
	! budget decomposition write_output prints: it is the one loss channel
	! of the H2 row that cannot be rebuilt from the output files.
	real*8, dimension(:), allocatable :: P_H2_eq
	real*8 :: a_lines_lw = 0.0d0
	logical :: warned_lw_fit_range = .false.
	! NH2_db96_max: Draine & Bertoldi (1996) sec. 5.2, the agreement of their
	! eq. (37) with an exact calculation is excellent "even out to the largest
	! column densities considered, N2 = 3 x 10^21 cm^-2".  They state no upper
	! limit on the fit itself, so a deeper column is undemonstrated rather
	! than outside a published range.  This fit sets the band share A only;
	! the photodissociation RATE is on the level-resolved table.
	real*8, parameter :: NH2_db96_max     = 3.0d21

	! Largest ratio, over the cells of the last equilibrium field solve, of
	! the star-ward H2 column to the column at which the self-shielding table
	! stops being a value and becomes an upper bound (lyman_werner.f90
	! sec. 2d).  Written by fuv_lw_photon_field, read by the warnings here and
	! in write_output.  Zero for a run with no Lyman-Werner field.
	real*8 :: lw_col_over_overlap = 0.0d0

	! Oxygen-chemistry state (the A2 option; all zero unless
	! thereis_oxychem).  Written by the equilibrium solve, read by
	! write_output for output/Ion_species.txt, output/Oxygen_chemistry.txt
	! and output/FUV_bands.txt, and by the heat breakdown.
	!   nox_eq   cols 1 OH, 2 H2O, 3 CO [cm^-3]
	!   n_o1d_eq O(1D) density from its local steady state [cm^-3]
	!   NH2O_col, NOH_col   star-ward columns of the two FUV absorbers
	!   j_h2o_fuv, j_oh_fuv band-resolved photodissociation rates [s^-1]
	!   tau_fuv             band-resolved star-ward optical depth
	!   heat_fuv            FUV photolysis heating [erg cm^-3 s^-1]
	real*8, dimension(:,:), allocatable :: nox_eq
	real*8, dimension(:),   allocatable :: n_o1d_eq
	real*8, dimension(:),   allocatable :: NH2O_col, NOH_col
	!   NCO_col             star-ward CO column [cm^-2]
	!   k_co_diss           CO photodissociation rate, the cell mean [s^-1]
	!   theta_co_shield     the Visser shielding function of that cell's own
	!                       two columns; 1 where nothing shields
	real*8, dimension(:),   allocatable :: NCO_col, k_co_diss
	real*8, dimension(:),   allocatable :: theta_co_shield
	real*8, dimension(:,:), allocatable :: j_h2o_fuv, j_oh_fuv, tau_fuv
	real*8, dimension(:),   allocatable :: heat_fuv
	!   heat_chem           chemical heat of the collisional molecular
	!                       reactions [erg cm^-3 s^-1]
	real*8, dimension(:),   allocatable :: heat_chem

	! The coefficient state each cell was left in by the last equilibrium
	! sweep, kept so the carrier transport (diffusive_photochemistry) can
	! evaluate the SAME chemistry rows this module solves without rebuilding
	! the photoionization, recombination and collisional rates. It is a copy
	! of the ieq_cell each cell was solved with, which is why a field added
	! to that type reaches the transport operator with no second edit here.
	! bg_ready is false until the first sweep has filled it: a transport step
	! before that would run on zeros.
	type(ion_rates), dimension(:), allocatable :: bg_cell
	logical :: bg_ready = .false.

	! THE CELL STATE OF THE SWEEP, KEPT SO THE ELIMINATED-SPECIES CLOSURE
	! CAN BE MEASURED WITHOUT SOLVING IT AGAIN.
	!
	! The closure of B1a section 2.4 -- the neutral stage of every element
	! and the electron density eliminated by the element totals and by
	! charge neutrality -- is a residual of the SAME system rows the cell
	! solve used, and until now it existed only inside the sweep, which also
	! rewrites f_sp and rho. A certification of a state cannot run that
	! sweep: it would change the state it is judging.
	!
	! What a residual of one cell needs is the cell's rate state and the
	! electron density the turnover scale is referred to. Both are recorded
	! here as the sweep passes each cell, so that
	! ionization_closure_residual_profile can reinstall a cell and call
	! normalized_reaction_residual on it. Nothing here feeds the solve; the
	! arrays are written and never read by ioniz_eq itself.
	!
	! The composition is NOT recorded with them: the evaluator reads it
	! from the caller's f_sp, so that a state whose composition has moved
	! since the sweep is measured as it stands rather than as the sweep
	! left it. What IS checked is that the element totals of the state are
	! the ones these rates were built with; a state whose totals have moved
	! is a different gas and the closure is unavailable for it.
	type(ion_rates), dimension(:), allocatable :: ieq_rate_cell
	real*8,  dimension(:),   allocatable :: ieq_ne_cell
	real*8,  dimension(:),   allocatable :: ieq_TK_cell
	real*8,  dimension(:),   allocatable :: ieq_ntot_cell
	! Metal coefficients of each cell, in the order set_metal_coeffs takes
	! them: ntot, g0, g1, b0, b1, a1, a2.
	real*8,  dimension(:,:,:), allocatable :: ieq_met_coef
	integer, parameter :: ieq_n_met_coef = 7
	logical :: ieq_rates_ready   = .false.
	integer :: ieq_neq_stored    = 0
	integer :: ieq_mbase_stored  = 0
	integer :: ieq_iox_stored    = 0

	! THE RATE STATE OF EVERY CELL AS THE LAST SWEEP LEFT IT, held aside
	! and put back. The arrays above are what the closure evaluators
	! measure a composition against, so the composition and the rates
	! beside it are ONE state: an operation that is undone has to put the
	! rates back together with the composition, or the closure of the
	! reinstated composition is read against the rates of the state that
	! was discarded (MEASURED on the carrier_retry column: 5.3e-4 against
	! 1.6e-16, twelve decades, docs/lhs1140b_stationary_D4a_20260918.md
	! table 2). The enumeration stands here, beside the declarations, so
	! that an array added to the rate state is added to the snapshot in
	! the same place. Allocation status is part of the state: an array
	! that was not allocated when the snapshot was taken is deallocated
	! again.
	type :: ieq_rate_state
		type(ion_rates), allocatable :: rate_cell(:)
		real*8,  allocatable :: ne_cell(:), TK_cell(:), ntot_cell(:)
		real*8,  allocatable :: met_coef(:,:,:)
		integer, allocatable :: nonroot_streak(:)
		logical :: rates_ready  = .false.
		integer :: neq_stored   = 0
		integer :: mbase_stored = 0
		integer :: iox_stored   = 0
	end type ieq_rate_state

	! The scratch of one cell that the closure evaluator has to install,
	! held aside and put back around every measurement (see the
	! evaluator's header).
	type :: ieq_closure_scratch
		type(ion_rates) :: cell
		integer :: nelem = 0
		logical :: met_alloc = .false.
		real*8  :: ntot(n_melem) = 0.0d0
		real*8  :: g0(n_melem)   = 0.0d0
		real*8  :: g1(n_melem)   = 0.0d0
		real*8  :: b0(n_melem)   = 0.0d0
		real*8  :: b1(n_melem)   = 0.0d0
		real*8  :: a1(n_melem)   = 0.0d0
		real*8  :: a2(n_melem)   = 0.0d0
		integer :: top(n_melem)  = 0
		! The molecular and oxygen rate coefficients of System_HeH_mol and
		! the reciprocal row turnover scales, in one packed vector.
		real*8  :: molc(27) = 0.0d0
		real*8  :: minv(n_mol_rows_max) = 1.0d0
	end type ieq_closure_scratch

	! THE FROZEN BACKGROUND OF A STATE THE STEADY SOLVER HOLDS.
	!
	! bg_cell is rewritten by EVERY molecular sweep, and the steady solver
	! sweeps on states it does not keep: Jacobian columns, Krylov products and
	! line-search trials (docs/Update_EXHALE_stage1.md section 121). When the solver
	! returns, the last sweep it ran was in general a trial the line search
	! REJECTED, and the best-iterate restore can additionally hand back a state
	! from several outer iterations earlier -- so the background left in
	! bg_cell need not be the background of the state that is handed back.
	! Everything else the solver returns is consistent: f_sp is restored with
	! the best iterate, and nmol_eq / nox_eq are rebuilt from f_sp at the top
	! of every sweep.
	!
	! These two copies close that gap. bg_cell_adopted is refreshed whenever
	! the solver evaluates a state it holds -- its current iterate, or a trial
	! the line search has just accepted -- and bg_cell_best follows the
	! best-iterate bookkeeping; the solver installs one of them before
	! returning. The invariant is that every cell state derived from the
	! returned state was evaluated AT the returned state.
	type(ion_rates), dimension(:), allocatable, save :: bg_cell_adopted
	type(ion_rates), dimension(:), allocatable, save :: bg_cell_best

	! WHICH STATE THE SWEEP IS EVALUATING, and the three ledgers that keep
	! them apart (docs/Update_EXHALE_stage1.md section 121).
	!
	! The marching loop calls ioniz_eq on the state the run holds. The steady
	! (JFNK / PTC) solver calls it on three different things: the current
	! iterate, whose equilibrium the solver's residual is built from; the
	! finite-difference and Krylov probes Y + eps*v, which are directional
	! derivative samples and carry no physical meaning of their own; and the
	! line-search trials, which are candidates until one is adopted. A single
	! set of run-wide counters records all of them as if the run had accepted
	! them: measured on backup/regression/heh_1_newton, the marching phase
	! reported not one non-root acceptance in 2002 steps and the JFNK phase
	! then added 18,339, largest residual 1.7e+97, none of which the run ever
	! adopted. The same mixing gave ieq_nonroot_streak -- whose limit stops
	! the run -- a meaning it does not have: probe sweeps both raised and
	! cleared it, so "consecutive sweeps on a non-root" counted evaluations
	! the run threw away.
	!
	! Each sweep is therefore tagged, the totals go to the ledger of that tag,
	! and the end-of-run report prints the three separately. Only the two
	! ledgers of states the run holds carry the non-root streak and its stop.
	integer, parameter :: ieq_state_marching          = 1
	integer, parameter :: ieq_state_steady_iterate    = 2
	integer, parameter :: ieq_state_steady_candidate  = 3
	integer, save :: ieq_sweep_state_kind = ieq_state_marching
	! EXHALE_MASS_PROJECTION, read once (mass_projection_of_the_sweep)
	logical, save :: mass_projection_on    = .true.
	logical, save :: mass_projection_known = .false.
	! The two thresholds of the projection (see the block that applies it)
	! and what it has done over the run: the largest |s - 1| seen and the
	! number of sweeps at which it was refused.  Read by the run's own
	! reports; nothing branches on them.
	real*8, parameter :: ieq_mass_projection_report  = 1.0d-12
	real*8, parameter :: ieq_mass_projection_refuse  = 1.0d-8
	real*8, save      :: ieq_mass_projection_worst   = 0.0d0
	integer, save     :: n_ieq_mass_projection_refused = 0
	public :: ieq_mass_projection_worst, n_ieq_mass_projection_refused
	public :: ieq_mass_projection_report, ieq_mass_projection_refuse

	! One ledger per sweep tag. The fields are the run-wide totals that used
	! to be separate module variables:
	!   n_sweep        equilibrium sweeps recorded in this ledger
	!   n_reseed       stored states rejected as a starting point
	!   n_retry        cell solves that needed a second or third starting point
	!   n_unphys       first roots that lay outside the physical simplex
	!   n_noroot       cells where no starting point produced an admissible
	!                  root (atomic branch: the uncoupled balance was handed
	!                  back)
	!   n_mol_clamped  molecular cells whose roots all left the simplex, so the
	!                  closest one was clamped onto the element budget
	!   n_mol_info     hybrd1 exit-code histogram of the molecular cell solves
	!                  (0 = improper input or iflag < 0, 1 = converged to tol,
	!                  2 = iteration limit, 3 = xtol too small, 4/5 = no
	!                  progress); every attempt of every molecular cell counts
	!   acc_n          acceptances by class 1-5 (see ieq_res_tol below)
	!   acc_resmax     largest normalized reaction residual of each class
	!   hist_conv      residual-decade histogram of the physical iterates that
	!   hist_uncv      the solver did / did not report as converged
	!   n_cce_*        constrained element-conserving continuation solve:
	!   cce_seconds    cells promoted, candidates accepted as roots, hybrd1
	!                  field solves spent, wall time
	!   n_nonfinite    cells whose accepted state has a normalized reaction
	!                  residual that is not a finite number -- the state broke
	!                  a balance row of the network and cannot be described
	!   n_offsimplex   cells at which NO starting point stayed inside the
	!                  element simplex (n_mol_clamped + the atomic handbacks)
	!   viol_worst     largest element-budget violation of any iterate that
	!                  left the element simplex, whether or not a later
	!                  starting point of the same cell then found a root
	!   streak_peak    longest run of CONSECUTIVE sweeps any one cell spent on
	!                  a non-root acceptance. Zero in the probe ledger by
	!                  construction: a probe neither raises the streak nor
	!                  clears it, so "consecutive" has no meaning there.
	! The last three are the admissibility signal the steady solver reads back
	! from a sweep; they are accumulated in every ledger, and it is the
	! candidate ledger's copy that decides whether a trial or probe state is
	! usable at all.
	type :: ioniz_eq_ledger
		integer :: n_sweep       = 0
		integer :: n_reseed      = 0
		integer :: n_retry       = 0
		integer :: n_unphys      = 0
		integer :: n_noroot      = 0
		integer :: n_mol_clamped = 0
		integer :: n_mol_info(0:5) = 0
		integer :: acc_n(6)      = 0
		real*8  :: acc_resmax(6) = 0.0d0
		integer :: hist_conv(0:15) = 0
		integer :: hist_uncv(0:15) = 0
		integer :: n_cce_attempt = 0
		integer :: n_cce_root    = 0
		integer :: n_cce_solve   = 0
		real*8  :: cce_seconds   = 0.0d0
		integer :: n_nonfinite   = 0
		integer :: n_offsimplex  = 0
		real*8  :: viol_worst    = 0.0d0
		integer :: streak_peak   = 0
	end type ioniz_eq_ledger

	! The marching ledger is kept once per LEDGER FAMILY
	! (docs/a0_run_mode_contract_20260906.md section 5): index
	! ledger_family_init holds the sweeps taken while the run was reaching a
	! state, index ledger_family_phys those taken inside accepted physical
	! steps. The fields are the same; what differs is what a number in them
	! means, and summing the two would produce an acceptance census belonging
	! to no single run state. An initialization run records everything in the
	! first and leaves the second empty, which is what it was before the
	! split.
	type(ioniz_eq_ledger), save :: ieq_marching_ledger(2)
	type(ioniz_eq_ledger), save :: ieq_steady_iterate_ledger
	type(ioniz_eq_ledger), save :: ieq_steady_candidate_ledger

	! The acceptance census of the LAST sweep, whatever its tag. The
	! elimination inside a steady residual evaluation is a Picard iteration
	! whose fixed point is the composition c*(Y) the residual R(Y) = L(Y) +
	! S(Y,c*(Y)) is defined at, and a cell accepted under the relaxation
	! amnesty (class 4) or above its cap (class 6, the entry composition
	! kept) is a cell that has no root at this state: the elimination then
	! carries the seed and the residual is not a function of Y. A caller that
	! has to know that about one sweep reads it here.
	type(ioniz_eq_ledger), save :: ieq_sweep_ledger_last

	! Acceptance tolerance on the normalized reaction residual of an
	! equilibrium state (docs/supersonic_molecular_base.md section 11.5-A;
	! docs/Update_EXHALE_stage1.md section 113). A state is accepted as a ROOT of
	! the ionization/chemical network only when it lies inside the element
	! bounds AND the largest reaction imbalance of any row, in units of the
	! row's turnover rate (normalized_reaction_residual), is at or below
	! this value; the solver exit code is neither sufficient nor necessary
	! (MINPACK info = 1 is an xtol statement about the step, and info = 4
	! routinely returns finished roots it cannot certify).
	!
	! The value is set from the measured populations of the full 8-case
	! regression matrix plus the He/H = 1 molecular case (2026-08-31):
	! across the matrix every solver-converged accepted root sits at
	! res <= 7.8e-8 and every root accepted without solver convergence at
	! res <= 2.7e-7 -- the scale is hybrd1's xtol = sqrt(eps) ~ 1.5e-8
	! times a scaled-Jacobian norm of order 1-10 -- while the smallest
	! above-tolerance acceptance of the matrix (a cold-start base cell
	! recovering over the next sweeps) is 1.1e-3 and the poisoned events
	! of the He/H = 1 case reach 2e-2..1.4e6. 1e-6 sits 3.7x above the
	! measured root tail, in a decade ([1e-6, 1e-5)) that is empty in
	! every matrix histogram, and 1100x below the matrix's smallest
	! non-root; the He/H = 1 case additionally carries one borderline
	! accepted iterate at 2.1e-6, which this value classifies (and marks)
	! as a non-root rather than stretching the root band to cover it.
	real*8, parameter :: ieq_res_tol = 1.0d-6

	! THE LOWER GHOST CELLS' OWN ACCEPTANCE, cell by cell, for the ghost
	! record (boundary_state_trace.f90). A measurement: nothing in the
	! solve, the boundary or the energy update reads any of it. The
	! reaction residual is the one the cell was ACCEPTED at, the partition
	! residual is that of x_H2 - q_H2,base (1 - x_ion) at the returned
	! state, and the pass count is the pass the closure ended on.
	real*8,  save :: ghost_acc_res_cell(1-Ng:0)      = 0.0d0
	real*8,  save :: ghost_closure_res_cell(1-Ng:0)  = 0.0d0
	integer, save :: ghost_closure_pass_cell(1-Ng:0) = 0

	! Said once for the run, where the ghost seed key is set, so that a run
	! without the key is the run without the code.
	logical, save :: ghost_seed_announced = .false.

	! THE COMPOSITION OF THE CELLS ABOVE THE LOWER GHOSTS, HELD AT THE ONE
	! THE SWEEP WAS GIVEN (EXHALE_INTERIOR_COMPOSITION_HELD, default off).
	! A MEASUREMENT: with the key set the sweep returns the composition it
	! was handed on every cell above the two lower ghosts and the
	! composition it solved on the ghosts themselves, so that a base row
	! measured after the call carries the ghost's own refresh and nothing
	! else. VALIDITY: a diagnostic of the base row, and of nothing else.
	! The heating and cooling the sweep returns are those of the
	! composition it SOLVED on every cell, so a run with the key set is not
	! a state anything may be certified on. With the key unset not one
	! number moves: the branch is not entered.
	logical, save :: interior_held_on   = .false.
	logical, save :: interior_held_read = .false.
	real*8, allocatable, save :: f_sp_interior_entry(:,:)

	! THE STOPPING TOLERANCE OF THE COMPOSITION'S INNER SOLVE, and the hook
	! that replaces it for a whole run (EXHALE_IEQ_TOL=<x>, default off).
	! The value the sweep uses is sqrt(dpmpar(1)) = 1.49e-8, MINPACK's own
	! scale: hybrd1 reads it as xtol, a RELATIVE STEP criterion, and
	! newton_dense reads it as ftol, a criterion on the residual norm of the
	! cell's network (newton_solver.f90 lines 84 to 140), so on either route
	! the returned composition is pinned only to about that, and the number
	! of inner iterations changes with the state. The composition is
	! therefore a PIECEWISE map of the state and the residual assembled from
	! it carries a non-smoothness floor far above its own rounding, which is
	! what the finite-difference Jacobian action of the stationary solve
	! divides by its probe arc (item N31). The hook exists to measure which
	! tolerance that floor follows; a non-positive or unreadable value
	! leaves the sweep at sqrt(dpmpar(1)).
	real*8,  save :: ieq_inner_tol      = 0.0d0
	logical, save :: ieq_inner_tol_read = .false.

	! Relaxation amnesty and its limit. The cold-start relaxation of a
	! healthy run genuinely passes through non-root acceptances and
	! recovers: measured on the regression matrix, the metals-on molecular
	! gates accept res ~ 28 states at step 0 in the shielded base and res
	! ~ 2.6e-3 states in one or two cells for up to 53 CONSECUTIVE sweeps
	! (mol_ir_bands, cell 260) before the solve lands on roots for the
	! rest of the run. An unconditional stop on the first non-root would
	! therefore kill working configurations, and no severity threshold
	! separates them (the healthy step-0 spike, 28, exceeds the poisoned
	! He/H = 1 step-0 event, 2e-2). What does separate a recovering
	! transient from a solution RESTING on a non-root is persistence: a
	! steady or pseudo-steady state re-evaluates the same cell every sweep,
	! so a wind built on a non-root fails the same cell indefinitely. A
	! cell may therefore carry a non-root acceptance -- loudly counted and
	! reported, never silent -- for at most this many consecutive sweeps
	! (~19x the largest healthy streak measured); one sweep more and the
	! run stops with full diagnostics (nonroot_equilibrium_stop).
	! Persistence is not the only bound. It says nothing about how far ONE
	! accepted non-root may move the cell in a single sweep, which is what
	! ieq_nonroot_res_cap below bounds; the streak counts both classes.
	integer, parameter :: ieq_nonroot_streak_stop = 1000
	! Consecutive sweeps each cell has spent on a non-root acceptance
	! (reset to zero by any accepted root), and the run-wide peak.
	integer, allocatable, save :: ieq_nonroot_streak(:)

	! THE SEVERITY BOUND OF THE AMNESTY. Persistence separates a recovering
	! transient from a solution resting on a non-root; it does not bound how
	! far one accepted non-root may move the cell in a single sweep. A
	! candidate whose normalized reaction residual is of order 1e4 is not a
	! nearly-solved state, and adopting it is a jump of decades in one sweep:
	! the measured event is oxygen_chemistry cell 7 at step 12, res 5.2e4
	! accepted, which took n_e from 1.1e3 to 3.5e10 cm^-3 and n(H+) from 0 to
	! 2.6e8 in one step and left the carrier system stiff. Above this cap the
	! candidate is therefore NOT adopted: the cell keeps the composition it
	! entered the sweep with (acceptance class 6), which is a state the run
	! already held, and the sweep goes on. It is not a stop, and the streak
	! that ends a run resting on a non-root counts class 6 exactly as it
	! counts class 4.
	! MEASURED on the sixteen default regression cases, 2026-09-06. Fifteen
	! of them accept NO non-root at all; every amnesty event of the matrix
	! belongs to the step-0 sweep of oxygen_chemistry, and its 115 events
	! fall in two populations with nearly four empty decades between them:
	!   res 1e-6 to 8.9   111 events, cells 1 to 185, the cold-start
	!                     transient (16 at 1e-6, 38 at 1e-5, 24 at 1e-4,
	!                     10 at 1e-3, 5 at 1e-2, 14 at 1e-1, 4 in [1,9))
	!   res 8.6e4 to 1.5e5  4 events, cells 85 to 88 -- and cells 85 to 88
	!                     are exactly the four the energy update then fails
	!                     on, asking for a temperature below the floor.
	! Nothing lies between 8.9 and 8.6e4. The cap is the smallest power of
	! ten at least 3x above the largest event of the first population
	! (3 x 8.9 = 27), which also clears the res ~28 amnesty events the
	! metals-on molecular gates were measured to accept in 2026-08 (those
	! gates accept no non-root at all today: the constrained continuation
	! reaches those cells). It sits 11x above the first population and 860x
	! below the second.
	real*8, parameter :: ieq_nonroot_res_cap = 1.0d2

	! The acceptance statistics of the He-branch equilibrium solves
	! (docs/Update_EXHALE_stage1.md section 113) live in the ledgers declared above:
	! how many cell states were accepted as (1) solver-converged roots,
	! (2) roots without solver convergence, (3) projected/handback states
	! whose rechecked residual still marks a root, (4) NON-ROOT states
	! accepted under the relaxation amnesty, (5) roots of the constrained
	! element-conserving continuation solve
	! (constrained_chemical_equilibrium), (6) NON-ROOT candidates above the
	! amnesty cap, for which the cell kept the composition it entered the
	! sweep with, with the largest normalized
	! reaction residual each class carried; the residual-decade histograms of
	! the physical iterates (converged and not), which locate the gap between
	! roots and non-roots the tolerance sits in; and the cost of the
	! constrained solve.

	! Print budget of the acceptance report lines: a pathological run must
	! not flood the log; the ledgers keep the full population. Only states
	! the run holds are reported line by line -- a steady-solver probe is a
	! directional-derivative sample, not an acceptance, and its non-root
	! events (18,339 in one measured JFNK phase) would bury the ones that
	! matter.
	integer, save :: ieq_acc_nprint = 0
	integer, parameter :: ieq_acc_nprint_max = 2000

	contains

	! THE STOPPING TOLERANCE THE COMPOSITION'S INNER SOLVE IS ASKED FOR.
	! Returns what EXHALE_IEQ_TOL names, or the argument when the variable
	! is unset, unreadable or non-positive. Pure of module state, so a
	! caller may cache the answer and a test may ask it twice with two
	! environments.
	real*8 function composition_solve_tolerance(tol_minpack)
	real*8, intent(in) :: tol_minpack
	character(len=32)  :: env_ieq
	real*8             :: tol_asked
	composition_solve_tolerance = tol_minpack
	call get_environment_variable('EXHALE_IEQ_TOL', env_ieq)
	if (len_trim(env_ieq) .eq. 0) return
	tol_asked = 0.0d0
	read(env_ieq,*,err=971,end=971) tol_asked
  971	continue
	if (tol_asked .gt. 0.0d0) composition_solve_tolerance = tol_asked
	end function composition_solve_tolerance


	! WHETHER THE CELLS ABOVE THE LOWER GHOSTS ARE HANDED BACK AS THEY CAME.
	! EXHALE_INTERIOR_COMPOSITION_HELD, read once for the run and announced
	! only when it is set, so that a run without the key is the run without
	! the code.
	logical function interior_composition_is_held()
	character(len=32) :: env_held
	if (.not. interior_held_read) then
		call get_environment_variable('EXHALE_INTERIOR_COMPOSITION_HELD',   &
		                              env_held)
		interior_held_on   = (len_trim(env_held) .gt. 0)
		interior_held_read = .true.
		if (interior_held_on) write(*,'(A)') ' (ioniz_eq) EXHALE_'//        &
		     'INTERIOR_COMPOSITION_HELD: the sweep returns the composition'//&
		     ' it was given on every cell above the two lower ghosts'
	endif
	interior_composition_is_held = interior_held_on
	end function interior_composition_is_held


	subroutine ioniz_eq_allocate_arrays
	! Allocate the grid-sized module arrays once the number of cells N is
	! known; called from EXHALE_main right after input_read. The values are
	! the initializers the declarations used to carry.

	allocate(nmol_eq(1-Ng:N+Ng,4))
	! bg_cell and ieq_rate_cell below come out of allocation as the zero
	! rate state -- no gas, no rates, no imposed carrier partition -- by the
	! default initializers of ion_rates: a cell the sweep has not reached
	! yet carries no rates, and the checkpoint that hashes, saves and
	! restores these arrays from the first step on therefore hashes a number
	! of the state and not of whatever block the allocator recycled.
	allocate(bg_cell(1-Ng:N+Ng))
	allocate(ieq_nonroot_streak(1-Ng:N+Ng))
	ieq_nonroot_streak = 0
	allocate(NH2_col_lw(1-Ng:N+Ng), f_shield_lw(1-Ng:N+Ng),               &
	         k_lw_diss(1-Ng:N+Ng), tr_lines_lw(1-Ng:N+Ng),                &
	         p_lw_single(1-Ng:N+Ng), p_lw_absorbed(1-Ng:N+Ng),            &
	         P_H2_eq(1-Ng:N+Ng))
	! The heating channels of the state each sweep returns; utils_ion_eq
	! owns the array and the name list, the sweep fills it.
	allocate(heat_channel_state(1-Ng:N+Ng,n_heat_channel))
	heat_channel_state = 0.0d0
	allocate(nox_eq(1-Ng:N+Ng,3), n_o1d_eq(1-Ng:N+Ng),                    &
	         NH2O_col(1-Ng:N+Ng), NOH_col(1-Ng:N+Ng),                     &
	         j_h2o_fuv(1-Ng:N+Ng,n_fuv_band),                             &
	         j_oh_fuv(1-Ng:N+Ng,n_fuv_band),                              &
	         tau_fuv(1-Ng:N+Ng,n_fuv_band), heat_fuv(1-Ng:N+Ng),           &
	         heat_chem(1-Ng:N+Ng))
	allocate(NCO_col(1-Ng:N+Ng), k_co_diss(1-Ng:N+Ng),                    &
	         theta_co_shield(1-Ng:N+Ng))
	NCO_col         = 0.0d0
	k_co_diss       = 0.0d0
	theta_co_shield = 1.0d0

	nmol_eq     = 0.0d0
	NH2_col_lw  = 0.0d0
	f_shield_lw = 1.0d0
	tr_lines_lw = 1.0d0
	k_lw_diss   = 0.0d0
	p_lw_single   = 0.0d0
	p_lw_absorbed = 0.0d0
	P_H2_eq     = 0.0d0
	nox_eq      = 0.0d0
	n_o1d_eq    = 0.0d0
	NH2O_col    = 0.0d0
	NOH_col     = 0.0d0
	j_h2o_fuv   = 0.0d0
	j_oh_fuv    = 0.0d0
	tau_fuv     = 0.0d0
	! The two heating channels the sweep copies out of heat_channel_state
	! (photolysis and collisional chemical heat): no gas has been swept
	! yet, so neither channel is depositing anything.
	heat_fuv    = 0.0d0
	heat_chem   = 0.0d0

	! The cell state the closure evaluator reads (see its declarations).
	! n_x_max of the sweep: the molecular layout (7) plus the triplet plus
	! the two oxygen carriers plus two unknowns per metal element.
	allocate(ieq_rate_cell(1-Ng:N+Ng), ieq_ne_cell(1-Ng:N+Ng),            &
	         ieq_TK_cell(1-Ng:N+Ng), ieq_ntot_cell(1-Ng:N+Ng),            &
	         ieq_met_coef(n_melem,ieq_n_met_coef,1-Ng:N+Ng))
	ieq_ne_cell        = 0.0d0
	ieq_TK_cell        = 0.0d0
	ieq_ntot_cell      = 0.0d0
	ieq_met_coef       = 0.0d0

	end subroutine ioniz_eq_allocate_arrays

	subroutine molecular_carrier_densities_from_state(rho, f_sp)
	! Number densities [cm^-3] of the molecular species and the oxygen
	! carriers implied by a composition (rho, f_sp), written into the module
	! arrays nmol_eq and nox_eq that the output routines read.
	!
	! WHY IT IS NEEDED. Those arrays are filled by the equilibrium sweep, so
	! before the first sweep of a run they are zero -- and write_output takes
	! the H2 / H2+ / H3+ / HeH+ / OH / H2O / CO columns from them. A state
	! written BEFORE any sweep therefore reports a molecular gas as atomic:
	! measured on the mol_base_handoff case reloaded with EXHALE_DUMP_IC=1,
	! all four molecular columns came out exactly zero against the 5.34e12
	! cm^-3 of H2 the restart file carried. The state that a diagnostic
	! writes has to be the state it was given, so the dump fills them from
	! f_sp first.
	!
	! It is the same product f_sp*rho*n0 the sweep's own extraction uses --
	! but the sweep does not refresh the arrays afterwards, and the state
	! moves on: after the energy step Apply_BC re-extrapolates the lower
	! ghost cells' rho, so at a write the atomic columns (f_sp*rho at the
	! write) and the molecular columns (the sweep's) describe the ghost at
	! two densities. Measured on the hot-Uranus gate while its base still
	! breathes: the ghost rows' H nuclei per unit mass were 3e-7 (40 steps)
	! to 3e-6 (1000 steps) above the interior's, He exact, and zero once the
	! base had settled -- and a restart of such a file failed the round-trip
	! identity in exactly those two rows. Every write of a state therefore
	! calls this first (docs/Update_EXHALE_stage1.md section 169).
	real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp

	if (thereis_mol) then
		nmol_eq(:,1) = f_sp(:,isp_H2)  *rho*n0
		nmol_eq(:,2) = f_sp(:,isp_H2p) *rho*n0
		nmol_eq(:,3) = f_sp(:,isp_H3p) *rho*n0
		nmol_eq(:,4) = f_sp(:,isp_HeHp)*rho*n0
	endif
	if (thereis_oxychem) then
		nox_eq(:,1) = f_sp(:,isp_OH) *rho*n0
		nox_eq(:,2) = f_sp(:,isp_H2O)*rho*n0
		nox_eq(:,3) = f_sp(:,isp_CO) *rho*n0
	endif

	end subroutine molecular_carrier_densities_from_state

	!----------------------------------!

	subroutine ioniz_eq(T_in,n_io,f_sp_io,heat_out,cool_out,q,sweep_ledger,&
	                    rho_recon)
      	 		  
	integer :: j,im
	logical :: usednt                     ! Task 2: Newton-vs-fallback flag
	! Electron count of the lower ghost at the composition this sweep
	! returns: a measurement handed to base_boundary, not a boundary input.
	real*8  :: ghost_ne_solved
	! The ghost record's own scratch: the charge the written composition
	! carries, its distance from the electron count, and the ghost's
	! elemental ratios against the first physical cell's.
	integer :: jg
	real*8  :: f_ghost_seed(1-Ng:0,n_species)
	character(len=256) :: ghost_seed_used
	real*8  :: q_written, chg_gap, elem_gap, ratio_g, ratio_1
	real*8  :: xg_h2, xg_h2p, xg_h3p, xg_hehp
	! The ghost composition fixed point (ghost_composition_fixed_point_tol):
	! whether the contract covers this run, whether it has been reached, how
	! many applications of the composition map have run, the state this sweep
	! was entered at, the ghost's particle plus electron count entering and
	! leaving an application, and the two moves they give.
	logical :: ghost_contract_on, ghost_fp_reached, bg_ready_entry
	integer :: ghost_fp_pass
	real*8  :: ghost_fp_move, ghost_thermal_move
	real*8  :: f_sp_sweep_entry(1-Ng:N+Ng,n_species)
	real*8  :: np_ghost_entry(1-Ng:0), np_ghost_return(1-Ng:0)

	real*8, dimension(1-Ng:N+Ng),   intent(in) :: T_in
	! CHEMISTRY PRESERVES THE DENSITY IT IS GIVEN
	! (docs/b1_target_system_20260906.md T2.1, replacing D0 C2).  Chemical
	! reactions rearrange nucleons and electrons among species and create no
	! mass, so the mass density is an INPUT of this step and not one of its
	! results: n_io is intent(in), the species fractions returned are
	! normalized by it, and the mass sum of the composition is a CHECK
	! against it (rho_recon below), never a correction applied to it.  The
	! species fractions are read on entry and overwritten with the
	! equilibrium result on exit (in place); callers must not alias them.
	real*8, dimension(1-Ng:N+Ng),   intent(in) :: n_io
	real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp_io
	! What THIS sweep found, for a caller that needs to know whether the state
	! it handed in can be described at all: the same fields the run-wide
	! ledgers accumulate, for this one sweep. The steady solver reads
	! n_nonfinite and n_offsimplex off it to decide whether a probe or a
	! line-search trial is usable (docs/Update_EXHALE_stage1.md section 121); the
	! marching loop does not ask.
	type(ioniz_eq_ledger), optional, intent(out) :: sweep_ledger
	! DIAGNOSTIC ONLY (T2.2, T2.3): the mass density reconstructed from the
	! species the sweep returns, in the same adimensional units as n_io.  The
	! declared mass convention is that the electron mass is carried with its
	! ion (decision 2), so an ionization moves no mass between species and
	! this sum reproduces n_io whenever the element totals are conserved.
	! The caller compares the two; nothing here uses the reconstruction.
	real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: rho_recon

	real*8, dimension(1-Ng:N+Ng) ::  T_K      ! Dimensional temperature
	! Share of an H2 vibrational excitation collisionally de-excited to heat
	real*8, dimension(1-Ng:N+Ng) ::  f_vib_quench, e_vib_bound
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhi,nhii,                   & ! Species densities
	                                 nhe,nhei,nheii,nheiii,nheiTR,  &
	                                 ne,n_in_dim,n_tot
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng) ::  xion
   ! Metal ion densities in canonical species_table order (col im maps
   ! to f_sp column mion_fsp(im)); used throughout in place of named
   ! scalars for each ion so the driver scales with the number of metals.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  nm
   ! Total density of each metal element (canonical element order),
   ! held constant across the cell sweep (sum of its three stages).
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_tot
   ! ELEMENT total of each metal element, i.e. every nucleus of it in the
   ! cell whatever carries it. It differs from nm_tot only for oxygen and
   ! carbon, and only with the oxygen chemistry on: nm_tot is what the metal
   ! ionization block is solved against (the FREE oxygen family, and the
   ! carbon not locked in CO), nm_el is the conserved element. Identical
   ! arrays without the option.
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_el
   ! Carbon monoxide density of each cell [cm^-3], from the CO <-> C + O
   ! chemical equilibrium, and the free oxygen family it leaves behind.
   real*8, dimension(1-Ng:N+Ng) ::  nCO_cell
   ! Row of the first oxygen-carrier unknown.
   integer :: iox
   ! Mass density of the returned composition, reconstructed from the species
   ! mass table.  A check against n_io, never a replacement for it (T2.2).
   real*8, dimension(1-Ng:N+Ng) ::  n_rec

   ! Metal recombination and collisional ionization rates for each ion from
   ! eval_cool (canonical order); bridged to the named rc*/a_ion_*
   ! scalars below for the params packing.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m,aion_m

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
   real*8, dimension(1-Ng:N+Ng) ::  P_H2      ! (molecular; zero unless mol)
   ! Dissociative part of P_H2 (H2 + hv -> H + H+ + e-): a SUBSET of it.
   real*8, dimension(1-Ng:N+Ng) ::  P_H2_di
   ! Double-ionization part (H2 + hv -> H+ + H+ + 2e-) and neutral-
   ! dissociation part (H2 + hv -> H + H) of the same P_H2: two further
   ! SUBSETS of it, disjoint from P_H2_di and from each other. Both are
   ! identically zero unless their options are on.
   real*8, dimension(1-Ng:N+Ng) ::  P_H2_dd, P_H2_nd
   ! Metal photoionization rates for each ion (canonical order) from PH_heat.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m
                       	
   ! Heating and cooling of the composition this sweep RETURNS: both are
   ! assembled after the cell sweep, from the post-sweep densities and the
   ! rates below (see the assembly block after the sweep).
   real*8, dimension(1-Ng:N+Ng) ::  heat,cool
   ! Photoheating rate of ONE particle of each absorber [erg s^-1] in the
   ! attenuated field of each cell: a property of the radiation field, built
   ! from the entry columns, and the piece of the photoheating that does NOT
   ! move when the sweep changes the composition.
   real*8, dimension(1-Ng:N+Ng) ::  h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  h1_m
   ! THE XUV FIELD, SOLVED WITH THE COMPOSITION IN ONE TRAVERSAL.
   ! N*_col carry the column of each absorber inward as the sweep goes, and
   ! hold the columns of the composition the sweep HAS ALREADY SOLVED;
   ! N*_face(j) is the column outside cell j, which is the depth its
   ! star-ward face sees.  nheiS_face is the He I ground singlet of the
   ! entry composition, formed once per cell for the column and the cell's
   ! own optical depth.  See the block loop below for why this is a sweep.
   real*8 ::  N1_col,N15_col,N2_col,NTR_col,NH2_col_xuv
   real*8, dimension(n_mphot) ::  Nm_col_xuv
   real*8, dimension(1-Ng:N+Ng) ::  N1_face,N15_face,N2_face,NTR_face,     &
                                    NH2col_face, nheiS_face, nh2_entry
   real*8, dimension(1-Ng:N+Ng,n_mphot) ::  Nm_face
   real*8 ::  N1_blk,N15_blk,N2_blk,NTR_blk,NH2_blk
   real*8, dimension(n_mphot) ::  Nm_blk
   integer :: jb_lo, jb_hi
   ! What the field routine returns for one cell that is not already a
   ! whole-grid array here.
   real*8 ::  Pm_row(n_mion), h1m_row(n_mion), chan_row(6)
   real*8 ::  heat_row, q_abs_row
   ! Constants of the molecule and the switches of the photoelectron
   ! partition, resolved once per call rather than per cell.
   real*8  ::  D0_H2_xuv, E_ker_H2_dd_xuv
   logical ::  sec_on_xuv, mol_sec_xuv, has_h2_xuv
   ! Coefficients returned by the post-sweep eval_cool call, which is made
   ! for the cooling SUM alone: the sweep was solved with the pre-sweep
   ! coefficients above, and these are discarded.
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB_post,rcheiiB_post,rcheiiiB_post
   real*8, dimension(1-Ng:N+Ng) ::  aHI_post,aHeI_post,aHeII_post,aTR_post
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m_post,aion_m_post

   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR

   ! He recombination radiation -> H ionization coupling scratch
   ! (use_he_rec_coupling; zero-effect when off).
   real*8, dimension(1-Ng:N+Ng) ::  rcheiiB_hrc,dP_HI_hrc,dheat_hrc
   real*8, dimension(1-Ng:N+Ng) ::  dP_H2_hrc
   ! Metal share of the He recombination photons (item P34), added to P_m
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  dP_m_hrc

	real*8, dimension(1-Ng:N+Ng) :: q13,q31a,q31b,Q31
	real*8 :: A31

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII,a_ion_HeITR

	! Equilibrium system setup
   real*8 :: tol,dpmpar
   real*8, dimension(60) :: params
   ! Each element's metal coefficients for ion_system_HeH_metals (canonical order),
   ! built per cell from the 2D rate arrays and handed to set_metal_coeffs.
   real*8, dimension(n_melem) :: meg_ntot,meg_g0,meg_g1,meg_b0,meg_b1, &
                                 meg_a1,meg_a2
   ! Highest stage per element handed to set_metal_coeffs (2 = three-stage,
   ! 1 = two-stage); see species_table::melem_top.
   integer, dimension(n_melem) :: meg_top
   integer :: i0,top
   ! Base index of the first metal element's X+ unknown in sys_x: 4 normally,
   ! 5 when the He triplet occupies x(4) (merged HeITR+metals system), and
   ! metal_row_base() = 8 or 9 in the molecular layout, where x(4..7) are the
   ! H nuclei bound into H2/H2+/H3+/HeH+ and x(8) the triplet.
   integer :: mbase
   ! Molecular solve (thereis_mol): count of cells whose roots all left the
   ! physical simplex, so the closest one was clamped back onto the element
   ! budget, and the local gas pressure [bar] the H2 dissociation
   ! equilibrium of the retry seed is evaluated at.
   integer :: n_mol_clamped
   ! Histogram of the molecular hybrd1 exit code over the sweep (see
   ! ieq_n_mol_info above); combined over threads by the reduction below.
   integer :: n_mol_info(0:5)
   real*8  :: pbar_loc

   ! Ionization solve: validation of the returned root.
   ! The equilibrium systems are polynomial and possess roots outside the
   ! physical simplex (negative stage fractions, or ionized stages of one
   ! element summing above its nucleus total); such a root makes the neutral
   ! density negative and is not a solution of the physical problem even when
   ! the algebraic residual vanishes. Each cell therefore tries up to three
   ! starting points and keeps the best admissible root; x_root_best holds it,
   ! ok_rank grades the attempt (0 none, 1 physical, 2 physical + converged),
   ! and the counters feed the one-line sweep summary.
   ! Largest N_eq over the layouts: molecular (7) + triplet (1) + the two
   ! oxygen carriers + two per metal element. The oxygen pair was added on
   ! 2026-08-30 with the oxygen chemistry, and leaving it out was a real
   ! out-of-bounds write, not a spare-capacity question: with every option
   ! on, N_eq = 30 against a bound of 28, and the copy of the root into
   ! x_root_best ran two elements past the end of a thread-private stack
   ! array. -O3 did not notice; run_fcheck.sh did.
   integer, parameter :: n_x_max = 10 + 2*n_melem
   integer :: iatt, info_ieq, ok_rank, best_rank
   integer :: n_ieq_reseed, n_ieq_retry, n_ieq_unphys, n_ieq_fail
   logical :: conv_ieq, phys_ieq
   real*8, dimension(n_x_max) :: x_root_best
   ! Molecular cell whose roots all left the physical simplex: how far the
   ! least-offending one lies outside it, so the closest to a state can be
   ! selected and clamped back onto the element budget.
   real*8  :: viol, viol_best
   ! Acceptance statistics of this sweep (section 113): counts and largest
   ! normalized reaction residual by acceptance class (1 = solver-converged
   ! root, 2 = root without solver convergence, 3 = projected/handback
   ! state that rechecks as a root, 4 = non-root under the relaxation
   ! amnesty, 5 = root of the constrained element-conserving continuation
   ! solve, 6 = non-root ABOVE the amnesty cap, for which the cell keeps
   ! the composition it entered the sweep with), and the residual-decade
   ! histograms of the physical iterates.
   integer :: n_acc(6), hist_conv(0:15), hist_uncv(0:15)
   real*8  :: acc_resmax(6)
   integer :: acc_class
   real*8  :: acc_res, res_att
   ! Root selection (section 113): residual of the kept rank-1/2 root, the
   ! first bounded iterate whose residual failed the tolerance (the state
   ! of last resort, matching the old rank-1 order), and whether any
   ! attempt stayed inside the simplex / supplied a clamp candidate.
   real*8  :: res_best, res_nonroot
   real*8, dimension(n_x_max) :: x_nonroot
   integer :: info_best, info_nonroot
   logical :: have_nonroot, have_clamp
   ! The composition the cell ENTERED the sweep with: the previous adopted
   ! step's, or the initial condition at step 0. It is the state a cell
   ! keeps when the only candidate the solve produced is a non-root whose
   ! residual is above ieq_nonroot_res_cap (acceptance class 6), and it is
   ! also the starting point of the solve for every step after the first,
   ! so it is built once per cell and used for both.
   real*8, dimension(n_x_max) :: x_entry

   ! THE CELL'S OWN FIELD, SOLVED WITH THE CELL'S OWN COMPOSITION
   ! (xuv_self_field_passes).  it_self counts the passes, x_self holds the
   ! stage fractions the previous one accepted, self_moved is how far the
   ! accepted state moved between two of them, and last_self marks the pass
   ! whose acceptance is the cell's, i.e. the one the ledger records.
   ! nheiS_it and nh2_it are the ground singlet and the molecular hydrogen
   ! of the iterate, the two absorber densities that are not stored as
   ! stage densities of their own.
   integer :: it_self
   logical :: last_self
   ! Whether the state this pass accepted advances the non-root streak
   ! (sweep_advances_the_streak).
   logical :: state_is_recorded
   real*8  :: self_moved, nheiS_it, nh2_it
   ! Passes this cell may take: the self-field passes for every cell, and
   ! the closure passes of the base handoff for a lower ghost that carries
   ! one (see the ghost reservoir block of the sweep).
   integer :: n_self_max
   ! Ionized hydrogen-nucleus fraction of a lower-boundary ghost, against
   ! which the base handoff's molecular partition is stated, the partition
   ! imposed from it, and the residual of the pair at the state the cell
   ! returned (see the ghost reservoir block of the sweep).
   real*8  :: x_ion_ghost, x_h2_imposed, ghost_closure_move
   ! The fixed-point map of that closure at the previous pass, for the
   ! secant step on it: the partition imposed then and the residual
   ! x2 (1 - x_ion) - x_H2 it returned (see the ghost reservoir block).
   real*8  :: x_h2_secant_prev, f_h2_secant_prev, x_h2_map, f_h2_map
   real*8  :: x_h2_secant
   logical :: have_h2_secant_prev
   real*8  :: ghost_closure_res
   logical :: ghost_closed, ghost_cell_closure
   real*8, dimension(n_x_max) :: x_self

   ! Promotion of a molecular cell with no root to the constrained
   ! element-conserving continuation solve (section 11.5-B): the candidate
   ! it returns, whether the continuation reached the cell's full radiation
   ! field, whether the same judge marked the candidate a root, its
   ! residual, and the cost -- hybrd1 calls and wall time. The clock counts
   ! are integer(8) because a nanosecond-resolution system_clock overflows
   ! the default integer in a few seconds.
   real*8, dimension(n_x_max) :: x_cce
   logical :: cce_full, cce_root
   integer :: n_cce_fs
   real*8  :: res_cce
   integer(8) :: clk_beg, clk_end, clk_rate
   ! Sweep totals of the same, combined over threads by the reduction below.
   integer :: n_cce_attempt, n_cce_root, n_cce_solve
   real*8  :: cce_seconds
   ! Admissibility of the state this sweep was handed (section 121):
   ! cells whose accepted state carries a normalized reaction residual that
   ! is not a finite number, cells at which no starting point stayed inside
   ! the element simplex, and the largest element-budget violation seen
   ! before a clamp. Combined over threads by the same reduction.
   integer :: n_res_nonfinite
   real*8  :: viol_sweep_worst
   ! This sweep's totals, handed to the caller and added to the ledger of
   ! whatever state kind the sweep was tagged with.
   type(ioniz_eq_ledger) :: ledger_sweep
   ! A1 elemental census of the state this sweep was handed (element_census).
   type(element_census_state) :: cen_ieq

   ! Heating and cooling of the state this routine RETURNS (n_io, f_sp_io),
   ! adimensional: assembled after the cell sweep from the post-sweep
   ! densities and the pre-sweep rates.
   real*8, dimension(1-Ng:N+Ng),intent(out) :: heat_out,cool_out
   ! Heating efficiency, deposited over absorbed photon energy. Both its
   ! numerator and its denominator are contractions with the composition
   ! the field of the cell was last built at, which is the entry
   ! composition at xuv_self_field_passes = 1 and the last self-field
   ! iterate above it. It is a diagnostic only (nothing in the solve or the
   ! energy update reads it), so it is left at that state rather than
   ! re-formed from the composition the sweep returns.
   real*8, dimension(1-Ng:N+Ng),intent(out) :: q

   !----------------------------------------------------------!      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))
   ! Read once for the run, and announced only when it is armed, so that a
   ! run without the hook is the run without the code.
   if (.not. ieq_inner_tol_read) then
      ieq_inner_tol      = composition_solve_tolerance(tol)
      ieq_inner_tol_read = .true.
      if (ieq_inner_tol .ne. tol)                                         &
         write(*,'(A,ES11.3,A,ES11.3)') ' (ioniz_eq) EXHALE_IEQ_TOL:'//   &
              ' the inner composition solve stops at ', ieq_inner_tol,    &
              ' in place of ', tol
   endif
   tol = ieq_inner_tol

	! ---- THE COMPOSITION THE LOWER GHOST CELLS ENTER THIS SWEEP WITH ----
	!
	! Default: the rows the caller installed, which for a restart are the
	! composition of the gas the reservoir holds at the base level, with the
	! partition within each element taken from the first physical cell
	! (docs/lhs1140b_stationary_D5b2_20260918.md section 7). With
	! EXHALE_GHOST_COMPOSITION_SEED set they are replaced by the stated
	! composition, so that the same ghost system is solved from a different
	! starting point and the seed dependence of the ghost it returns can be
	! measured. VALIDITY: a measurement of the ghost solve; a run with the
	! key set is entered at a composition the boundary model does not state.
	! With the key unset lower_ghost_seed_rows returns the rows it was given
	! and the branch is not entered, so no number moves.
	!
	! It stands before the element census because it changes the ghost's
	! elemental ratios, which the census asserts are invariant ACROSS the
	! sweep and not across a seeding.
	if (ghost_composition_seed_armed()) then
		call lower_ghost_seed_rows(f_sp_io(1-Ng:0,:), f_ghost_seed,       &
		                           ghost_seed_used)
		f_sp_io(1-Ng:0,:) = f_ghost_seed
		if (.not. ghost_seed_announced) then
			write(*,'(A)') ' (ioniz_eq) the lower ghost cells enter the'//&
			     ' composition sweep at '//trim(ghost_seed_used)
			ghost_seed_announced = .true.
		endif
	endif
	call write_ghost_seed_rows(f_sp_io(1-Ng:0,:))

	! The composition the cells above the two lower ghosts were handed, kept
	! where the diagnostic policy above asks for it to be returned.
	if (interior_composition_is_held()) then
		if (.not. allocated(f_sp_interior_entry))                           &
			allocate(f_sp_interior_entry(1:N+Ng,n_species))
		f_sp_interior_entry = f_sp_io(1:N+Ng,:)
	endif

	! A1 element-budget assertion around the whole sweep.  The sweep moves
	! nuclei between stages and carriers and creates none, so every ratio
	! n_El/n_H is invariant across it.  Since T2.1 the density is an input
	! the sweep may not write, so the ABSOLUTE nucleus densities are
	! invariant too and the verify below asserts them (rho_is_fixed).  The label
	! carries the sweep tag, so a failure names which caller (marching, the
	! steady iterate, or a steady candidate) broke the budget.
	! Off unless EXHALE_ELEMENT_ASSERT is set.
	if (ieq_sweep_state_kind .eq. ieq_state_marching) then
		call element_census_take('ioniz_eq [marching]', n_io, f_sp_io,    &
		                         cen_ieq)
	else if (ieq_sweep_state_kind .eq. ieq_state_steady_iterate) then
		call element_census_take('ioniz_eq [steady iterate]', n_io,       &
		                         f_sp_io, cen_ieq)
	else
		call element_census_take('ioniz_eq [steady candidate]', n_io,     &
		                         f_sp_io, cen_ieq)
	endif

	! ---- THE GHOST COMPOSITION IS A FIXED POINT OF THIS SWEEP ----
	!
	! The contract (ghost_composition_fixed_point_tol). A sweep is a map from
	! the composition it is entered at to the composition it returns: the
	! entry composition sets the columns, the cells' own depths, the electron
	! and particle counts, the secondary-ionization partition, the shielding
	! and the ionization the base handoff's partition is stated against, and
	! the cells are then solved in that context. The lower boundary is built
	! on the two lower ghost cells, so a ghost that is not a fixed point of
	! this map is a boundary that is a function of the composition the run was
	! entered at, and MEASURED it is: two entry compositions return ghosts
	! 4.7e-3 to 5.5e-2 apart in their trace ions, both accepted by every test
	! the sweep applies (READ, docs/lhs1140b_p1_step1b_20260919.md section 1).
	!
	! The map is therefore applied again, at the ghost it returned and at the
	! composition the cells above the ghosts were HANDED, until the ghost and
	! its particle plus electron count stop moving within the accuracy. The
	! interior is restored at each application, so the composition this
	! routine returns for the physical cells is ONE application of the map,
	! as it has always been, and the steady residual built on it is one
	! application of one operator.
	!
	! A run whose ghost is already at its fixed point -- the usual case, since
	! the ghost a sweep is entered at is the one the previous sweep returned
	! -- takes ONE application and costs nothing.
	!
	! THE APPLICATIONS ARE ONE SWEEP, so nothing that separates one sweep from
	! the next may advance between them. The frozen background's readiness is
	! such a thing: it is false until a sweep has filled it, and it is what
	! turns the transported carriers from something this sweep INITIALIZES
	! into something handed to it (the imposed-fraction block of the cell
	! loop). Letting it turn over inside one sweep would make the second
	! application a later sweep: MEASURED on mol_carrier at step 0, 185 cells
	! then leave the physical simplex and every cell of the wind rests on a
	! non-root. It is held at the value the sweep was entered with and takes
	! the sweep's own value when the sweep ends.
	ghost_contract_on = thereis_mol .and. thereis_He
	bg_ready_entry    = bg_ready
	ghost_fp_pass     = 0
	ghost_fp_move     = 0.0d0
	ghost_thermal_move = 0.0d0
	ghost_fp_reached  = .not. ghost_contract_on
	if (ghost_contract_on) f_sp_sweep_entry = f_sp_io

	ghost_fixed_point: do
	ghost_fp_pass            = ghost_fp_pass + 1
	sweep_advances_the_streak = (ghost_fp_pass .eq. 1)
	bg_ready                 = bg_ready_entry

	!----------------------------------!
	
	! Preliminary profiles exctraction
	
	! Dimensional total number density profile and temperature
	n_in_dim = n_io*n0
	T_K      = T_in*T0
		
	! Extract species profiles
	nhi    = f_sp_io(:,1)*n_in_dim    ! HI
	nhii   = f_sp_io(:,2)*n_in_dim    ! HII
	if (thereis_He) then
		nhei   = f_sp_io(:,3)*n_in_dim    ! HeI
		nheii  = f_sp_io(:,4)*n_in_dim    ! HeII
		nheiii = f_sp_io(:,5)*n_in_dim    ! HeIII
		nheiTR = f_sp_io(:,6)*n_in_dim    ! HeITR
	else
		! Enforce condition of zero helium
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0
	endif

    ! Metal ion densities in canonical species_table order (col im maps to
    ! f_sp column mion_fsp(im)).
    do im = 1,n_mion
       nm(:,im) = f_sp_io(:,mion_fsp(im))*n_in_dim
    enddo

	! molecular species (zero when thereis_mol is off)
	if (thereis_mol) then
		nmol_eq(:,1) = f_sp_io(:,isp_H2)  *n_in_dim
		nmol_eq(:,2) = f_sp_io(:,isp_H2p) *n_in_dim
		nmol_eq(:,3) = f_sp_io(:,isp_H3p) *n_in_dim
		nmol_eq(:,4) = f_sp_io(:,isp_HeHp)*n_in_dim
	endif

	! oxygen carriers (zero when thereis_oxychem is off)
	if (thereis_oxychem) then
		nox_eq(:,1) = f_sp_io(:,isp_OH) *n_in_dim
		nox_eq(:,2) = f_sp_io(:,isp_H2O)*n_in_dim
		nox_eq(:,3) = f_sp_io(:,isp_CO) *n_in_dim
	endif

    ! Total density of each metal element (sum of its three stages, in the
    ! same neutral+singly+doubly order as the original nC=nci+ncii+nciii).
    do im = 1,n_melem
       i0 = melem_i0(im)
       if (melem_top(im) .ge. 2) then
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1) + nm(:,i0+2)
       else
          ! Two-stage element: neutral + singly ionized only.
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1)
       endif
    enddo
    nm_el    = nm_tot
    nCO_cell = 0.0d0

    ! ---- oxygen and carbon element totals, and the CO reservoir ----
    ! With the oxygen chemistry on, the O I / C I columns hold FREE ATOMIC
    ! oxygen and carbon: the rest of each element sits in OH, H2O and CO.
    ! The element total -- the quantity that is conserved and that the
    ! abundance melem_ab fixes -- is therefore the stage sum PLUS the
    ! molecular carriers, counted through bsp_nO / bsp_nC of the species
    ! table. Without this the closure
    !   n_O,tot = n(OI)+n(OII)+n(OIII)+n(OH)+n(H2O)+n(CO)
    ! is not closed and half the element leaks away one sweep at a time.
    !
    ! CO is then taken out of both elements at the CO <-> C + O chemical
    ! equilibrium of the local (n, T) and frozen there for the solve
    ! (decision D4: CO reacts with nothing in the audited set, but it holds
    ! 45-46% of the oxygen at every measured level, so a network that gave
    ! the whole oxygen abundance to the water family would over-supply the
    ! OH cycle by about a factor two). What is left is the FREE OXYGEN
    ! FAMILY shared by O I/II/III, OH and H2O, and the carbon left for the
    ! metal block's C I/II/III.
    if (thereis_oxychem) then
       nm_el(:,iel_O) = nm_tot(:,iel_O)                                  &
                      + nox_eq(:,1) + nox_eq(:,2) + nox_eq(:,3)
       nm_el(:,iel_C) = nm_tot(:,iel_C) + nox_eq(:,3)
       if (carrier_transport) then
          ! CO is a TRANSPORTED reservoir: the carrier operator moved it,
          ! and this sweep takes the CO the cell has rather than the CO the
          ! local (n, T) would make. The difference is not cosmetic -- the
          ! equilibrium form re-forms CO in a cool outer wind, where a
          ! transported one would have been carried out of the molecular
          ! layer and quenched.
          nCO_cell = max(nox_eq(:,3), 0.0d0)
          nCO_cell = min(nCO_cell, nm_el(:,iel_C))
          nCO_cell = min(nCO_cell, nm_el(:,iel_O))
       else
          do j = 1-Ng,N+Ng
             nCO_cell(j) = co_equilibrium_density(nm_el(j,iel_C),        &
                                                  nm_el(j,iel_O), T_K(j))
          enddo
       endif
       nm_tot(:,iel_O) = max(nm_el(:,iel_O) - nCO_cell, 0.0d0)
       nm_tot(:,iel_C) = max(nm_el(:,iel_C) - nCO_cell, 0.0d0)
       nox_eq(:,3)     = nCO_cell
    endif
	
	! Total number densities: H and He NUCLEI totals, every molecular and
	! oxygen carrier counted with its nucleus multiplicity (the chemistry
	! conserves elements).  ONE definition, shared with the heating dump and
	! the advection-corrected post-process -- see utils.
	call hydrogen_helium_nuclei_density(nhi,nhii,nhei,nheii,nheiii,       &
	                                    nh,nhe,nmol_eq,nox_eq)

	! Free electron density (assuming overall neutrality; nm adds the
	! metal electrons under the eos_metals policy, nmol_eq the molecular-ion
	! electrons -- it is zero for an atomic run, so the sum is unchanged there)
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol_eq)

	! Ionized fraction handed to the photoelectron partition. Dalgarno, Yan &
	! Liu (1999) section 7 define it as "the number density ratio of the
	! electrons to the hydrogen and helium nuclei", so the numerator is the
	! TOTAL free electron density -- molecular ions and metal electrons
	! included, and He III counted twice -- not a count of H and He ions.
	! Metals are an extension beyond both sources' composition and are
	! marked as such at the code site (electron_energy_degradation.f90);
	! they can push the ratio above one, which the clamp absorbs.
	xion = min(max(ne/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Total gas-particle density (electrons excluded), i.e. the density of
	! third bodies M for the three-body molecular reactions R12/R13/R15, and
	! -- with the electrons added back -- the gas pressure p = (n_tot+n_e)kT
	! of the chemical-equilibrium retry seed below.  It is NOT n_in_dim:
	! n_in_dim = rho/m_H is a MASS density in m_H units (calc_rho weights
	! each species by bsp_mass), so it over-counts the particles by the mean
	! particle mass -- 2.3x at an H2-rich base, up to 4x in a He-dominated
	! one.  Only the molecular path reads it, so an atomic run is unchanged.
	! The GAS PARTICLE DENSITY of every cell, from the one stoichiometric
	! sum of utilities.  It is a property of the gas and not of the
	! molecular option: an atomic mixture has a particle density too, the
	! pressure the retries below form is (n_tot + n_e) kB T in any gas, and
	! the frozen background the transport-chemistry operator reads carries
	! it cell by cell.  The molecular and oxygen densities enter the sum
	! where they exist and are zero where they do not.
	if (thereis_oxychem) then
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq,      &
		               nox_eq)
	else
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq)
	endif
	! The ghost's heavy-particle plus electron count at the composition this
	! application was ENTERED at, against which the count it returns is
	! measured (the thermodynamic leg of the ghost composition contract).
	if (ghost_contract_on) np_ghost_entry = n_tot(1-Ng:0) + ne(1-Ng:0)

	! Cell-by-cell pressure-broadening factor for the opacity ('P' model).
	! opacity_pT_factor returns 1.0 for all other models, so opa_pf=1
	! and the column densities are unchanged (bit-identical).
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

	!----------------------------------!

	! The FUV photon field of the molecular layer, for ALL of its absorbers
	! at once: H2 in the Lyman and Werner lines (Draine & Bertoldi 1996) and
	! H2O and OH in their continua (the A2 oxygen option).  Over 912-1110 A
	! the three share one beam, so they are solved together --
	! fuv_lw_photon_field in utils_ion_eq states why and how the shared beam
	! splits.  The columns use the same radial integration and opa_pf
	! weighting as every other absorber column, and the incoming (pre-solve)
	! densities, like the photoionization columns built inside PH_heat_HHe.
	call fuv_lw_photon_field(nmol_eq(:,1), nox_eq(:,2), nox_eq(:,1),      &
	                         nox_eq(:,3), nhi,                            &
	                         T_K, nh,                                     &
	                         NH2_col_lw, NH2O_col, NOH_col, NCO_col,      &
	                         f_shield_lw, tr_lines_lw, tau_fuv,           &
	                         k_lw_diss, p_lw_single, p_lw_absorbed,       &
	                         k_co_diss, theta_co_shield,                  &
	                         j_h2o_fuv, j_oh_fuv, a_lines_lw,             &
	                         lw_col_over_overlap)
	if (thereis_mol .and. F_LW_star .gt. 0.0d0 .and.                      &
	    .not. warned_lw_fit_range) then
		if (lw_col_over_overlap .gt. 1.0d0 .or.                            &
		    maxval(NH2_col_lw) .gt. NH2_db96_max) then
			warned_lw_fit_range = .true.
			write(*,'(a,es9.2,a,f6.2,a,es9.2,a)') ' (ioniz_eq)'//         &
			  ' WARNING: the star-ward H2 column reaches ',               &
			  maxval(NH2_col_lw), ' cm^-2, which is ',                    &
			  lw_col_over_overlap, ' times the top of the column axis'// &
			  ' of the overlapping-line self-shielding table, above'//    &
			  ' which its edge value is returned rather than a'//         &
			  ' calculated one, and the band share the H2 lines remove'//&
			  ' is on a fit demonstrated only to ', NH2_db96_max,         &
			  ' cm^-2 (Draine & Bertoldi 1996).'
		endif
	endif

	!----------------------------------!

    !---- The inputs of the photoelectron partition ----!

	! f_vib_quench is the share of an H2 vibrational excitation that is
	! collisionally de-excited into heat rather than radiated away, and
	! e_vib_bound the mean internal energy of the level a B or C
	! fluorescence lands on. The photoelectron partition needs both and does
	! not carry the temperature, so they are formed here, from T and the
	! pre-solve densities, for the whole grid; the field routine reads them
	! cell by cell. Both are zero in a run without molecular hydrogen, which
	! is what switches the two vibrational heat channels off.
	f_vib_quench = 0.0d0
	e_vib_bound  = 0.0d0
	if (thereis_He .and. thereis_mol) then
		! The three colliders of the cascade are atomic hydrogen, helium and
		! H2, each with its own published coefficient
		! (h2_vibrational_relaxation).  The helium one is the GROUND-SINGLET
		! neutral: the quantum calculation behind that coefficient is for
		! neutral He + H2, and a helium ion interacts with the molecule
		! through a different potential and carries no coefficient here.
		f_vib_quench = h2_vibrational_heat_fraction(T_K, nhi,             &
		                   nmol_eq(:,1),                                  &
		                   he_ground_singlet_density(nhei, nheiTR))
		e_vib_bound  = h2_energy_per_bound_fluorescence_eV(T_K)
	endif

	! Switches and molecular constants of the field, resolved once.
	sec_on_xuv      = use_sec_ion .and. sec_ion_active
	has_h2_xuv      = thereis_He .and. thereis_mol
	mol_sec_xuv     = sec_on_xuv .and. has_h2_xuv
	D0_H2_xuv       = h2_dissociation_energy_eV()
	E_ker_H2_dd_xuv = h2_double_fragment_kinetic_energy()

	nh2_entry = 0.0d0
	if (has_h2_xuv) nh2_entry = nmol_eq(:,1)

	! The one-particle photoheating rates of the absorbers this run does not
	! carry stay zero, and so does the metal block until the sweep fills it.
	P_m     = 0.0d0
	h1_HeI  = 0.0d0
	h1_HeII = 0.0d0
	h1_HeTR = 0.0d0
	h1_H2   = 0.0d0
	h1_m    = 0.0d0


   ! The recombination and collisional-ionization rate coefficients the
   ! cell sweep below solves with. They depend on T alone, which is what
   ! lets them be built before the sweep and used unchanged by it.
   !
   ! ONLY the coefficients. The cooling is a contraction of them with the
   ! COMPOSITION, and the sweep is about to replace the composition, so the
   ! cooling this routine returns is assembled after the sweep and from the
   ! state the sweep produced (the eval_cool call below the sweep). The
   ! cooling of the entry composition is not needed by anything and is not
   ! formed: nothing but that second evaluation contributes to the returned
   ! `cool`, because the marching temperature update rebuilds the cooling
   ! from eval_cool alone, so any term added on top would be a source the
   ! marching relaxes without and the steady residual demands.
    call chemical_rate_coefficients(T_K,                                  &
	                    rchiiB,rcheiiB,rcheiiiB, rec_m,                   &
	                    a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,            &
	                    a_ion_HeITR=a_ion_HeITR)

	! Capture the ground-state H proton-budget coefficients on
	! every pass (the converged pass is the one read out by write_excited_H).
	if (use_excited_H) then
		cion_HI  = a_ion_HI          ! collisional ionization [cm^3 s^-1]
		arec_HII = rchiiB            ! case-B recombination  [cm^3 s^-1]
	endif

	! Charge-exchange rate coefficients are evaluated per cell below
	! (cx_set_cell) before each metal ionization solve.

	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten
	endif

	! He recombination radiation ionizing H I and H2 (Draine 2011 emission,
	! absorbed locally: the cell keeps the fraction 1 - exp(-tau_c) of each
	! channel's photons at its own optical depth and the rest leave).
	! Uses the pre-solve (lagged) densities, like the other lagged rate terms;
	! corrects the He II recombination coefficient rcheiiB, adds an H I
	! photoionization rate to P_HI and an H2 one to P_H2. Its photoelectron
	! HEATING is not taken from this call: it carries the absorber densities,
	! so it is re-evaluated after the sweep. Which species takes the photon is
	! the ratio of the absorption coefficients at the energy of the channel,
	! so in a molecular layer, where H2 outnumbers H I by three decades and
	! absorbs 1.2-3.7 times as strongly, essentially all of it goes to H2.
	! The He II recombination *cooling* is kT times the coefficient this
	! balance removes He+ with as far as that coefficient is a function of
	! temperature alone (Cool_coeff: alpha_rec_HeII_total, the two capture
	! channels rcheiiB + rcheiTR when the triplet is tracked). What this
	! coupling adds on top of that -- the weight y_net on the ground-term
	! channel and the 0.25 alpha_B excited-singlet capture -- is a function
	! of the absorber densities, and eval_cool is also evaluated at perturbed
	! temperatures by the semi-implicit energy update, whose dC/dT divides
	! the difference of two eval_cool calls by 1e-5 T; the reason is at
	! alpha_rec_HeII_total. At T = 1e4 K and y_net = 0.81 that leaves 10 per
	! cent of this channel uncharged.
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nmol_eq(:,1), nhei, nheii, nheiTR,  &
		                     ne, nm, A31, q31a, q31b,                      &
		                     rcheiiB_hrc, dP_HI_hrc, dP_H2_hrc,            &
		                     dP_m_hrc, dheat_hrc)
		rcheiiB = rcheiiB_hrc
		! The photoionization rates this coupling adds -- to H I, to H2 and
		! to each metal ion -- are added to the stellar ones cell by cell in
		! the sweep below, where the stellar field of that cell has just been
		! evaluated.
		!
		! THE H2 SHARE IS NOT SPLIT INTO THE DISSOCIATIVE CHANNEL,
		! deliberately. These photons sit at 19.8-24.6 eV, where the
		! dissociative branching of frac_H2_dissociative_ionization is
		! 1.9-2.3%; splitting off that much of a channel that is itself a
		! correction would be well inside the +/-4-5% the measured branching
		! carries. The whole of dP_H2_hrc therefore makes H2+. The same holds
		! for the two channels added later: these photons are far below the
		! 51.4 eV double-ionization threshold and outside the 33-41 eV
		! neutral window, so of the four channels only the H2+ one can
		! receive them, and P_H2_dd / P_H2_nd stay untouched.
	endif


   !----------------------------------!

   ! Ionization equilibrium system solution

	! THE FIELD AND THE COMPOSITION, IN ONE TRAVERSAL OR IN TWO.
	!
	! A cell's photoionization rates are set by the columns of the cells
	! OUTSIDE it and by its own, and by nothing inside it
	! (photoionization_field_at_cell_HHe states the dependence). The sweep
	! runs outside in, so the columns it hands each cell CAN be the columns
	! of the composition it has already solved, which makes the field a cell
	! is solved at the field of that composition in one traversal instead of
	! the field of the composition the previous traversal was handed.
	!
	! xuv_field_block_cells says which. Zero, the default, is one block: the
	! whole grid at the composition the sweep was handed. A positive width
	! carries the solved columns from block to block; inside a block the
	! columns still come from the entry composition, so its cells do not
	! depend on each other and are solved together. The measured cost of
	! each width, and why zero, are at that constant's declaration.
	!
	! THE CELL'S OWN DEPTH IS A SEPARATE LAG, and a separate knob. A cell's
	! dtau is built from the densities the cell is being solved for, so one
	! field build and one solve leave the returned composition implying a
	! slightly different field from the one it was solved in, whatever the
	! block width is. xuv_self_field_passes says how many times the pair is
	! solved against each other inside the cell.
	!
	! THE FIXED POINT IS THE SAME ONE for every setting of both. At the
	! fixed point the composition a cell enters with and the one it returns
	! are the same state, so every width and every pass count builds the
	! same columns and the same field, and the coupled source step of
	! EXHALE_main converges to the same pair. What differs is how many of
	! its passes that takes.
	N1_col      = 0.0d0
	N15_col     = 0.0d0
	N2_col      = 0.0d0
	NTR_col     = 0.0d0
	NH2_col_xuv = 0.0d0
	Nm_col_xuv  = 0.0d0
	jb_hi       = N+Ng
	if (.not.thereis_He) then ! If no helium

		xuv_block_H: do while (jb_hi .ge. 1-Ng)

		jb_lo = 1-Ng
		if (xuv_field_block_cells .gt. 0)                                  &
			jb_lo = max(1-Ng, jb_hi - xuv_field_block_cells + 1)

		! The column each cell of this block sees above it: the running
		! column of the composition already solved, carried on across the
		! block with the entry composition of the cells outside it.
		N1_blk  = N1_col
		N15_blk = 0.0d0
		N2_blk  = 0.0d0
		NTR_blk = 0.0d0
		NH2_blk = 0.0d0
		Nm_blk  = 0.0d0
		call advance_starward_columns(jb_hi, jb_lo, nhi, nhei, nheii,      &
		         nheiTR, nh2_entry, nm, .false.,                           &
		         N1_blk,N15_blk,N2_blk,NTR_blk,NH2_blk,Nm_blk,             &
		         N1_face=N1_face)

	! The H2 thermochemistry table is built here, SERIALLY, if no caller
	! built it yet (the main program does at startup; a test driver may
	! not): the sweep below is the first parallel region that reads it,
	! and the lazy build inside keq_H_H_to_H2 was removed on 2026-09-13
	! (review P2) because its readiness read outside the critical region
	! was a data race.
	if (thereis_mol .and. .not. h2_thermochemistry_ready())               &
		call h2_thermochemistry_init
		!$omp parallel do default(shared) schedule(static)                 &
		!$omp   private(j, heat_row) if(marching_step > 0)
		do j = jb_hi, jb_lo, -1
			call photoionization_field_at_cell_H(j, N1_face(j), nhi(j),    &
			         xion(j), sec_on_xuv, P_HI(j), h1_HI(j), heat_row,     &
			         q(j))
		enddo
		!$omp end parallel do

		! The Balmer continuum of H(n=2) adds an effective H I
		! photoionization rate [s^-1] (a proton source). gph_balmer_HI is
		! filled from the previous converged outer pass by
		! excited_hydrogen::excited_H_update, so the coupling is decoupled
		! from the hydro sub-step and carries that documented lag; its
		! heating is added after the sweep with every other heating term.
		! Default-off => no change.
		if (use_excited_H) then
			do j = jb_hi, jb_lo, -1
				gph_ground_HI(j) = P_HI(j)   ! the pure ground-state rate
				P_HI(j) = P_HI(j) + gph_balmer_HI(j)
			enddo
		endif

		! Ionization solves in each cell are independent (the marching_step>0 warm-start uses
		! this cell's own previous-step value), so the sweep is OpenMP-parallel
		! over cells. sys_x/wa/info are threadprivate (global_parameters); only
		! the subroutine-local scratch is private. marching_step==0 runs serial (the
		! clause) because its first-step warm-start reads the neighbor cell.
		!$omp parallel do default(shared) schedule(dynamic,8)                  &
		!$omp   private(params, usednt, it_self, last_self, self_moved,        &
		!$omp           x_self, heat_row) if(marching_step > 0)
		do j = jb_hi, jb_lo, -1

			! Lazily allocate this thread's threadprivate NL scratch.
			if (.not. allocated(sys_x)) allocate(sys_x(N_eq))
			if (.not. allocated(wa))    allocate(wa(lwa))

			self_field_H: do it_self = 1, xuv_self_field_passes

			! The cell's own attenuation at the composition the previous
			! pass returned (xuv_self_field_passes).  The column outside the
			! cell is the block's and does not move here.
			if (it_self .gt. 1) then
				call photoionization_field_at_cell_H(j, N1_face(j), nhi(j), &
				         xion(j), sec_on_xuv, P_HI(j), h1_HI(j), heat_row,  &
				         q(j))
				if (use_excited_H) then
					gph_ground_HI(j) = P_HI(j)
					P_HI(j) = P_HI(j) + gph_balmer_HI(j)
				endif
			endif

			! Ionization equilibrium system setup: named-field cell state
			! (Inc 4). System_H reads these; params is now only the MINPACK
			! transport argument (unread by the converted system).
			ieq_cell%P_HI     = P_HI(j)
			ieq_cell%rchiiB   = rchiiB(j)
			ieq_cell%nh       = nh(j)
			ieq_cell%jcell    = j
			ieq_cell%a_ion_HI = a_ion_HI(j)

			 ! Initial guess
			if (marching_step.le.0) then
				if(r(j).le.(1.5))then
					sys_x(1) = r(j)-0.5
				else
					sys_x(1) = 1.0
				endif
			else

				sys_x(1) = nhii(j)/nh(j)

			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			call solve_ieq(ion_system_H,jac_system_H,N_eq,sys_x,    &
				      params,tol,wa,lwa,usednt)

			! Has the cell's own field stopped moving its own composition.
			self_moved = 1.0d0
			if (it_self .gt. 1)                                            &
				self_moved = maxval(abs(sys_x(1:N_eq) - x_self(1:N_eq)))
			x_self(1:N_eq) = sys_x(1:N_eq)
			last_self = (it_self .ge. xuv_self_field_passes)                &
			            .or. (self_moved .le. xuv_self_field_tol)

			! Extract solution profiles
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)

			! The cell state, kept for the closure evaluator (see
			! ieq_rate_cell). Written here, never read by the sweep.
			ieq_rate_cell(j)      = ieq_cell
			ieq_ne_cell(j)        = ne(j)
			ieq_TK_cell(j)        = T_K(j)
			ieq_ntot_cell(j)      = 0.0d0

			if (last_self) exit self_field_H
			enddo self_field_H

		enddo
		!$omp end parallel do

		! Carry the column past this block with the composition just solved.
		call advance_starward_columns(jb_hi, jb_lo, nhi, nhei, nheii,      &
		         nheiTR, nh2_entry, nm, .false.,                           &
		         N1_col,N15_col,N2_col,NTR_col,NH2_col_xuv,Nm_col_xuv)

		jb_hi = jb_lo - 1
		enddo xuv_block_H

		ieq_neq_stored   = N_eq
		ieq_mbase_stored = 0
		ieq_iox_stored   = 0
		ieq_rates_ready  = .true.

		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0

	else  ! If there is helium

		! Metal unknowns start at sys_x(4) normally, but shift to sys_x(5)
		! when the He triplet occupies sys_x(4) (merged HeITR+metals system).
		! In the molecular layout the molecular unknowns own sys_x(4..7) and
		! the metals follow at metal_row_base() (single definition, shared
		! with the merged molecular+metals residual).
		if (thereis_mol) then
			mbase = metal_row_base()
			iox   = oxygen_row_base()
		else
			mbase = 4
			if (thereis_HeITR) mbase = 5
			iox   = 0
		endif

		! OpenMP-parallel cell sweep (see the no-He branch above). The cell-by-cell
		! metal coefficients (met_*, System_HeH_metals) and charge-exchange rates
		! (cx_kc, cx_metal_base) are threadprivate, so each thread keeps its own;
		! cx_metal_base is broadcast (copyin) and toggled 4<->5 per cell. All the
		! subroutine-local scratch is private. marching_step==0 stays serial (neighbor
		! warm-start).
		! Failed-cell counters (combined via the reduction below); reset once
		! per equilibrium sweep. n_mol_clamped counts molecular cells whose roots
		! were all outside the physical simplex; the n_ieq_* counters are the
		! atomic root validation (see the declarations above).
		n_mol_clamped = 0
		n_mol_info(:) = 0
		n_ieq_reseed = 0
		n_ieq_retry  = 0
		n_ieq_unphys = 0
		n_ieq_fail   = 0
		n_acc(:)      = 0
		acc_resmax(:) = 0.0d0
		hist_conv(:)  = 0
		hist_uncv(:)  = 0
		n_cce_attempt = 0
		n_cce_root    = 0
		n_cce_solve   = 0
		cce_seconds   = 0.0d0
		n_res_nonfinite  = 0
		viol_sweep_worst = 0.0d0
		! What the base handoff's closure reaches in this sweep, reduced
		! over the cells of every block.
		ghost_closure_res_sweep  = 0.0d0
		ghost_closure_pass_sweep = 0

		xuv_block_HHe: do while (jb_hi .ge. 1-Ng)

		jb_lo = 1-Ng
		if (xuv_field_block_cells .gt. 0)                                  &
			jb_lo = max(1-Ng, jb_hi - xuv_field_block_cells + 1)

		! The He I GROUND SINGLET of the entry composition, n(1^1S) =
		! n(He I) - n(2^3S): what the He I opacity, the He I photoelectric
		! heating and the He I secondary-ionization target all act on. Formed
		! here, with the columns, before the cells of the block are solved.
		do j = jb_hi, jb_lo, -1
			nheiS_face(j) = nhei(j)
			if (thereis_HeITR) nheiS_face(j) =                             &
			         he_ground_singlet_density(nhei(j), nheiTR(j))
			if (has_h2_xuv) nh2_entry(j) = nmol_eq(j,1)
		enddo

		! The columns each cell of this block sees above it: the running
		! columns of the composition already solved, carried on across the
		! block with the entry composition of the cells outside it.
		N1_blk  = N1_col
		N15_blk = N15_col
		N2_blk  = N2_col
		NTR_blk = NTR_col
		NH2_blk = NH2_col_xuv
		Nm_blk  = Nm_col_xuv
		call advance_starward_columns(jb_hi, jb_lo, nhi, nheiS_face,       &
		         nheii, nheiTR, nh2_entry, nm, has_h2_xuv,                 &
		         N1_blk,N15_blk,N2_blk,NTR_blk,NH2_blk,Nm_blk,             &
		         N1_face,N15_face,N2_face,NTR_face,NH2col_face,Nm_face)

		!$omp parallel do default(shared) schedule(static)                 &
		!$omp   private(j, Pm_row, h1m_row, chan_row, heat_row, q_abs_row) &
		!$omp   if(marching_step > 0)
		do j = jb_hi, jb_lo, -1
			call photoionization_field_at_cell_HHe(j,                      &
			         N1_face(j),N15_face(j),N2_face(j),NTR_face(j),        &
			         NH2col_face(j),Nm_face(j,:),                              &
			         nhi(j),nheiS_face(j),nheii(j),nheiTR(j),nh2_entry(j),     &
			         nm(j,:), xion(j), f_vib_quench(j), e_vib_bound(j),    &
			         has_h2_xuv, sec_on_xuv, mol_sec_xuv,                  &
			         D0_H2_xuv, E_ker_H2_dd_xuv,                           &
			         P_HI(j),P_HeI(j),P_HeII(j),P_HeITR(j), Pm_row,        &
			         P_H2(j),P_H2_di(j),P_H2_dd(j),P_H2_nd(j),             &
			         h1_HI(j),h1_HeI(j),h1_HeII(j),h1_HeTR(j),h1_H2(j),    &
			         h1m_row, heat_row, chan_row, q(j), q_abs_row)
			P_m(j,:)  = Pm_row
			h1_m(j,:) = h1m_row
		enddo
		!$omp end parallel do

		! The two rates that are added to the stellar field of this cell,
		! both carrying a lag of one sweep and both formed above for the
		! whole grid. Serial, because the second keeps a diagnostic count.
		do j = jb_hi, jb_lo, -1

			! The Balmer continuum of H(n=2): an effective H I
			! photoionization rate [s^-1], i.e. a proton source.
			! gph_balmer_HI is filled from the previous converged outer pass
			! by excited_hydrogen::excited_H_update, so the coupling is
			! decoupled from the hydro sub-step and carries that documented
			! lag; its heating is added after the sweep with every other
			! heating term. Default-off => no change.
			if (use_excited_H) then
				gph_ground_HI(j) = P_HI(j)   ! the pure ground-state rate
				P_HI(j) = P_HI(j) + gph_balmer_HI(j)
			endif

			! He recombination radiation, absorbed locally.
			if (use_he_rec_coupling) then
				! Count where the coupling out-ionizes the stellar field by
				! three decades, before P_HI absorbs it (diagnostic; silent
				! when zero).
				! Counted once per sweep, at the first application of the
				! composition map: the count is of cells in a sweep, not of
				! cells in an application (sweep_advances_the_streak).
				if (j .ge. 1 .and. j .le. N .and.                          &
				    sweep_advances_the_streak) then
					if (dP_HI_hrc(j) .gt.                                  &
					    he_rec_dominant_ratio*P_HI(j)) then
						n_cells_he_rec_photoionization_dominant =          &
							n_cells_he_rec_photoionization_dominant + 1
						he_rec_photoionization_ratio_max =                 &
							max(he_rec_photoionization_ratio_max,          &
							    dP_HI_hrc(j)/max(P_HI(j),1.0d-99))
					endif
				endif
				P_HI(j) = P_HI(j) + dP_HI_hrc(j)
				! The same photons absorbed by H2: an addition to the H2
				! photoionization rate, which drives the H2 destruction row
				! and the H2+ formation row of the molecular system.
				if (thereis_mol) P_H2(j) = P_H2(j) + dP_H2_hrc(j)
				! And the share each metal ion takes of the same photons,
				! an addition to its photoionization rate (item P34).
				if (thereis_metals) P_m(j,:) = P_m(j,:) + dP_m_hrc(j,:)
			endif

		enddo

		!$omp parallel do default(shared) schedule(dynamic,8) copyin(cx_metal_base) &
		!$omp   private(params, usednt, i0, top, im, meg_ntot, meg_g0, meg_g1,      &
		!$omp           meg_b0, meg_b1, meg_a1, meg_a2, meg_top,                     &
		!$omp           pbar_loc, viol, viol_best,                                  &
		!$omp           iatt, info_ieq, ok_rank, best_rank, conv_ieq, phys_ieq,      &
		!$omp           x_root_best, x_entry, acc_class, acc_res, res_att,          &
		!$omp           res_best, res_nonroot, x_nonroot, info_best, info_nonroot,   &
		!$omp           have_nonroot, have_clamp,                                    &
		!$omp           x_cce, cce_full, cce_root, n_cce_fs, res_cce,                &
		!$omp           clk_beg, clk_end, clk_rate,                                  &
		!$omp           it_self, last_self, state_is_recorded, self_moved,           &
		!$omp           x_self, n_self_max,                                         &
		!$omp           nheiS_it, nh2_it, x_ion_ghost, x_h2_imposed,                 &
		!$omp           ghost_closure_move, ghost_closure_res, ghost_closed,         &
		!$omp           ghost_cell_closure, x_h2_secant_prev, f_h2_secant_prev,      &
		!$omp           x_h2_map, f_h2_map, x_h2_secant, have_h2_secant_prev,        &
		!$omp           Pm_row, h1m_row, chan_row, heat_row, q_abs_row)              &
		!$omp   reduction(+:n_mol_clamped,n_mol_info,n_ieq_reseed,n_ieq_retry, &
		!$omp               n_ieq_unphys,n_ieq_fail,n_acc,hist_conv,hist_uncv, &
		!$omp               n_cce_attempt,n_cce_root,n_cce_solve,cce_seconds,  &
		!$omp               n_res_nonfinite)                                   &
		!$omp   reduction(max:acc_resmax,viol_sweep_worst,ghost_closure_res_sweep,  &
		!$omp               ghost_closure_pass_sweep) if(marching_step > 0)
		do j = jb_hi, jb_lo, -1

			! Lazily allocate this thread's threadprivate NL scratch.
			if (.not. allocated(sys_x))   allocate(sys_x(N_eq))
			if (.not. allocated(sys_sol)) allocate(sys_sol(N_eq))
			if (.not. allocated(wa))      allocate(wa(lwa))

			! THE COMPOSITION THE CELL ENTERS THE SWEEP WITH, in the
			! fraction layout of the system this run solves: the previous
			! adopted step's composition, or the initial condition at step
			! 0. Two things read it. It is the starting point of the solve
			! for every step after the first, and it is the state a cell
			! KEEPS when the only candidate the solve produced is a
			! non-root whose residual is above ieq_nonroot_res_cap
			! (acceptance class 6). Built here, once, so that the guess and
			! the kept state are one object -- and OUTSIDE the self-field
			! loop below, so that the state a cell keeps is the one it
			! entered the sweep with and not an intermediate iterate.
			x_entry(1:N_eq) = 0.0d0
			x_entry(1) = nhii(j)/nh(j)
			x_entry(2) = nheii(j)/nhe(j)
			x_entry(3) = nheiii(j)/nhe(j)
			if (thereis_HeITR) x_entry(4) = nheiTR(j)/nhe(j)
			if (thereis_mol) then
				x_entry(4) = 2.0d0*nmol_eq(j,1)/nh(j)
				x_entry(5) = 2.0d0*nmol_eq(j,2)/nh(j)
				x_entry(6) = 3.0d0*nmol_eq(j,3)/nh(j)
				x_entry(7) = nmol_eq(j,4)/nh(j)
				if (thereis_HeITR) x_entry(8) = nheiTR(j)/nhe(j)
				if (thereis_oxychem) then
					x_entry(iox)   = nox_eq(j,1)                          &
					             /max(nm_tot(j,iel_O),1.0d-30)
					x_entry(iox+1) = nox_eq(j,2)                          &
					             /max(nm_tot(j,iel_O),1.0d-30)
				endif
			endif
			if (thereis_metals) then
				do im = 1,n_melem
					i0 = melem_i0(im)
					x_entry(mbase+2*(im-1)) = nm(j,i0+1)                  &
					                  /max(nm_tot(j,im),1.0d-30)
					if (melem_top(im) .ge. 2) then
						x_entry(mbase+1+2*(im-1)) = nm(j,i0+2)            &
						                  /max(nm_tot(j,im),1.0d-30)
					else
						x_entry(mbase+1+2*(im-1)) = 0.0d0
					endif
				enddo
			endif

			! THE PASSES THIS CELL MAY TAKE. Every cell takes the
			! self-field passes the run asked for. A lower ghost whose
			! molecular partition is stated by the base handoff takes as
			! many more as its closure needs: the imposed partition and the
			! ghost's own ionization balance are one system (see the ghost
			! reservoir block below and base_ghost_closure_tol), and a cell
			! that leaves after one pass leaves it unclosed. The radiation
			! field is recomputed on the self-field passes alone, so the
			! closure passes solve the same cell against the same field.
			ghost_cell_closure = thereis_mol .and. j .le. 0 .and.         &
			                     base_h2_composition_imposed()
			n_self_max = xuv_self_field_passes
			if (ghost_cell_closure) n_self_max =                          &
			     max(n_self_max, base_ghost_closure_passes)
			ghost_closed       = .not. ghost_cell_closure
			ghost_closure_move = 0.0d0
			ghost_closure_res  = 0.0d0
			x_h2_imposed       = 0.0d0
			have_h2_secant_prev = .false.
			x_h2_secant_prev    = 0.0d0
			f_h2_secant_prev    = 0.0d0

			self_field_HHe: do it_self = 1, n_self_max

			! THE CELL'S OWN ATTENUATION AT THE COMPOSITION IT RETURNED.
			! The column outside the cell is the block's and does not move
			! here; what moves is this cell's own dtau, which the field
			! routine forms from the absorber densities handed to it (see
			! xuv_self_field_passes).  The two rates the sweep adds on top
			! of the stellar field -- the Balmer continuum of H(n=2) and the
			! locally absorbed He recombination radiation -- are whole-grid
			! quantities of the previous outer pass, so they are re-added in
			! the same order the first pass added them.  Their diagnostic
			! counters are taken on the first pass alone, where the serial
			! loop above took them.
			!
			! ON A SELF-FIELD PASS ALONE. The further passes a base ghost
			! takes to close its imposed molecular partition are passes of
			! the same cell against the SAME field: the quantity they close
			! is the pair (partition, ionization), and letting the field
			! move with them would make the closure a different iteration
			! from the one the run asked for. A run whose self-field passes
			! are the default one therefore sees exactly the field it saw.
			if (it_self .gt. 1 .and. it_self .le. xuv_self_field_passes) then
				nheiS_it = nhei(j)
				if (thereis_HeITR) nheiS_it =                                  &
				         he_ground_singlet_density(nhei(j), nheiTR(j))
				nh2_it = 0.0d0
				if (has_h2_xuv) nh2_it = nmol_eq(j,1)
				call photoionization_field_at_cell_HHe(j,                  &
				         N1_face(j),N15_face(j),N2_face(j),NTR_face(j),    &
				         NH2col_face(j),Nm_face(j,:),                      &
				         nhi(j),nheiS_it,nheii(j),nheiTR(j),nh2_it,        &
				         nm(j,:), xion(j), f_vib_quench(j),                &
				         e_vib_bound(j),                                   &
				         has_h2_xuv, sec_on_xuv, mol_sec_xuv,              &
				         D0_H2_xuv, E_ker_H2_dd_xuv,                       &
				         P_HI(j),P_HeI(j),P_HeII(j),P_HeITR(j), Pm_row,    &
				         P_H2(j),P_H2_di(j),P_H2_dd(j),P_H2_nd(j),         &
				         h1_HI(j),h1_HeI(j),h1_HeII(j),h1_HeTR(j),         &
				         h1_H2(j),                                         &
				         h1m_row, heat_row, chan_row, q(j), q_abs_row)
				P_m(j,:)  = Pm_row
				h1_m(j,:) = h1m_row
				if (use_excited_H) then
					gph_ground_HI(j) = P_HI(j)
					P_HI(j) = P_HI(j) + gph_balmer_HI(j)
				endif
				if (use_he_rec_coupling) then
					P_HI(j) = P_HI(j) + dP_HI_hrc(j)
					if (thereis_mol) P_H2(j) = P_H2(j) + dP_H2_hrc(j)
					if (thereis_metals) P_m(j,:) = P_m(j,:) + dP_m_hrc(j,:)
				endif
			endif

			! System coefficients: named-field cell state (Inc 4). The He
			! equilibrium systems read these; params is now only the MINPACK
			! transport argument (unread by the converted systems).
			ieq_cell%P_HI       = P_HI(j)
			ieq_cell%P_HeI      = P_HeI(j)
			ieq_cell%P_HeII     = P_HeII(j)
			ieq_cell%rchiiB     = rchiiB(j)
			ieq_cell%rcheiiB    = rcheiiB(j)
			ieq_cell%rcheiiiB   = rcheiiiB(j)
			ieq_cell%nh         = nh(j)
			ieq_cell%jcell      = j
			ieq_cell%nhe        = nhe(j)
			ieq_cell%a_ion_HI   = a_ion_HI(j)
			ieq_cell%a_ion_HeI  = a_ion_HeI(j)
			ieq_cell%a_ion_HeII = a_ion_HeII(j)

			! He <-> H charge-exchange rate coefficients (Huang Table 4 group
			! B): read by he_h_cx_fvec/he_h_cx_jac in every He system. Depends
			! only on T, so evaluate once per cell here (cheap). The residual
			! routines add nothing when he_h_charge_exchange is off, so this is
			! harmless (and bit-identical) in that case.
			call he_h_cx_rates(T_K(j), ieq_cell%kcx_He0_Hp,               &
			                           ieq_cell%kcx_Hep_H0)

			! The He 2^3S channels: the cell's rates where the level is
			! tracked, and zero where it is not, so that every system and
			! every operator that reads them describes the same gas. The
			! rows of the transport-chemistry operator read them out of the
			! frozen background of the step, which is kept for an atomic run
			! as well as for a molecular one.
			if (thereis_HeITR) then
				ieq_cell%rcheiTR = rcheiTR(j)
				ieq_cell%A31     = A31
				ieq_cell%P_HeITR = P_HeITR(j)
				ieq_cell%q13     = q13(j)
				ieq_cell%q31a    = q31a(j)
				ieq_cell%q31b    = q31b(j)
				ieq_cell%Q31     = Q31(j)
				ieq_cell%a_ion_HeITR = a_ion_HeITR(j)
			else
				ieq_cell%rcheiTR = 0.0d0     ! no triplet channels
				ieq_cell%A31     = 0.0d0
				ieq_cell%P_HeITR = 0.0d0
				ieq_cell%q13     = 0.0d0
				ieq_cell%q31a    = 0.0d0
				ieq_cell%q31b    = 0.0d0
				ieq_cell%Q31     = 0.0d0
				ieq_cell%a_ion_HeITR = 0.0d0
			endif

			! The temperature and the gas particle density of the cell. Both
			! are read by every operator that evaluates a chemistry row on
			! the frozen background of a step, so they belong to the cell
			! state of any gas and not to the molecular layout alone.
			ieq_cell%T_K  = T_K(j)
			ieq_cell%ntot = n_tot(j)        ! M for the 3-body rates

			! molecular cell state (System_HeH_mol layout)
			if (thereis_mol) then
				ieq_cell%P_H2    = P_H2(j)
				ieq_cell%P_H2_di = P_H2_di(j)
				ieq_cell%P_H2_dd = P_H2_dd(j)
				ieq_cell%P_H2_nd = P_H2_nd(j)
				P_H2_eq(j)    = P_H2(j)
				ieq_cell%k_LW = k_lw_diss(j)   ! 0 without a LW band flux
				! Compute the molecular rate coefficients that are invariant
				! across this cell's Newton solve (they depend only on T and
				! n_tot); the residual then reads them, like set_metal_coeffs.
				! n_tot, not n_in_dim: R13 and R15 fold the third body M into
				! their coefficient, so they need the same particle density
				! the R12 term of the residual reads (see the n_tot comment
				! above).
				call set_mol_coeffs(T_K(j), n_tot(j))
				! Oxygen-chemistry cell state: the free oxygen family the
				! carrier fractions are measured against, the frozen CO
				! reservoir, and the rate coefficients and band photolysis
				! rates of the carriers.
				if (thereis_oxychem) then
					ieq_cell%n_ofam = nm_tot(j,iel_O)
					ieq_cell%n_co   = nCO_cell(j)
					call set_oxygen_coeffs(T_K(j), j_h2o_fuv(j,:),        &
					                               j_oh_fuv(j,:))
				else
					ieq_cell%n_ofam = 0.0d0
					ieq_cell%n_co   = 0.0d0
				endif
				! Turnover scale of each balance row, from the coefficients
				! just built and this cell's densities; the residual divides
				! by it so the helium and molecular blocks reach hybrd1 with
				! the same weight.
				! photo_scale = 1: the cell's own radiation field. Only
				! the continuation of constrained_chemical_equilibrium
				! ever passes anything else, and it restores this on
				! every exit.
				call set_mol_turnover_rates(ne(j), 1.0d0)
				! The oxygen rows get theirs after that reset, and it also
				! adds the oxygen terms to the H2 row's bound.
				if (thereis_oxychem)                                  &
					call set_oxygen_turnover_rates(iox, 1.0d0)
				! Molecular partitions this cell does not own. Default: it
				! owns all of them and every row is a balance.
				ieq_cell%x_h2_fixed     = .false.
				ieq_cell%x_ox_fixed     = .false.
				ieq_cell%x_h2_fix       = 0.0d0
				ieq_cell%x_oh_fix       = 0.0d0
				ieq_cell%x_h2o_fix      = 0.0d0
				! With the carriers transported, their partition is not a
				! local root any more: the transport-chemistry solve owns
				! H2, OH and H2O and this sweep is handed the answer. The
				! remaining rows are then solved against it, which is what
				! keeps the ionization stages consistent with the
				! transported composition.
				! bg_ready OR a restart: on a cold start the very first
				! sweep is what INITIALIZES the carriers, from the local
				! equilibrium of the seeded state, and imposing a seed on
				! it would be imposing nothing. On a restart the loaded
				! state already carries them and re-solving it locally would
				! throw the transported partition away -- which is what the
				! round-trip gate G6 measures.
				if (thereis_mol .and. carrier_transport                   &
				    .and. (bg_ready .or. do_load_IC)) then
					ieq_cell%x_h2_fixed = .true.
					ieq_cell%x_h2_fix  = 2.0d0*nmol_eq(j,1)/nh(j)
					! The oxygen carriers exist only with the oxygen cycle;
					! without it there is nothing to hand this sweep.
					if (thereis_oxychem) then
						ieq_cell%x_ox_fixed = .true.
						ieq_cell%x_oh_fix  = nox_eq(j,1)                  &
						                   /max(nm_tot(j,iel_O),1.0d-30)
						ieq_cell%x_h2o_fix = nox_eq(j,2)                  &
						                   /max(nm_tot(j,iel_O),1.0d-30)
					endif
				endif
				! THE LOWER-BOUNDARY RESERVOIR. The ghost cells below the
				! base are not a piece of atmosphere this cell's radiation
				! field gets to determine: they are the gas flowing IN, and
				! its molecular partition was set in the lower atmosphere,
				! where the photodissociating field and the mixing that the
				! base cell cannot see both act. A shielded base cell left
				! to its own balance drives H2 to the fully molecular limit
				! -- measured at 0.9997 of the H nuclei against the 0.9251
				! the handoff states -- and the equation of state and the
				! species state then describe different gas at the same
				! ghost. So the handoff value is imposed here and the
				! ionization stages are solved against it.
				!
				! Both lower ghosts are pinned: they are the same reservoir,
				! and Rec_BC builds face states from each of them.
				!
				! Only when a handoff actually stated the value. Without one
				! the partition comes from the chemical-equilibrium fit,
				! which is this cell's own local estimate rather than
				! upstream information -- imposing it would be imposing the
				! solver's own answer on the solver.
				!
				! WHAT THE HANDOFF PRESCRIBES IS THE PARTITION OF THE
				! NON-IONIZED HYDROGEN. The lower-atmosphere model states
				! q_H2 for gas that its own field acts on, and that gas is
				! neutral; the ionization of this reservoir is produced by
				! the wind's field, which that model never saw, so the
				! reservoir's H+ is determined here by the ionization
				! balance and the handoff fraction x2 applies to the
				! hydrogen nuclei that are left. The row this sweep solves
				! is x(4) = 2 n(H2)/n_H over ALL the hydrogen nuclei, so
				! what is imposed on it is x2 (1 - x_ion).
				!
				! x2 OF EVERY NUCLEUS IS AN OVER-PRESCRIPTION, and one that
				! nothing checked. It leaves 1 - x2 of the hydrogen for H+,
				! H2+, H3+ and HeH+ together -- 1.7e-5 at the Koskinen 2022
				! handoff x2 = 0.99998 -- while the ionization balance of a
				! reservoir the wind has warmed asks for its own x(H+),
				! measured at 3.2e-4 there, 19 times the room. The system
				! is then infeasible: the solve breaks the pinned row to
				! pay for the ionization and the sweep reports the cell as
				! resting on a non-root. x2 (1 - x_ion) instead leaves
				! (1 - x2)(1 - x_ion) + x_ion, which is never less than the
				! ionized fraction itself, so the prescription and the
				! ionization always fit in the hydrogen the cell has.
				!
				! x_ion is the ionized hydrogen-nucleus fraction of the
				! state this sweep was handed, so the pair (partition,
				! ionization) is reached by the alternation and is exact at
				! a stationary state. The molecular ions carry hydrogen
				! nuclei too -- two in H2+, three in H3+, one in HeH+ --
				! and they are ionized, so they count here with the same
				! multiplicities the unknowns x(5), x(6), x(7) carry.
				if (base_h2_composition_imposed() .and. j .le. 0) then
					x_ion_ghost = (nhii(j) + 2.0d0*nmol_eq(j,2)           &
					            + 3.0d0*nmol_eq(j,3) + nmol_eq(j,4))      &
					            /nh(j)
					if (.not. (x_ion_ghost .ge. 0.0d0)) x_ion_ghost = 0.0d0
					if (x_ion_ghost .gt. 1.0d0) x_ion_ghost = 1.0d0
					ieq_cell%x_h2_fixed = .true.
					! THE FIXED POINT IS REACHED BY A SECANT STEP, NOT BY
					! SUBSTITUTION. Imposing x2 (1 - x_ion) at the x_ion
					! the previous pass returned is the map
					! x -> g(x) = x2 (1 - x_ion(x)), and its slope
					! g' = -x2 dx_ion/dx is positive (more H2 at fixed T
					! leaves less ionized hydrogen): measured 0.74 on
					! the ghost of the LHS 1140 b molecular seed at
					! He/H = 2.13 and 0.90 at He/H = 9.7, a monotone
					! contraction that substitution would need about 60
					! and 180 passes to take to 1e-10 (extrapolated
					! from those ratios at pass 30). The root of
					! F(x) = g(x) - x is a scalar equation, and the
					! secant through the last two evaluations solves it
					! superlinearly. It is used only when both
					! evaluations were made in the same radiation field
					! (it_self - 1 past the self-field passes), only
					! where F falls through the pair (the slope of F is
					! g' - 1 < 0, so the root is on the side the two
					! residuals point to) and only when its point is a
					! partition the prescription can take, 0 <= x <= x2;
					! otherwise the pass substitutes. The room the
					! solve needs is not checked against the x_ion of
					! the previous pass: at the secant point the linear
					! model of g gives x_ion = 1 - x/x2, and x plus that
					! is 1 - x (1 - x2)/x2 <= 1 for every x2 <= 1.
					x_h2_map = base_h2_nuclei_fraction()                  &
					         *(1.0d0 - x_ion_ghost)
					ieq_cell%x_h2_fix = x_h2_map
					if (it_self .gt. 1) then
						f_h2_map = x_h2_map - x_h2_imposed
						if (have_h2_secant_prev .and.                     &
						    (f_h2_map - f_h2_secant_prev)                 &
						    *(x_h2_imposed - x_h2_secant_prev)            &
						    .lt. 0.0d0) then
							x_h2_secant = x_h2_imposed - f_h2_map         &
							     *(x_h2_imposed - x_h2_secant_prev)       &
							     /(f_h2_map - f_h2_secant_prev)
							if (x_h2_secant .ge. 0.0d0 .and.              &
							    x_h2_secant .le.                          &
							         base_h2_nuclei_fraction())           &
								ieq_cell%x_h2_fix = x_h2_secant
						endif
						if (it_self - 1 .ge. xuv_self_field_passes) then
							have_h2_secant_prev = .true.
							x_h2_secant_prev    = x_h2_imposed
							f_h2_secant_prev    = f_h2_map
						endif
					endif
					! The partition this pass imposes, kept so that the
					! closure below can measure what the solve made of it
					! and how far one more pass would move it.
					ghost_closure_move = abs(ieq_cell%x_h2_fix            &
					                         - x_h2_imposed)
					x_h2_imposed       = ieq_cell%x_h2_fix
					! AND THE PRESCRIPTION IS CHECKED AGAINST THE CELL'S
					! OWN HYDROGEN, which is what the over-prescription
					! above never was. The room the pinned row leaves is
					! 1 - x_h2_fix and the hydrogen already held in ions is
					! x_ion; a prescription that leaves less room than that
					! cannot be satisfied by any composition, and is
					! refused here by name instead of reaching the solver
					! and being reported as a streak of non-roots.
					if (x_h2_map + x_ion_ghost .gt. 1.0d0) then
						!$omp critical (ieq_acc_report)
						write(*,'(A)') ' (ioniz_eq) STOP: the base '//    &
							'handoff prescribes more H2 than the '//      &
							'reservoir has non-ionized hydrogen'
						write(*,'(A,I0,A,ES12.5)') '   ghost cell ', j,   &
							'   imposed 2 n(H2)/n_H ', x_h2_map
						write(*,'(A,ES12.5,A,ES12.5)')                    &
							'   ionized hydrogen fraction it must make '//&
							'room for ', x_ion_ghost,                     &
							'   room left by the prescription ',          &
							1.0d0 - x_h2_map
						write(*,'(A,ES12.5)')                             &
							'   handoff x2 = base_h2_nuclei_fraction() ', &
							base_h2_nuclei_fraction()
						flush(6)
						!$omp end critical (ieq_acc_report)
						error stop 'ioniz_eq: infeasible base H2 handoff'
					endif
				endif
			endif

			! Ionization stages this cell does not own. Default: it owns
			! all three and every stage row is a balance.
			ieq_cell%x_hp_fixed     = .false.
			ieq_cell%x_heii_fixed   = .false.
			ieq_cell%x_heiii_fixed  = .false.
			ieq_cell%x_hp_fix       = 0.0d0
			ieq_cell%x_heii_fix     = 0.0d0
			ieq_cell%x_heiii_fix    = 0.0d0
			! THE IONIZATION STAGES THE FLOW CARRIES, on their own
			! gate. They are stages of an element, not molecular
			! carriers: their transport is the element nucleus flux of
			! hydrogen and of helium, and the condition for imposing
			! them is the key and the availability of a transported
			! state, not the molecular carrier configuration they
			! share an operator with, and the same statement holds
			! of the gas: a hydrogen and helium mixture has these
			! stages whether or not it has molecules, so the block
			! stands beside the molecular cell state and not
			! inside it.
			! nhii, nheii and nheiii here are f_sp*n as this sweep was
			! handed them, i.e. what the transport operator's
			! write-back left, so the fractions imposed are the
			! transported ones.
			! ONLY THE CELLS THE OPERATOR OWNS. The transport
			! write-back fills 1..N+Ng and leaves the two lower
			! ghosts alone -- they are the inflow reservoir, and
			! the only carrier anything below the base states a
			! partition for there is H2, through the molecular
			! handoff (carrier_base_composition_imposed); no
			! ionized fraction is one of them, so this operator
			! states none for the reservoir. Imposing the
			! transported value there would pin the ghosts at
			! whatever the last sweep left and never let them be
			! re-solved again.
			if (ionization_transport .and. j .ge. 1                   &
			    .and. (bg_ready .or. do_load_IC)) then
				ieq_cell%x_hp_fixed = .true.
				ieq_cell%x_hp_fix   = nhii(j)/nh(j)
				! The two ionized helium stages, per helium NUCLEUS,
				! which is the variable their rows are written in
				! (rows 2 and 3 of every system that reaches them).
				if (thereis_He) then
					ieq_cell%x_heii_fixed  = .true.
					ieq_cell%x_heiii_fixed = .true.
					ieq_cell%x_heii_fix    = nheii(j) /nhe(j)
					ieq_cell%x_heiii_fix   = nheiii(j)/nhe(j)
				endif
			endif

			! Keep this cell's coefficient state for the transport-chemistry
			! operator, which evaluates the same rows on the same background
			! between sweeps. Each iteration writes its own element, so the
			! parallel sweep needs no guard. It is kept wherever that
			! operator has a row to solve, which in an atomic gas is the
			! three ionization stages alone.
			if (transported_rows_exist()) bg_cell(j) = ieq_cell

			! Each element's metal coefficients are handed to
			! ion_system_HeH_metals via set_metal_coeffs; the charge-
			! exchange rate coefficients are stored for this cell by
			! cx_set_cell (used inside the residual by cx_add_to_fvec).
			if (thereis_metals) then
				call cx_set_cell(T_K(j))

				! Build the metal coefficients for each element in canonical order
				! from the 2D rate arrays (col i0 = neutral, i0+1 = +,
				! i0+2 = ++) and store them for the residual.
				do im = 1,n_melem
					i0  = melem_i0(im)
					top = melem_top(im)
					meg_top(im)  = top
					meg_ntot(im) = nm_tot(j,im)
					meg_g0(im)   = P_m(j,i0)
					meg_b0(im)   = aion_m(j,i0)
					meg_a1(im)   = rec_m(j,i0+1)
					if (top .ge. 2) then
						meg_g1(im) = P_m(j,i0+1)
						meg_b1(im) = aion_m(j,i0+1)
						meg_a2(im) = rec_m(j,i0+2)
					else
						! No second ionization stage for this element.
						meg_g1(im) = 0.0d0
						meg_b1(im) = 0.0d0
						meg_a2(im) = 0.0d0
					endif
				enddo
				call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1, &
				                    meg_b0, meg_b1, meg_a1, meg_a2,     &
				                    meg_top)

				! Metal rows of the molecular system get their turnover
				! scale too, after set_mol_turnover_rates has reset the
				! array and set_metal_coeffs has filled met_*.
				if (thereis_mol)                                     &
					call set_mol_metal_turnover_rates(ne(j), 1.0d0)
			endif

			! The cell state this solve is about to run on, kept for the
			! closure evaluator (see ieq_rate_cell). Written here, never
			! read by the sweep.
			ieq_rate_cell(j) = ieq_cell
			ieq_ne_cell(j)   = ne(j)
			ieq_TK_cell(j)   = T_K(j)
			ieq_ntot_cell(j) = n_tot(j)
			if (thereis_metals) then
				ieq_met_coef(:,1,j) = meg_ntot
				ieq_met_coef(:,2,j) = meg_g0
				ieq_met_coef(:,3,j) = meg_g1
				ieq_met_coef(:,4,j) = meg_b0
				ieq_met_coef(:,5,j) = meg_b1
				ieq_met_coef(:,6,j) = meg_a1
				ieq_met_coef(:,7,j) = meg_a2
			endif

			! Initial guess.  The second and later passes of the self-field
			! loop start from the state the previous one accepted, which is
			! the composition their field was built from.
			if (it_self .gt. 1) then
				sys_x(1:N_eq) = x_self(1:N_eq)
			else if (marching_step .eq. 0) then
				if (j .eq. N+Ng) then
					sys_x(1) = 1.0
					sys_x(2) = 1.0
					sys_x(3) = 1.0
					if (thereis_HeITR) sys_x(4) = 0.01
					if (thereis_mol) then
						! start fully ionized aloft; molecules negligible
						sys_x(4:7) = 1.0d-10
						if (thereis_HeITR) sys_x(8) = 0.01
						! Oxygen carriers likewise negligible aloft: the
						! oxygen there is atomic and ionized.
						if (thereis_oxychem) then
							sys_x(iox)   = 1.0d-12
							sys_x(iox+1) = 1.0d-12
						endif
					endif
					if (thereis_metals) then
						! Each metal starts fully singly ionized.
						do im = 1,n_melem
							sys_x(mbase + 2*(im-1))   = 1.0  ! X+  frac
							sys_x(mbase+1 + 2*(im-1)) = 0.0  ! X++ frac
						enddo
					endif
				else
					sys_x(1) = nhii(j+1)/nh(j+1)
					sys_x(2) = nheii(j+1)/nhe(j+1)
					sys_x(3) = nheiii(j+1)/nhe(j+1)
					if (thereis_HeITR) &
						sys_x(4) = nheiTR(j+1)/nhe(j+1)
					if (thereis_mol) then
						sys_x(4) = max(2.0d0*nmol_eq(j+1,1)/nh(j+1), 1.0d-10)
						sys_x(5) = max(2.0d0*nmol_eq(j+1,2)/nh(j+1), 1.0d-12)
						sys_x(6) = max(3.0d0*nmol_eq(j+1,3)/nh(j+1), 1.0d-12)
						sys_x(7) = max(nmol_eq(j+1,4)/nh(j+1), 1.0d-14)
						if (thereis_HeITR) sys_x(8) = nheiTR(j+1)/nhe(j+1)
						if (thereis_oxychem) then
							sys_x(iox)   = max(nox_eq(j+1,1)              &
							         /max(nm_tot(j+1,iel_O),1.0d-30),     &
							                   1.0d-14)
							sys_x(iox+1) = max(nox_eq(j+1,2)              &
							         /max(nm_tot(j+1,iel_O),1.0d-30),     &
							                   1.0d-14)
						endif
					endif
					if (thereis_metals) then
						do im = 1,n_melem
							i0 = melem_i0(im)
							sys_x(mbase+2*(im-1)) = nm(j+1,i0+1)        &
							                  /max(nm_tot(j+1,im),1.0d-30)
							if (melem_top(im) .ge. 2) then
								sys_x(mbase+1+2*(im-1)) = nm(j+1,i0+2)        &
								                  /max(nm_tot(j+1,im),1.0d-30)
							else
								sys_x(mbase+1+2*(im-1)) = 0.0d0
							endif
						enddo
					endif
				endif
			else
				sys_x(1:N_eq) = x_entry(1:N_eq)

				! The guess just built is this cell's previous state. If that
				! state is not physical it must not seed the solve again: an
				! unphysical root fed back as the next step's guess reproduces
				! itself indefinitely (the self-sticking negative H II / Fe I
				! state at the base). Start from the local ionization balance
				! instead, which lies inside the simplex by construction. The
				! molecular network keeps its own chemical-equilibrium retry.
				if (.not.thereis_mol) then
					if (.not.ionization_fractions_physical(sys_x,N_eq,mbase)) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_reseed = n_ieq_reseed + 1
					endif
				endif
			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			! The He metastable-triplet systems keep the MINPACK solve (no
			! analytic Jacobian written for the triplet kinetics).
			if (thereis_mol) then
				! Molecular network, with the trace metals solved in the same
				! system when they are present (System_HeH_mol_metals: the two
				! blocks share the free electron density, which the metals
				! dominate in the shielded molecular base).
				!
				! hybrd1 is bistable in its initial guess, because the cell
				! itself is: a molecular basin (the dense, optically-thick
				! base, where the true root is strongly molecular) and an
				! atomic basin (the wind above the H2 -> H front). From a
				! zero-molecular warm start the solve fails at the base; at
				! the front cells themselves neither warm start need land in
				! the right basin. One starting point in each basin is
				! therefore tried in turn:
				!   1  the guess built above (previous state / neighbor),
				!   2  molecular basin -- the H2 dissociation equilibrium of
				!      the local (p, T) with every element in its own
				!      ionization balance at the incoming n_e,
				!   3  atomic basin -- the same balance with no H2,
				! and an iterate is accepted as a ROOT only if it is physical
				! (every stage fraction >= 0, each element's tracked stages
				! summing to at most its nuclei) AND its normalized reaction
				! residual is within ieq_res_tol. The solver exit code is
				! neither sufficient nor necessary for that judgement (see
				! ieq_res_tol above).
				!
				! A cell that finds such a root takes it (a solver-converged
				! one immediately). A cell that does not keeps the first
				! bounded iterate as a MARKED non-root under the relaxation
				! amnesty (nonroot_streak_update), and when NONE of the
				! attempts even stayed inside the simplex it clamps the one
				! lying closest onto the element budget and RECHECKS the
				! reaction residual, which the projection changes. Every
				! outcome is computed at THIS cell's state; none of them is
				! the cell's previous composition, which would make the
				! equilibrium -- and the steady residual built on it -- a
				! function of the sequence of evaluations rather than of the
				! state being evaluated.
				best_rank = 0
				viol_best = 0.0d0
				res_best  = 0.0d0
				info_best = 0
				have_nonroot = .false.
				have_clamp   = .false.
				! Gas pressure p = (n_tot + n_e) kB T [bar], the same ideal-gas
				! law the EOS uses (comp_p_from_T).
				pbar_loc  = (n_tot(j) + ne(j))*kb_erg*T_K(j)/1.0d6
				do iatt = 1,3
					if (iatt .eq. 2) then
						call dissociation_ionization_balance_at_fixed_ne(  &
						                 sys_x,N_eq,mbase,ne(j),pbar_loc,T_K(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 3) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					endif

					if (thereis_metals) then
						! Metals appended above the molecular unknowns; point
						! charge exchange at their rows for this solve.
						cx_metal_base = mbase
						call hybrd1(ion_system_HeH_mol_metals,N_eq,sys_x,  &
						            sys_sol,tol,info,wa,lwa,params)
						cx_metal_base = 4
					else
						call hybrd1(ion_system_HeH_mol,N_eq,sys_x,sys_sol, &
						            tol,info,wa,lwa,params)
					endif
					conv_ieq = (info .eq. 1)
					n_mol_info(max(0,min(5,info))) =                     &
						n_mol_info(max(0,min(5,info))) + 1

					phys_ieq = ionization_fractions_physical(sys_x,N_eq,mbase)
					if (iatt .eq. 1 .and. .not.phys_ieq)                &
						n_ieq_unphys = n_ieq_unphys + 1
					! Reaction residual of a physical iterate (section 113):
					! the histograms measure where converged roots sit and
					! where the bounded-but-unconverged iterates do.
					res_att = 0.0d0
					if (phys_ieq) then
						res_att = normalized_reaction_residual(sys_x,     &
						                              N_eq,mbase,ne(j))
						if (conv_ieq) then
							hist_conv(residual_decade(res_att)) =         &
								hist_conv(residual_decade(res_att)) + 1
						else
							hist_uncv(residual_decade(res_att)) =         &
								hist_uncv(residual_decade(res_att)) + 1
						endif
					endif
					ok_rank = 0
					if (phys_ieq .and. res_att .le. ieq_res_tol) ok_rank = 1
					if (ok_rank .eq. 1 .and. conv_ieq) ok_rank = 2
					! Strict improvement only, so attempt 1 wins any tie.
					if (ok_rank .gt. best_rank) then
						best_rank = ok_rank
						x_root_best(1:N_eq) = sys_x(1:N_eq)
						res_best  = res_att
						info_best = info
					else if (best_rank .eq. 0) then
						if (phys_ieq) then
							! A bounded iterate that is not a root: the state
							! of last resort, first one kept (the old rank-1
							! order).
							if (.not. have_nonroot) then
								have_nonroot = .true.
								x_nonroot(1:N_eq) = sys_x(1:N_eq)
								res_nonroot  = res_att
								info_nonroot = info
							endif
						else
							! Outside the simplex: remember the iterate that
							! lies closest, to be clamped onto it below.
							viol = element_budget_violation(sys_x,N_eq,mbase)
							if (.not. have_clamp .or. viol .lt. viol_best) then
								have_clamp = .true.
								viol_best  = viol
								x_root_best(1:N_eq) = sys_x(1:N_eq)
							endif
						endif
					endif
					if (ok_rank .eq. 2) exit
				enddo

				! No starting point produced a root. Before the degraded
				! paths below -- the marked non-root of last resort and the
				! projection onto the element budget -- solve the CONSTRAINED
				! problem: positive species densities as unknowns, element
				! conservation as explicit rows (so that HeH+ and the oxygen
				! carriers satisfy BOTH element budgets they belong to at
				! once, which a projection cannot do), charge neutrality
				! identically, and a continuation in the radiation field from
				! the dense molecular limit instead of three unconnected
				! guesses (docs/supersonic_molecular_base.md section 11.5-B).
				!
				! Its result is a candidate, not an acceptance: it faces the
				! SAME judge every other candidate faces, bounds and the
				! normalized reaction residual against ieq_res_tol. Reaching
				! the full field is a statement about the continuation, not
				! about the cell's network.
				!
				! Only the molecular branch is promoted. The atomic matrix
				! cases carry no measured class-3/4 population (section 113),
				! and the constrained formulation is aimed at what the
				! molecular network has and the atomic one does not: species
				! that carry two elements at once, whose two budgets a
				! fraction layout can only satisfy one at a time.
				cce_root = .false.
				res_cce  = 0.0d0
				if (best_rank .eq. 0) then
					n_cce_attempt = n_cce_attempt + 1
					x_cce(1:N_eq) = 0.0d0
					call system_clock(count=clk_beg, count_rate=clk_rate)
					call equilibrium_from_molecular_limit(x_cce,N_eq,  &
					              mbase,iox,ne(j),pbar_loc,cce_full,   &
					              n_cce_fs)
					call system_clock(count=clk_end)
					if (clk_rate .gt. 0)                               &
						cce_seconds = cce_seconds                      &
						  + dble(clk_end - clk_beg)/dble(clk_rate)
					n_cce_solve = n_cce_solve + n_cce_fs
					if (cce_full) then
						phys_ieq =                                     &
						  ionization_fractions_physical(x_cce,N_eq,mbase)
						if (.not.phys_ieq) then
							! The conservation rows are residuals, so
							! they leave at most tolerance-level slack;
							! removing it is a round-off correction of
							! the solve, and the residual is measured
							! AFTER it, which is the section 11.5-A
							! recheck rule. A larger violation is not
							! round-off and is not clamped away.
							if (element_budget_violation(x_cce,N_eq,   &
							             mbase) .le. 1.0d-6) then
								call clamp_fractions_to_element_budget(&
								             x_cce,N_eq,mbase)
								phys_ieq =                             &
								  ionization_fractions_physical(x_cce, &
								                        N_eq,mbase)
							endif
						endif
						if (phys_ieq) then
							res_cce =                                  &
							  normalized_reaction_residual(x_cce,N_eq, &
							                        mbase,ne(j))
							if (res_cce .le. ieq_res_tol)              &
								cce_root = .true.
						endif
					endif
					if (cce_root) n_cce_root = n_cce_root + 1
				endif

				! Accept the state (section 113). Classes: 1 = solver-
				! converged root, 2 = root without solver convergence,
				! 3 = projected state whose RECHECKED residual still marks a
				! root, 4 = non-root under the relaxation amnesty
				! (nonroot_streak_update reports it and stops the run if it
				! persists), 5 = root of the constrained element-conserving
				! continuation solve just above, 6 = non-root whose residual
				! is above ieq_nonroot_res_cap, for which the cell keeps the
				! composition it entered the sweep with.
				! THE CONTRACT EVERY CONSUMER READS: classes 1, 2, 3 and 5
				! are ROOTS of the requested equations -- each one passed the
				! same two tests, a state inside the element simplex and a
				! normalized reaction residual at or below ieq_res_tol -- and
				! differ only in how the root was found. Classes 4 and 6 are
				! the classes that are NOT roots. A consumer counting states
				! its chemistry does not describe counts 4 and 6 and nothing
				! else (n_cells_without_chemical_root, steady_residual.f90).
				if (best_rank .eq. 2) then
					acc_class = 1
					acc_res   = res_att   ! sys_x is the attempt just measured
				else if (best_rank .eq. 1) then
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					acc_class = 2
					acc_res   = res_best
					info      = info_best
				else if (cce_root) then
					sys_x(1:N_eq) = x_cce(1:N_eq)
					acc_class = 5
					acc_res   = res_cce
					! info stays the exit code of the last fraction-layout
					! attempt: it is what the promotion followed, and the
					! continuation's own status is not a solver code (its
					! cost is in the end-of-run constrained-solve line).
					call report_acceptance_event(                      &
						'constrained-continuation root accepted',      &
						j,marching_step,r(j),T_K(j),info,viol_best,acc_res)
				else if (have_nonroot) then
					sys_x(1:N_eq) = x_nonroot(1:N_eq)
					acc_class = 4
					acc_res   = res_nonroot
					info      = info_nonroot
				else
					! Every attempt left the physical simplex: project the
					! closest one onto the element budget -- the violations
					! are the solver's own resolution of a root sitting on a
					! face (a fully dissociated or fully ionized species) --
					! and recheck the reaction residual, which the nonlinear
					! projection changes.
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					call clamp_fractions_to_element_budget(sys_x,N_eq,mbase)
					n_mol_clamped = n_mol_clamped + 1
					acc_res   = normalized_reaction_residual(sys_x,N_eq,   &
					                                         mbase,ne(j))
					acc_class = 3
					if (acc_res .gt. ieq_res_tol) acc_class = 4
				endif
				! A residual that is not a finite number is not a small
				! one: the state has driven a balance row of the network
				! out of the reals, and no comparison against a tolerance
				! is meaningful for it (NaN fails every one of them, so the
				! class-3 recheck above would have let it pass as a root).
				! It is a non-root, it is counted as the sweep's
				! admissibility signal, and it is kept out of the
				! largest-residual statistic, which a single NaN would
				! otherwise destroy for the whole run.
				if (.not. finite_real(acc_res)) then
					acc_class = 4
					n_res_nonfinite = n_res_nonfinite + 1
				endif
				! THE SEVERITY BOUND OF THE AMNESTY
				! (nonroot_acceptance_class). A candidate too far from
				! solving the network is not adopted: the cell keeps
				! the composition it entered the sweep with, a state
				! the run already held, and the heating and cooling
				! built after the sweep are those of THAT composition
				! (docs/Update_EXHALE_stage1.md section 170).
				if (acc_class .eq. 4) then
					acc_class = nonroot_acceptance_class(acc_res)
					if (acc_class .eq. 6)                                  &
						sys_x(1:N_eq) = x_entry(1:N_eq)
				endif
				! Has the cell's own field stopped moving its own composition
				! (xuv_self_field_passes).  The comparison is between two
				! ACCEPTED states of the same cell, in the stage fractions the
				! coupled loop of EXHALE_main measures its own convergence in.
				self_moved = 1.0d0
				if (it_self .gt. 1)                                        &
					self_moved = maxval(abs(sys_x(1:N_eq)                  &
					                        - x_self(1:N_eq)))
				x_self(1:N_eq) = sys_x(1:N_eq)
				! THE BASE HANDOFF'S OWN EQUATION, AT THE STATE THIS PASS
				! RETURNED. The ghost owes x_H2 = x2 (1 - x_ion) with the
				! x_ion of its own ionization balance, and the residual of
				! that equation is measured on the solved state: x(4) is
				! the solved partition and x(1), x(5), x(6), x(7) hold the
				! hydrogen nuclei that are ionized, with the multiplicities
				! their rows carry. The pair is closed when that residual
				! is below the tolerance.
				if (ghost_cell_closure) then
					x_ion_ghost = sys_x(1) + sys_x(5) + sys_x(6)          &
					            + sys_x(7)
					if (.not. (x_ion_ghost .ge. 0.0d0)) x_ion_ghost = 0.0d0
					if (x_ion_ghost .gt. 1.0d0) x_ion_ghost = 1.0d0
					ghost_closure_res = abs(sys_x(4)                      &
					     - base_h2_nuclei_fraction()*(1.0d0 - x_ion_ghost))
					ghost_closed = (ghost_closure_res .le.                &
					                base_ghost_closure_tol)
				endif
				last_self = (it_self .ge. n_self_max)                      &
				            .or. (self_moved .le. xuv_self_field_tol       &
				                  .and. ghost_closed)
				! WHICH PASS ADVANCES THE NON-ROOT STREAK. The streak
				! counts CONSECUTIVE SWEEPS a cell rests on a non-root, and
				! a sweep under the ghost composition fixed point applies
				! the composition map more than once
				! (ghost_composition_fixed_point_tol), so it is advanced by
				! the first application alone. Every other ledger row is
				! reset at each application and therefore describes the
				! application whose state the sweep returns.
				state_is_recorded = last_self .and. sweep_advances_the_streak
				! A BOUNDARY NOTHING STANDS BEHIND IS REFUSED. The ghost
				! composition is the state the lower boundary is built on,
				! and an unclosed pair leaves it an arbitrary iterate of
				! the alternation rather than the solution of a stated
				! equation.
				if (ghost_cell_closure .and. .not. ghost_closed .and.     &
				    it_self .ge. n_self_max) then
					!$omp critical (ieq_acc_report)
					write(*,'(A)') ' (ioniz_eq) STOP: the base handoff'// &
						' partition of a lower ghost did not close'//     &
						' against its own ionization balance'
					write(*,'(A,I0,A,I0,A)') '   ghost cell ', j,         &
						'   passes ', it_self, ''
					write(*,'(A,ES12.5,A,ES12.5)') '   last move of the'//&
						' imposed 2 n(H2)/n_H ', ghost_closure_move,      &
						'   tolerance ', base_ghost_closure_tol
					write(*,'(A,ES12.5)') '   residual of x_H2 -'//       &
						' x2 (1 - x_ion) at the returned state ',         &
						ghost_closure_res
					flush(6)
					!$omp end critical (ieq_acc_report)
					error stop 'ioniz_eq: base ghost H2 closure failed'
				endif
				if (ghost_cell_closure .and. last_self) then
					ghost_closure_res_sweep =                             &
					     max(ghost_closure_res_sweep, ghost_closure_res)
					ghost_closure_pass_sweep =                            &
					     max(ghost_closure_pass_sweep, it_self)
					ghost_closure_res_cell(j)  = ghost_closure_res
					ghost_closure_pass_cell(j) = it_self
				endif

				! THE LEDGER RECORDS THE STATE THE CELL ACCEPTED, and not
				! every iterate on the way to it: the acceptance class, its
				! largest residual, the element-budget excursion and the
				! consecutive-sweep non-root streak all describe the
				! composition this cell leaves the sweep with.  The
				! counters of each attempt above (retries, exit codes,
				! residual histograms) count solves, and count every one of
				! them.
				if (last_self) then
					n_acc(acc_class) = n_acc(acc_class) + 1
					if (finite_real(acc_res))                              &
						acc_resmax(acc_class) =                            &
							max(acc_resmax(acc_class),acc_res)
					viol_sweep_worst = max(viol_sweep_worst, viol_best)
					if (state_is_recorded)                                 &
					call nonroot_streak_update(acc_class,j,marching_step,r(j),      &
					                           T_K(j),n_in_dim(j),ne(j),    &
					                           info,viol_best,acc_res,      &
					                           sys_x,N_eq)
					if (j .le. 0) ghost_acc_res_cell(j) = acc_res
				endif
				if (last_self .and. j .eq. ieq_report_cell)                &
					call report_accepted_cell_state(j, acc_class, acc_res, &
					                                info, sys_x, N_eq)
				! The chemical decay rates of the accepted state, for the
				! cells a diagnostic named: the coefficients of this cell
				! are live only here.
				if (last_self .and. ieq_decay_ncell .gt. 0)                 &
					call measure_molecular_decay_rates(j, sys_x, N_eq,     &
					                                   mbase, ne(j))
			else
				! Atomic H/He (+ He 2^3S) (+ metals). Acceptance as a root
				! requires the same two conditions as the molecular block:
				! a physical state (every stage fraction non-negative, the
				! ionized stages of each element summing to at most its
				! nucleus total) AND a normalized reaction residual within
				! ieq_res_tol -- a vanishing solver status is neither
				! sufficient nor necessary. The cell is solved from up to
				! three starting points and the best root is kept:
				!   1  the guess built above (previous state / neighbor),
				!   2  the uncoupled ionization balance at the incoming n_e,
				!   3  the optically thick limit, all nuclei neutral.
				! Attempt 1 is accepted as it stands whenever it converges to
				! a physical root within tolerance, so a healthy cell takes
				! exactly the same solver path as before.
				best_rank = 0
				res_best  = 0.0d0
				info_best = 0
				viol_best = 0.0d0
				have_nonroot = .false.
				do iatt = 1,3
					if (iatt .eq. 2) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 3) then
						sys_x(1:N_eq) = 0.0d0
					endif

					if (thereis_HeITR .and. thereis_metals) then
						! Merged He-triplet + metals: triplet at sys_x(4),
						! metals at sys_x(5..). The cell-by-cell metal
						! coefficients (set_metal_coeffs) and charge-exchange
						! rates (cx_set_cell) were already loaded above in the
						! thereis_metals block; here we only point charge
						! exchange at the shifted metal rows for the merged
						! solve.
						cx_metal_base = 5
						call hybrd1(ion_system_HeH_TR_metals,N_eq,sys_x,  &
						            sys_sol,tol,info,wa,lwa,params)
						cx_metal_base = 4
						conv_ieq = (info .eq. 1)
						info_ieq = info
					else if (thereis_HeITR) then
						call hybrd1(ion_system_HeH_TR,N_eq,sys_x,sys_sol, &
						            tol,info,wa,lwa,params)
						conv_ieq = (info .eq. 1)
						info_ieq = info
					else if (thereis_metals) then
						call solve_ieq(ion_system_HeH_metals,             &
						            jac_system_HeH_metals,                &
						            N_eq,sys_x,params,tol,wa,lwa,usednt,  &
						            info_ieq)
						conv_ieq = (info_ieq .eq. 1)
					else
						call solve_ieq(ion_system_HeH,jac_system_HeH,     &
						            N_eq,sys_x,params,tol,wa,lwa,usednt,  &
						            info_ieq)
						conv_ieq = (info_ieq .eq. 1)
					endif

					phys_ieq = ionization_fractions_physical(sys_x,N_eq,mbase)
					if (iatt .eq. 1 .and. .not.phys_ieq)                &
						n_ieq_unphys = n_ieq_unphys + 1
					! Reaction residual of a physical iterate (section 113);
					! same convention as the molecular block above.
					res_att = 0.0d0
					if (phys_ieq) then
						res_att = normalized_reaction_residual(sys_x,     &
						                              N_eq,mbase,ne(j))
						if (conv_ieq) then
							hist_conv(residual_decade(res_att)) =         &
								hist_conv(residual_decade(res_att)) + 1
						else
							hist_uncv(residual_decade(res_att)) =         &
								hist_uncv(residual_decade(res_att)) + 1
						endif
					endif
					ok_rank = 0
					if (phys_ieq .and. res_att .le. ieq_res_tol) ok_rank = 1
					if (ok_rank .eq. 1 .and. conv_ieq) ok_rank = 2
					! Strict improvement only, so attempt 1 wins any tie and
					! an already healthy cell is untouched.
					if (ok_rank .gt. best_rank) then
						best_rank = ok_rank
						x_root_best(1:N_eq) = sys_x(1:N_eq)
						res_best  = res_att
						info_best = info_ieq
					else if (best_rank .eq. 0 .and. phys_ieq              &
					         .and. .not. have_nonroot) then
						! A bounded iterate that is not a root: the state of
						! last resort, first one kept (the old rank-1 order).
						have_nonroot = .true.
						x_nonroot(1:N_eq) = sys_x(1:N_eq)
						res_nonroot  = res_att
						info_nonroot = info_ieq
					endif
					if (ok_rank .eq. 2) exit
				enddo

				! Accept the state; same classes as the molecular block. A
				! handback state's residual is measured under the FULL
				! coupled system -- the quantity that decides whether the
				! uncoupled balance is a root of it at all.
				if (best_rank .eq. 2) then
					acc_class = 1
					acc_res   = res_att   ! sys_x is the attempt just measured
				else if (best_rank .eq. 1) then
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					acc_class = 2
					acc_res   = res_best
					info_ieq  = info_best
				else if (have_nonroot) then
					sys_x(1:N_eq) = x_nonroot(1:N_eq)
					acc_class = 4
					acc_res   = res_nonroot
					info_ieq  = info_nonroot
				else
					! Every attempt left the physical simplex. Rather than
					! propagate negative densities, hand back the uncoupled
					! ionization balance, which is admissible by
					! construction, and recheck its residual under the
					! coupled system.
					! How far outside the simplex the rejected iterate
					! lay, measured before it is replaced: the size of the
					! excursion is what a caller probing this state needs
					! to know, and it was previously reported as zero.
					viol_best = element_budget_violation(sys_x,N_eq,mbase)
					call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
					                                    mbase,ne(j))
					n_ieq_fail = n_ieq_fail + 1
					acc_res   = normalized_reaction_residual(sys_x,N_eq,   &
					                                         mbase,ne(j))
					acc_class = 3
					if (acc_res .gt. ieq_res_tol) acc_class = 4
				endif
				! Same rule as the molecular branch above: a non-finite
				! residual is a non-root and stays out of the statistic.
				if (.not. finite_real(acc_res)) then
					acc_class = 4
					n_res_nonfinite = n_res_nonfinite + 1
				endif
				! THE SEVERITY BOUND OF THE AMNESTY
				! (nonroot_acceptance_class). A candidate too far from
				! solving the network is not adopted: the cell keeps
				! the composition it entered the sweep with, a state
				! the run already held, and the heating and cooling
				! built after the sweep are those of THAT composition
				! (docs/Update_EXHALE_stage1.md section 170).
				if (acc_class .eq. 4) then
					acc_class = nonroot_acceptance_class(acc_res)
					if (acc_class .eq. 6)                                  &
						sys_x(1:N_eq) = x_entry(1:N_eq)
				endif
				! Has the cell's own field stopped moving its own composition
				! (xuv_self_field_passes).  The comparison is between two
				! ACCEPTED states of the same cell, in the stage fractions the
				! coupled loop of EXHALE_main measures its own convergence in.
				self_moved = 1.0d0
				if (it_self .gt. 1)                                        &
					self_moved = maxval(abs(sys_x(1:N_eq)                  &
					                        - x_self(1:N_eq)))
				x_self(1:N_eq) = sys_x(1:N_eq)
				last_self = (it_self .ge. n_self_max)                      &
				            .or. (self_moved .le. xuv_self_field_tol)
				! Which pass advances the non-root streak, as in the
				! molecular branch above.
				state_is_recorded = last_self .and. sweep_advances_the_streak

				! THE LEDGER RECORDS THE STATE THE CELL ACCEPTED, and not
				! every iterate on the way to it: the acceptance class, its
				! largest residual, the element-budget excursion and the
				! consecutive-sweep non-root streak all describe the
				! composition this cell leaves the sweep with.  The
				! counters of each attempt above (retries, exit codes,
				! residual histograms) count solves, and count every one of
				! them.
				if (last_self) then
					n_acc(acc_class) = n_acc(acc_class) + 1
					if (finite_real(acc_res))                              &
						acc_resmax(acc_class) =                            &
							max(acc_resmax(acc_class),acc_res)
					viol_sweep_worst = max(viol_sweep_worst, viol_best)
					if (state_is_recorded)                                 &
					call nonroot_streak_update(acc_class,j,marching_step,r(j),      &
					                           T_K(j),n_in_dim(j),ne(j),    &
					                           info_ieq,viol_best,acc_res,  &
					                           sys_x,N_eq)
					if (j .le. 0) ghost_acc_res_cell(j) = acc_res
				endif
				if (last_self .and. j .eq. ieq_report_cell)                &
					call report_accepted_cell_state(j, acc_class, acc_res, &
					                                info_ieq, sys_x, N_eq)
			endif

			! Extract solution profiles
			if (thereis_mol) then
				! guard tiny negatives from the NL solve
				sys_x(1:N_eq) = max(sys_x(1:N_eq), 0.0d0)
				nhii(j)   = nh(j)*sys_x(1)
				nmol_eq(j,1) = 0.5d0*sys_x(4)*nh(j)
				nmol_eq(j,2) = 0.5d0*sys_x(5)*nh(j)
				nmol_eq(j,3) = sys_x(6)*nh(j)/3.0d0
				nmol_eq(j,4) = sys_x(7)*nh(j)
				if (thereis_oxychem) then
					nox_eq(j,1) = sys_x(iox)  *nm_tot(j,iel_O)
					nox_eq(j,2) = sys_x(iox+1)*nm_tot(j,iel_O)
					n_o1d_eq(j) = excited_oxygen_density(nmol_eq(j,1),    &
					                                     nox_eq(j,2))
				endif
				nhi(j)    = nh(j)*max(1.0d0 - sys_x(1) - sys_x(4)      &
				              - sys_x(5) - sys_x(6) - sys_x(7), 0.0d0)
				! The H nuclei bound in OH and H2O are not atomic H.
				if (thereis_oxychem) nhi(j) = max(nhi(j)                  &
				              - nox_eq(j,1) - 2.0d0*nox_eq(j,2), 0.0d0)
				nheii(j)  = nhe(j)*sys_x(2)
				nheiii(j) = nhe(j)*sys_x(3)
				nhei(j)   = max(nhe(j)*(1.0d0 - sys_x(2) - sys_x(3))   &
				              - nmol_eq(j,4), 0.0d0)
				if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(8)
			else
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)
			nhei(j)   = nhe(j)*(1.0 - sys_x(2) - sys_x(3))
			nheii(j)  = nhe(j)*sys_x(2)
			nheiii(j) = nhe(j)*sys_x(3)
			if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(4)
			endif
			if (thereis_metals) then
				do im = 1,n_melem
					i0 = melem_i0(im)
					nm(j,i0+1) = nm_tot(j,im)*sys_x(mbase+2*(im-1))
					if (melem_top(im) .ge. 2) then
						nm(j,i0+2) = nm_tot(j,im)*sys_x(mbase+1+2*(im-1))
						nm(j,i0)   = nm_tot(j,im)                            &
						           *(1.0 - sys_x(mbase+2*(im-1)) - sys_x(mbase+1+2*(im-1)))
					else
						! Two-stage element: neutral = total - singly ionized.
						nm(j,i0)   = nm_tot(j,im)*(1.0 - sys_x(mbase+2*(im-1)))
					endif
				enddo
				! O I means FREE ATOMIC neutral oxygen with the oxygen
				! chemistry on: the closure above still holds the oxygen
				! bound in OH and H2O, which the residual removed from the
				! same quantity (System_HeH_mol_metals). Removing it here as
				! well is what keeps the written column and the solved one
				! the same object.
				if (thereis_oxychem)                                       &
					nm(j,melem_i0(iel_O)) = max(nm(j,melem_i0(iel_O))      &
					        - nox_eq(j,1) - nox_eq(j,2), 0.0d0)
			endif

			if (last_self) exit self_field_HHe
			enddo self_field_HHe

		enddo
		!$omp end parallel do

		! Carry the columns past this block with the composition just
		! solved: the ground singlet again, and the molecular hydrogen.
		do j = jb_hi, jb_lo, -1
			nheiS_face(j) = nhei(j)
			if (thereis_HeITR) nheiS_face(j) =                             &
			         he_ground_singlet_density(nhei(j), nheiTR(j))
			if (has_h2_xuv) nh2_entry(j) = nmol_eq(j,1)
		enddo
		call advance_starward_columns(jb_hi, jb_lo, nhi, nheiS_face,       &
		         nheii, nheiTR, nh2_entry, nm, has_h2_xuv,                 &
		         N1_col,N15_col,N2_col,NTR_col,NH2_col_xuv,Nm_col_xuv)

		jb_hi = jb_lo - 1
		enddo xuv_block_HHe

		ieq_neq_stored   = N_eq
		ieq_mbase_stored = mbase
		ieq_iox_stored   = iox
		ieq_rates_ready  = .true.

		! The frozen background is now a complete sweep old at worst, so the
		! transport-chemistry operator may run. It is filled wherever the
		! molecular network is on, whether or not its carriers are
		! transported: the local chemical root of the H2 row
		! (carrier_h2_chemical_root) reads it to seed a molecular run from an
		! atomic state, and a molecular run without carrier transport is
		! exactly such a run. The ionization stages add the atomic runs that
		! transport them. Every site that IMPOSES a transported fraction
		! tests its own transport key as well, so readiness imposes nothing.
		if (thereis_mol .or. transported_rows_exist()) bg_ready = .true.

		! One summary line per sweep when a molecular cell's roots were all
		! outside the physical simplex, so the closest one was clamped onto the
		! element budget.
		if (thereis_mol .and. n_mol_clamped .gt. 0) then
			write(*,'(A,I0,A)') ' (ioniz_eq) WARNING: every molecular '//     &
				'equilibrium root left the physical simplex at ',              &
				n_mol_clamped, ' cell(s), clamped onto the element budget'
		endif

		! One summary line per sweep whenever a state left the physical
		! simplex: how many cells had their stored state rejected as a
		! starting point, how many first roots were outside the simplex, and
		! how many ended on the ionization balance because no starting point
		! produced an admissible root. A cell that is merely retried because
		! the solver did not reach its tolerance is routine and stays silent
		! here; it is counted in the run-wide totals reported at the end.
		if (n_ieq_reseed + n_ieq_unphys + n_ieq_fail .gt. 0) then
			write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                             &
				' (ioniz_eq) step ', marching_step,                                &
				': ionization roots - ', n_ieq_reseed,                     &
				' stored state(s) rejected, ', n_ieq_unphys,               &
				' root(s) outside the simplex, no admissible root at ',    &
				n_ieq_fail, ' cell(s)'
		endif

		! This sweep's totals, then into the ledger of the state kind the
		! sweep was tagged with. A steady-solver probe never reaches the
		! ledgers of the states the run holds.
		ledger_sweep%n_sweep       = 1
		ledger_sweep%n_reseed      = n_ieq_reseed
		ledger_sweep%n_retry       = n_ieq_retry
		ledger_sweep%n_unphys      = n_ieq_unphys
		ledger_sweep%n_noroot      = n_ieq_fail
		ledger_sweep%n_mol_clamped = n_mol_clamped
		ledger_sweep%n_mol_info(:) = n_mol_info(:)
		ledger_sweep%acc_n(:)      = n_acc(:)
		ledger_sweep%acc_resmax(:) = acc_resmax(:)
		ledger_sweep%hist_conv(:)  = hist_conv(:)
		ledger_sweep%hist_uncv(:)  = hist_uncv(:)
		ledger_sweep%n_cce_attempt = n_cce_attempt
		ledger_sweep%n_cce_root    = n_cce_root
		ledger_sweep%n_cce_solve   = n_cce_solve
		ledger_sweep%cce_seconds   = cce_seconds
		ledger_sweep%n_nonfinite   = n_res_nonfinite
		! Every cell at which no starting point stayed inside the element
		! simplex: the molecular clamp and the atomic handback are the two
		! ways that ends.
		ledger_sweep%n_offsimplex  = n_mol_clamped + n_ieq_fail
		ledger_sweep%viol_worst    = viol_sweep_worst

		! The non-root streak counts CONSECUTIVE sweeps of a state the run
		! holds; a probe sweep leaves it alone (nonroot_streak_update), so
		! its peak belongs to the two ledgers of held states and stays zero
		! in the probe ledger.
		if (ieq_sweep_state_kind .ne. ieq_state_steady_candidate)          &
			ledger_sweep%streak_peak = maxval(ieq_nonroot_streak)

		call add_to_ioniz_eq_ledger(ledger_sweep)

	endif

	if (present(sweep_ledger)) sweep_ledger = ledger_sweep
	ieq_sweep_ledger_last = ledger_sweep


	!----------------------------------!

	! THE COMPOSITION THIS SWEEP RETURNS WEIGHS THE DENSITY IT WAS GIVEN.
	!
	! The mass fractions of a state are f_i = n_i m_i / rho, so
	! sum_i f_i A_i = 1 is the DEFINITION of f_sp and not a tolerance. The
	! sweep holds the density fixed (T2.1: it is an input the sweep may not
	! write) and solves the chemistry cell by cell, and what it returns
	! weighs that density only to the arithmetic of the cell solves: the
	! element totals come back conserved to about 4e-16 and the mass is a
	! weighted sum of them, so the closure leaves 1 by a few units in the
	! last place at every sweep and in the same direction.
	!
	! MEASURED (item L19, 2026-09-15,
	! docs/lhs1140b_stationary_L19_L20_20260915.md): without this
	! projection the departure RATCHETS, about 1e-14 per outer pass of a
	! stationary solve, reaching 5.1e-13 after nineteen passes and 6.2e-13
	! after twenty-four. Nothing bounded it -- the element operator's own
	! test compares a step's closure with the ENTRY closure of that step
	! and is blind to a monotone drift by construction -- and a state that
	! carries it is a state whose two halves describe two gases, which is
	! what a restart then has to choose between (item L18).
	!
	! THE DENSITY IS THE CONSERVED QUANTITY AND THE COMPOSITION IS THE
	! PARTITION ON IT, so the composition is put on the density and never
	! the other way round: one factor per cell over every species, which
	! leaves every element ratio, every ionization split and every
	! metal-to-hydrogen ratio exactly where the cell solve put them and
	! moves only the normalization. It is the same rule load_IC applies to
	! a restart file, so a state written, read back and re-measured passes
	! through one projection and not two different ones.
	!
	! It runs before the heating and the cooling below, so those are
	! contracted with the composition that is returned and written, and
	! before the calc_rho that normalizes f_sp_io, which then reads one to
	! the last bit. EXHALE_MASS_PROJECTION=0 restores the old behaviour for
	! a control experiment; nothing else reads that variable.
	!
	! THE PROJECTION IS BOUNDED, AND IT REFUSES RATHER THAN HIDES.  A
	! renormalization that is allowed to be any size is a way of making a
	! broken closure invisible: whatever the cell solve did to the mass, the
	! composition would come back weighing the density exactly and nothing
	! would say it had been moved.  So the factor is MEASURED and reported,
	! and it has two thresholds (Codex review, 2026-09-15):
	!
	!   |s - 1| <= 1e-12   the rounding of the cell solves, which is what
	!                      this projection exists to remove; silent.
	!   1e-12 < |s-1| <= 1e-8   reported, loudly, with the cell: the sweep
	!                      is moving mass by more than its own arithmetic.
	!   |s - 1| > 1e-8     NOT PROJECTED.  A closure error of that size is a
	!                      failure of the sweep and is reported as one; the
	!                      composition is handed back as the solve left it,
	!                      so the certification's own mass-closure line
	!                      (item L19) shows the failure instead of a
	!                      projection hiding it.
	!
	! The factor is also refused where it is not a usable number at all: a
	! non-finite or non-positive mixture mass is not a composition and
	! multiplying by it would put the NaN into every species column.
	!
	! WHAT THE PROJECTION MOVES, reported with it: every element total moves
	! by exactly |s - 1| in ABSOLUTE nucleus density, because the factor is
	! common to every species -- which is why no element RATIO moves at all,
	! and why the ratio is the wrong thing to report here.
	if (mass_projection_of_the_sweep()) then
		block
		real*8, dimension(1-Ng:N+Ng) :: m_of_comp, s_close
		real*8  :: s_worst
		integer :: im_p, j_worst
		logical :: s_ok
		if (thereis_mol) then
			if (thereis_oxychem) then
				call calc_rho(nhi,nhii,nhei,nheii,nheiii,m_of_comp,   &
				              nm,nmol_eq,nox_eq)
			else
				call calc_rho(nhi,nhii,nhei,nheii,nheiii,m_of_comp,   &
				              nm,nmol_eq)
			endif
		else
			call calc_rho(nhi,nhii,nhei,nheii,nheiii,m_of_comp,nm)
		endif
		s_ok    = .true.
		s_worst = 0.0d0
		j_worst = 0
		do im_p = 1-Ng, N+Ng
			if (m_of_comp(im_p) .gt. 0.0d0 .and.                      &
			    m_of_comp(im_p) .eq. m_of_comp(im_p) .and.            &
			    n_in_dim(im_p)  .gt. 0.0d0 .and.                      &
			    n_in_dim(im_p)  .eq. n_in_dim(im_p)) then
				s_close(im_p) = n_in_dim(im_p)/m_of_comp(im_p)
			else
				s_close(im_p) = 1.0d0
			endif
			if (s_close(im_p) .ne. s_close(im_p) .or.                 &
			    s_close(im_p) .le. 0.0d0) then
				s_close(im_p) = 1.0d0
				s_ok = .false.
			endif
		enddo
		do im_p = 1, N
			if (abs(s_close(im_p) - 1.0d0) .gt. s_worst) then
				s_worst = abs(s_close(im_p) - 1.0d0)
				j_worst = im_p
			endif
		enddo
		ieq_mass_projection_worst = max(ieq_mass_projection_worst,        &
		                                s_worst)
		if (.not. s_ok) then
			n_ieq_mass_projection_refused =                           &
			   n_ieq_mass_projection_refused + 1
			write(*,'(A)') ' (ioniz_eq) mass projection REFUSED: the'//&
			   ' mixture mass of a cell is not a usable number, so'// &
			   ' the composition is handed back unprojected.'
			s_close = 1.0d0
		else if (s_worst .gt. ieq_mass_projection_refuse) then
			n_ieq_mass_projection_refused =                           &
			   n_ieq_mass_projection_refused + 1
			write(*,'(A,ES10.3,A,I0,A,ES10.3,A)')                     &
			   ' (ioniz_eq) mass projection REFUSED: the composition'//&
			   ' this sweep returns weighs its own density to ',      &
			   s_worst, ' at cell ', j_worst, ', above ',             &
			   ieq_mass_projection_refuse, '.'
			write(*,'(A)') '   That is a failure of the cell solve'// &
			   ' and not the rounding this projection removes, so'//  &
			   ' the composition is handed'
			write(*,'(A)') '   back as the solve left it and the'//   &
			   ' certification''s mass-closure line reports it.'
			s_close = 1.0d0
		else if (s_worst .gt. ieq_mass_projection_report) then
			write(*,'(A,ES10.3,A,I0,A)') ' (ioniz_eq) mass'//         &
			   ' projection: the composition this sweep returns'//    &
			   ' weighs its own density to ', s_worst, ' at cell ',   &
			   j_worst, '; every element total is moved by that'//    &
			   ' much in absolute'
			write(*,'(A)') '   nucleus density, and no element ratio'//&
			   ' by anything (the factor is common to every species).'
		endif
		nhi    = nhi   *s_close
		nhii   = nhii  *s_close
		nhei   = nhei  *s_close
		nheii  = nheii *s_close
		nheiii = nheiii*s_close
		nheiTR = nheiTR*s_close
		do im_p = 1, n_mion
			nm(:,im_p) = nm(:,im_p)*s_close
		enddo
		nmol_eq(:,1) = nmol_eq(:,1)*s_close
		nmol_eq(:,2) = nmol_eq(:,2)*s_close
		nmol_eq(:,3) = nmol_eq(:,3)*s_close
		nmol_eq(:,4) = nmol_eq(:,4)*s_close
		nox_eq(:,1)  = nox_eq(:,1) *s_close
		nox_eq(:,2)  = nox_eq(:,2) *s_close
		nox_eq(:,3)  = nox_eq(:,3) *s_close
		end block
	endif

	! HEATING AND COOLING OF THE STATE THIS SWEEP RETURNS.
	!
	! Everything above the sweep is a RATE or a COEFFICIENT: the attenuated
	! field and its columns, the photoionization rates, the photoelectron
	! partition, the recombination and collisional ionization coefficients,
	! the He recombination corrections, the Lyman-Werner dissociation rate
	! and the FUV band photolysis rates. Those are properties of the
	! radiation field and of T, they are what the cell sweep was solved
	! with, and they keep the documented lag of one sweep.
	!
	! The heating and the cooling are not: each of them is a contraction of
	! those rates with the COMPOSITION, and the sweep has just replaced the
	! composition. Assembling them before the sweep described the previous
	! state -- the written heat/cool columns, the explicit energy update and
	! the steady residual all took a heating that belonged to a composition
	! this routine no longer returns. They are therefore assembled here,
	! from the post-sweep densities and the pre-sweep rates.

	! Electron and gas-particle densities of the post-sweep composition. The
	! heating channels below and the cooling need these, not the entry ones.
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol_eq)
	if (thereis_oxychem) then
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq,      &
		               nox_eq)
	else
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq)
	endif
	! The ghost's heavy-particle plus electron count at the composition this
	! application RETURNED: the thermodynamic leg of the contract. It is what
	! turns the ghost's temperature into the pressure and the internal energy
	! the boundary reads, so the contract is not closed until it has stopped
	! moving with the composition.
	if (ghost_contract_on) np_ghost_return = n_tot(1-Ng:0) + ne(1-Ng:0)

	! The heating of the post-sweep composition, channel by channel, from
	! the ONE assembly in utils_ion_eq. The rates it contracts are the
	! pre-sweep ones described above; the densities are the post-sweep ones
	! extracted just now. The channel array it fills is heat_channel_state,
	! which is what output/Heating_breakdown.txt writes, so that file and
	! the heat column of Hydro_ioniz.txt describe one state.
	call heating_of_composition(T_K,                                      &
	         nhi,nhii,nhei,nheii,nheiii,nheiTR, nm, nmol_eq, nox_eq,      &
	         ne, n_tot,                                                   &
	         h1_HI,h1_HeI,h1_HeII,h1_HeTR,h1_H2,h1_m,                     &
	         A31,q31a,q31b,Q31,                                           &
	         k_lw_diss, p_lw_single, k_co_diss, j_h2o_fuv, j_oh_fuv,      &
	         thereis_mol, thereis_oxychem, heat, heat_channel_state)

	! Two channels are also carried as module arrays of their own, because
	! write_output writes the FUV photolysis heating beside the band ledger
	! and the step controller checkpoints both. The channel array is the
	! definition; these are copies of it.
	heat_chem = heat_channel_state(:,15)
	heat_fuv  = heat_channel_state(:,16)

	! Cooling of the post-sweep composition. eval_cool builds its own
	! electron density from the densities it is given, so this is the same
	! ne as above. The coefficients it returns are those of the state it is
	! evaluated at and are NOT the ones the sweep used, so they go into
	! their own arrays and are discarded.
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,                    &
	               rchiiB_post,rcheiiB_post,rcheiiiB_post, rec_m_post,     &
	               aHI_post,aHeI_post,aHeII_post, aion_m_post,             &
	               cool, nheiTR=nheiTR, a_ion_HeITR=aTR_post,              &
	               nmol=nmol_eq, nox=nox_eq)

	!----------------------------------!

	! Density with atomic numbers (nm adds the metal mass under the
	! eos_metals policy)
   if (thereis_mol) then
      if (thereis_oxychem) then
         call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_rec,nm,nmol_eq,nox_eq)
      else
         call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_rec,nm,nmol_eq)
      endif
   else
      call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_rec,nm)
   endif
   if (present(rho_recon)) rho_recon = n_rec/n0

   ! Abundancies profiles
   ! Normalized by the density the step was GIVEN (T2.1), not by the mass
   ! sum of the composition: the sum is the check above.
   f_sp_io(:,1) = nhi/n_in_dim
   f_sp_io(:,2) = nhii/n_in_dim
   f_sp_io(:,3) = nhei/n_in_dim
   f_sp_io(:,4) = nheii/n_in_dim
   f_sp_io(:,5) = nheiii/n_in_dim
   f_sp_io(:,6) = nheiTR/n_in_dim
   ! molecular abundances
   if (thereis_mol) then
      f_sp_io(:,isp_H2)   = nmol_eq(:,1)/n_in_dim
      f_sp_io(:,isp_H2p)  = nmol_eq(:,2)/n_in_dim
      f_sp_io(:,isp_H3p)  = nmol_eq(:,3)/n_in_dim
      f_sp_io(:,isp_HeHp) = nmol_eq(:,4)/n_in_dim
   else
      f_sp_io(:,isp_H2:isp_HeHp) = 0.0d0
   endif

   ! oxygen-carrier abundances
   if (thereis_oxychem) then
      f_sp_io(:,isp_OH)  = nox_eq(:,1)/n_in_dim
      f_sp_io(:,isp_H2O) = nox_eq(:,2)/n_in_dim
      f_sp_io(:,isp_CO)  = nox_eq(:,3)/n_in_dim
   else
      f_sp_io(:,isp_OH:isp_CO) = 0.0d0
   endif

   ! Metal abundances in canonical order (col mion_fsp(im) of f_sp_io).
   do im = 1,n_mion
      f_sp_io(:,mion_fsp(im)) = nm(:,im)/n_in_dim
   enddo

   ! THE DECLARED DIAGNOSTIC POLICY OF THE CELLS ABOVE THE LOWER GHOSTS:
   ! they are handed back as they came (interior_composition_is_held, default
   ! off). The two lower ghost cells keep the composition the sweep solved,
   ! so a base row measured after this call carries the ghost's refresh and
   ! not the interior's.
   if (interior_composition_is_held())                                    &
      f_sp_io(1:N+Ng,:) = f_sp_interior_entry

   ! ---- HAS THE GHOST STOPPED MOVING UNDER ITS OWN MAP? ----
   !
   ! The composition the two lower ghost cells were handed against the one
   ! this application returned, and the same for their particle plus
   ! electron count. WHAT IS MEASURED is the move the map made,
   ! |S(g) - g| <= the accuracy, which is the stopping statement of a fixed
   ! point: the returned state's own move, |S(S(g)) - S(g)|, is that move
   ! times the map's gain, MEASURED at 4.9e-3 in the direction the seed
   ! spans (READ, docs/lhs1140b_p1_step1b_20260919.md section 6.1), so a
   ! ghost accepted here is three decades better than the accuracy as a
   ! fixed point of its own solve. The composition leg is taken in the
   ! composition's own units and the count leg against its own value, for
   ! the reasons stated at the two accuracies.
   ! Above the accuracy the map is applied again at the ghost just returned,
   ! with the cells above the ghosts restored to the composition they were
   ! handed, so that what this routine returns for them stays one
   ! application of one map.
   if (.not. ghost_fp_reached) then
      ghost_fp_move = ghost_composition_distance(                         &
                           f_sp_sweep_entry(1-Ng:0,:), f_sp_io(1-Ng:0,:))
      ghost_thermal_move = ghost_count_distance(np_ghost_entry,           &
                                                np_ghost_return)
      if (ghost_fp_move    .le. ghost_composition_fixed_point_move .and. &
          ghost_thermal_move .le. ghost_count_fixed_point_move) then
         ghost_fp_reached = .true.
      else if (ghost_fp_pass .ge. ghost_composition_fixed_point_passes)   &
      then
         call ghost_fixed_point_refused(ghost_fp_pass, ghost_fp_move,     &
              ghost_thermal_move,                                         &
              'the ghost composition did not reach a fixed point of'//    &
              ' the sweep that returns it')
      else
         f_sp_sweep_entry(1-Ng:0,:) = f_sp_io(1-Ng:0,:)
         f_sp_io = f_sp_sweep_entry
      endif
   endif
   if (ghost_fp_reached) exit ghost_fixed_point
   enddo ghost_fixed_point
   if (ghost_contract_on) then
      ghost_fixed_point_move_sweep   = ghost_fp_move
      ghost_thermal_move_sweep       = ghost_thermal_move
      ghost_fixed_point_passes_sweep = ghost_fp_pass
   endif

   ! Adimensional heating and cooling rates
   heat_out = heat/q0
   cool_out = cool/q0
      
   ! ---- THE SOLVED GHOST'S COUNTS, WHICH ARE NOT THE RESERVOIR'S ----
   !
   ! The electron and heavy-particle counts of the lower ghost at the
   ! composition this sweep returned. They are a MEASUREMENT of the state and
   ! are handed to base_boundary as such; the PRESCRIBED reservoir
   ! (ntot_bc + dp_bc as input_read resolved them, and the p, T, particles
   ! per unit mass and level radius set_base_reservoir was called with) keeps
   ! its own names and is not written from here.
   !
   ! WHY THEY MUST STAY APART. The sweep used to overwrite ntot_bc and dp_bc
   ! with these two numbers, and the comment that justified it described a
   ! base pressure boundary condition p = (ntot_bc + dp_bc) T0 at a pinned
   ! density -- the closure base_boundary replaced. Under the present
   ! boundary nothing downstream reads them back, and where the two are fed
   ! into one another the base ceases to be the base the state was solved at:
   ! restating the reservoir at the solved count moves the base pressure by
   ! 4.45 per cent and the cell-1 continuity row of a certified hot-Uranus
   ! state from 1.2e-08 to 1.0 (MEASURED,
   ! docs/lhs1140b_stationary_D5a_20260918.md section 6.2).
   !
   ! WHAT THE ELECTRON COUNT IS OF. It is the free electron density of the
   ! composition this sweep returned, counted by calc_ne and by nothing else:
   ! the array ne was filled by that routine from the post-sweep densities
   ! just above, so every charge carrier the equation of state and the sweep
   ! count is counted here with the same charge -- H II, He II, He III, the
   ! molecular ions H2+, H3+ and HeH+, and, under eos_metals, each metal
   ! stage. A second sum written out here would be a second policy: on the
   ! LHS 1140 b hot-Uranus ghosts the three molecular ions carry 1.02 and
   ! 0.78 per cent of the electron density of the two ghost cells (MEASURED,
   ! docs/PLAN_20260919_review.md section 4), which is not a rounding
   ! difference. VALIDITY: any network the sweep solves, molecular or
   ! atomic; in an atomic gas the molecular terms of calc_ne are absent and
   ! the value is the H/He (and metal) count it always was. ne is in cm^-3
   ! and n0 is the density scale, so the count handed on is in units of n0,
   ! as the heavy-particle count beside it is.
   !
   ! The heavy-particle count is taken from calc_ntot's n_tot rather than
   ! summed here, so the particle-count policy keeps its single definition;
   ! with an atomic network every nucleus is its own particle and n_tot is
   ! that count either way.
   ghost_ne_solved = ne(1-Ng)/n0
   call set_base_ghost_counts(n_tot(1-Ng)/n0, ghost_ne_solved)

   ! ---- THE GHOST THE SWEEP RETURNED, cell by cell, for the record ----
   !
   ! A MEASUREMENT. Nothing in the solve, the boundary or the energy update
   ! reads any of it; the ghost record (boundary_state_trace.f90) writes it
   ! where EXHALE_GHOST_RECORD asks for it, and with the key unset it is
   ! stored and never read.
   !
   ! THE CHARGE BUDGET is the charge the composition this sweep WROTE OUT
   ! carries, read back through the mass fractions f_sp_io and the charges
   ! of the species table, against the electron count calc_ne made of the
   ! same state. It closes the round trip every consumer of the state takes,
   ! so it is zero to the rounding of that trip and a nonzero value means the
   ! written composition and the reported electron count are two gases.
   !
   ! THE ELEMENTAL BUDGET is the largest relative departure, over helium and
   ! each trace element, of the ghost's nucleus ratio n_El/n_H from the first
   ! physical cell's. The reservoir row states those ratios, so it is the
   ! distance between the gas below the base and the gas above it.
   do jg = 1-Ng, 0
      q_written = f_sp_io(jg,isp_HII)
      if (thereis_He) q_written = q_written + f_sp_io(jg,isp_HeII)        &
                                + 2.0d0*f_sp_io(jg,isp_HeIII)
      if (thereis_mol) q_written = q_written + f_sp_io(jg,isp_H2p)        &
                                 + f_sp_io(jg,isp_H3p)                    &
                                 + f_sp_io(jg,isp_HeHp)
      if (eos_include_metals .and. thereis_metals) then
         do im = 1,n_mion
            if (mion_stage(im) .gt. 0)                                    &
               q_written = q_written                                      &
                         + dble(mion_stage(im))*f_sp_io(jg,mion_fsp(im))
         enddo
      endif
      q_written = q_written*n_in_dim(jg)
      chg_gap = 0.0d0
      if (ne(jg) .gt. 0.0d0) chg_gap = (q_written - ne(jg))/ne(jg)

      elem_gap = 0.0d0
      if (thereis_He .and. nh(jg) .gt. 0.0d0 .and. nh(1) .gt. 0.0d0) then
         ratio_g = nhe(jg)/nh(jg)
         ratio_1 = nhe(1) /nh(1)
         if (ratio_1 .gt. 0.0d0)                                          &
            elem_gap = max(elem_gap, abs(ratio_g - ratio_1)/ratio_1)
      endif
      if (thereis_metals .and. nh(jg) .gt. 0.0d0 .and. nh(1) .gt. 0.0d0)  &
      then
         do im = 1,n_melem
            ratio_g = nm_tot(jg,im)/nh(jg)
            ratio_1 = nm_tot(1,im) /nh(1)
            if (ratio_1 .gt. 0.0d0)                                       &
               elem_gap = max(elem_gap, abs(ratio_g - ratio_1)/ratio_1)
         enddo
      endif

      xg_h2   = 0.0d0
      xg_h2p  = 0.0d0
      xg_h3p  = 0.0d0
      xg_hehp = 0.0d0
      if (thereis_mol .and. nh(jg) .gt. 0.0d0) then
         xg_h2   = 2.0d0*nmol_eq(jg,1)/nh(jg)
         xg_h2p  = 2.0d0*nmol_eq(jg,2)/nh(jg)
         xg_h3p  = 3.0d0*nmol_eq(jg,3)/nh(jg)
         xg_hehp =       nmol_eq(jg,4)/nh(jg)
      endif
      call set_base_ghost_state_record(jg, n_tot(jg), ne(jg),             &
           xg_h2, xg_h2p, xg_h3p, xg_hehp,                                &
           nhii(jg)/max(nh(jg),1.0d-99),                                  &
           nheii(jg)/max(nhe(jg),1.0d-99),                                &
           nheiii(jg)/max(nhe(jg),1.0d-99),                               &
           chg_gap, elem_gap, ghost_acc_res_cell(jg),                     &
           ghost_closure_res_cell(jg), ghost_closure_pass_cell(jg))
   enddo

   ! The ghost this sweep returned, kept so that a later sweep of the same
   ! run can be seeded with it (EXHALE_GHOST_COMPOSITION_SEED, default off).
   call set_previous_sweep_ghost_rows(f_sp_io(1-Ng:0,:))

   ! And what the ghost's own molecular partition closed to, where the base
   ! handoff states one.
   if (base_h2_composition_imposed())                                   &
      call set_base_ghost_closure(ghost_closure_res_sweep,              &
                                  ghost_closure_pass_sweep)

	! What the ghost composition fixed point reached in this sweep, for the
	! boundary report: the last application's move in the ghost's species
	! densities and in its particle plus electron count, the applications it
	! took, and the accuracy they were held to.
	if (ghost_contract_on)                                               &
		call set_base_ghost_fixed_point(ghost_fixed_point_move_sweep,     &
		     ghost_thermal_move_sweep, ghost_fixed_point_passes_sweep,    &
		     ghost_composition_fixed_point_move,                          &
		     ghost_count_fixed_point_move)

	call element_census_verify(cen_ieq, n_io, f_sp_io, rho_is_fixed=.true.)

	contains


	!----------------------------------!

	real*8 function ghost_composition_distance(f_a, f_b) result(d)
	! HOW FAR TWO GHOST COMPOSITIONS STAND APART, in the units the
	! composition is carried and read in: the arguments are the f_sp layout,
	! the species' number density per unit mass of the cell, and the largest
	! move over both ghost cells and every species is taken. A species is not
	! divided by its own value: that weight belongs to no quantity the
	! boundary builds, and the solve that returns the composition pins it
	! absolutely, so a stage at 1e-10 of the cell has no resolved value of
	! its own to be measured against (see the accuracy above). It also needs
	! no floor: the metal columns of a run without metals sit at a value of
	! order 1e-301 and move by at most that.
	real*8, intent(in) :: f_a(1-Ng:0,n_species), f_b(1-Ng:0,n_species)
	integer :: jc, kc
	d = 0.0d0
	do jc = 1-Ng, 0
		do kc = 1, n_species
			d = max(d, abs(f_b(jc,kc) - f_a(jc,kc)))
		enddo
	enddo
	end function ghost_composition_distance

	!----------------------------------!

	real*8 function ghost_count_distance(np_a, np_b) result(d)
	! The same comparison for the ghost's heavy-particle plus electron count.
	real*8, intent(in) :: np_a(1-Ng:0), np_b(1-Ng:0)
	integer :: jc
	real*8  :: scale_c
	d = 0.0d0
	do jc = 1-Ng, 0
		scale_c = max(abs(np_a(jc)), abs(np_b(jc)))
		if (scale_c .le. 0.0d0) cycle
		d = max(d, abs(np_b(jc) - np_a(jc))/scale_c)
	enddo
	end function ghost_count_distance

	!----------------------------------!

	subroutine ghost_fixed_point_refused(npass, move, thermal_move, what)
	! A BOUNDARY NOTHING STANDS BEHIND IS REFUSED, with the numbers, as the
	! ghost's H2 partition closure is refused: a ghost composition that is
	! not a fixed point of the solve that returned it is an iterate of that
	! solve, and the base face state, the base continuity row and every
	! certificate taken from them would be functions of the composition the
	! run was entered at.
	integer, intent(in) :: npass
	real*8,  intent(in) :: move, thermal_move
	character(len=*), intent(in) :: what
	write(*,'(A)') ' (ioniz_eq) STOP: '//trim(what)
	write(*,'(A,I0,A,I0)') '   applications of the ghost composition'//   &
	     ' map ', npass, '   bound ', ghost_composition_fixed_point_passes
	write(*,'(A,ES12.5,A,ES12.5)') '   last move of the ghost species'//  &
	     ' densities ', move, '   of its particle plus electron count ',  &
	     thermal_move
	write(*,'(A,ES12.5,A,ES12.5)') '   accuracy the boundary states:'//   &
	     ' composition ', ghost_composition_fixed_point_move,             &
	     '   count ', ghost_count_fixed_point_move
	flush(6)
	error stop 'ioniz_eq: the ghost composition is not a fixed point'
	end subroutine ghost_fixed_point_refused

	! End of subroutine
	end subroutine ioniz_eq

	!----------------------------------!

	logical function ionization_fractions_physical(x,n,mbase) result(ok)
	! Is a root of the equilibrium system a physically admissible state?
	!
	! Every unknown of these systems is the fraction of one element's nuclei
	! found in one ionization stage, so a physical state has
	!   (i)  every fraction >= 0, and
	!   (ii) for each element, the tracked ionized stages summing to <= 1,
	!        the neutral stage being the remainder 1 - sum.
	! The residuals are polynomials in the fractions and do admit roots
	! outside this simplex; such a root gives a negative neutral density and
	! is not a solution of the physical problem, however small the residual.
	! Downstream it poisons the photoionization integrals (a negative
	! absorber column gives negative photoheating) and the element budget.
	!
	! The comparisons carry a small tolerance so that a root sitting on a
	! face of the simplex (a fully neutral or fully ionized element) is not
	! rejected for round-off; it is far below the violations this test is
	! meant to catch (fractions of ~1e-6 to ~1e-1 outside the simplex).
	!
	! Widening it to the solver's own xtol = sqrt(machine epsilon) was tried
	! and measured: it lets a cell stop on a root that sits just outside a
	! face instead of retrying from another starting point, and the retried
	! root is the better one. On examples/15 (HD 209458 b, molecular, cold)
	! the wider band turns info = 0 at ||R|| = 6.105e-6 in 89 outer
	! iterations into info = 2 at 1.112e-4 in 133, so the band stays at
	! round-off.
	!
	! Layout: x(1) = H II / H, x(2) = He II / He, x(3) = He III / He,
	! x(4) = He 2^3S / He when the triplet is tracked (x(8) in the molecular
	! layout, where x(4..7) are the H nuclei bound in H2, H2+, H3+, HeH+),
	! and the metal stages X+ / X++ from x(mbase) upwards, two per element.

	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	real*8, parameter   :: ftol = 1.0d-10
	integer :: im, ix, iy
	real*8  :: s

	ok = .false.

	! Hydrogen nuclei: H II, plus the H bound in molecules where tracked,
	! plus the H bound in the oxygen carriers (OH one nucleus, H2O two).
	! Those two are fractions of the OXYGEN family, so they are converted
	! with this cell's n_O(free)/n_H, the same way HeH+ is converted into
	! the helium budget below. Leaving them out lets a root whose free
	! atomic H is negative pass the test, and the extraction then clips it
	! to zero and loses the nuclei.
	if (x(1) .lt. -ftol) return
	s = x(1)
	if (thereis_mol) then
		do ix = 4,7
			if (x(ix) .lt. -ftol) return
			s = s + x(ix)
		enddo
	endif
	if (thereis_oxychem .and. ieq_cell%nh .gt. 0.0d0) then
		iy = oxygen_row_base()
		s = s + (x(iy) + 2.0d0*x(iy+1))*ieq_cell%n_ofam/ieq_cell%nh
	endif
	if (s .gt. 1.0d0 + ftol) return

	! Helium nuclei: He II, He III, the 2^3S metastable (which the systems
	! carry as a separate level inside the neutral stage) and, in the
	! molecular layout, the He nucleus bound in HeH+. The molecular systems
	! close the neutral He as n_He(1 - x2 - x3) - n_HeH+, so a root leaving
	! that negative is inadmissible even when x2 + x3 <= 1; x(7) is HeH+ per
	! H nucleus, hence the n_H/n_He conversion from this cell's state.
	if (thereis_He) then
		if (x(2) .lt. -ftol) return
		if (x(3) .lt. -ftol) return
		s = x(2) + x(3)
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			if (x(ix) .lt. -ftol) return
			s = s + x(ix)
		endif
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0)                &
			s = s + x(7)*ieq_cell%nh/ieq_cell%nhe
		if (s .gt. 1.0d0 + ftol) return
	endif

	! Each metal element separately: X+ (+ X++ for the three-stage elements).
	! With the oxygen chemistry on the oxygen element carries two more
	! fractions, n_OH and n_H2O over the free oxygen family, and the four of
	! them share one budget: free atomic O is the remainder.
	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			if (x(ix) .lt. -ftol) return
			s = x(ix)
			if (melem_top(im) .ge. 2) then
				if (x(ix+1) .lt. -ftol) return
				s = s + x(ix+1)
			endif
			if (thereis_oxychem .and. im .eq. iel_O) then
				iy = oxygen_row_base()
				if (x(iy)   .lt. -ftol) return
				if (x(iy+1) .lt. -ftol) return
				s = s + x(iy) + x(iy+1)
			endif
			if (s .gt. 1.0d0 + ftol) return
		enddo
	endif

	ok = .true.

	end function ionization_fractions_physical

	!----------------------------------!

	subroutine ionization_balance_at_fixed_ne(x,n,mbase,n_e)
	! Stage fractions of each element in its OWN ionization balance at a
	! given electron density: photoionization plus electron-impact ionization
	! against radiative recombination,
	!
	!    n_k (gamma_k + beta_k n_e) = alpha_{k+1} n_e n_{k+1},
	!
	! with the couplings between elements (charge exchange, and the
	! dependence of n_e on the unknowns themselves) dropped and n_e taken
	! from the incoming state. Writing u_k = gamma_k + beta_k n_e for the
	! rate out of stage k and d_k = alpha_k n_e for the rate back into stage
	! k-1, the three-stage solution is
	!
	!    (n_0, n_1, n_2) proportional to (d_1 d_2, u_0 d_2, u_0 u_1),
	!
	! non-negative and normalized to one, so the result always lies inside
	! the physical simplex whatever state it is asked to replace. It is used
	! as a starting point for the coupled solve, and as the state of last
	! resort for a cell where no starting point produced an admissible root.
	!
	! The He 2^3S fraction is the steady state of the metastable level at
	! those populations: recombination and collisional excitation feed it,
	! radiative decay, collisional de-excitation, photoionization and Penning
	! ionization on H0 drain it. It is capped by the neutral He fraction.
	!
	! Rates come from the same cell state the residuals read (ieq_cell and
	! the met_* coefficients of System_HeH_metals), so this is the incoming
	! cell's own physics, not a generic guess. Molecules are not part of this
	! balance: it is the molecule-free limit of the state, so in the molecular
	! layout the H2/H2+/H3+/HeH+ fractions x(4..7) are left at zero and the
	! metastable sits at x(8).

	integer, intent(in)  :: n, mbase
	real*8,  intent(in)  :: n_e
	real*8,  intent(out) :: x(n)
	integer :: itr
	real*8  :: u0,u1,d1,d2,w1,w2,s
	real*8  :: xneu, n_hi_loc

	x(1:n) = 0.0d0

	! Hydrogen
	u0 = ieq_cell%P_HI + ieq_cell%a_ion_HI*n_e
	d1 = ieq_cell%rchiiB*n_e
	s  = u0 + d1
	if (s .gt. 0.0d0) x(1) = u0/s

	if (thereis_He) then
		! Helium. Both He+ recombination channels return to neutral He: the
		! metastable capture rcheiTR is a branch of the recombination, so it
		! adds to the rate back into He I when the triplet is tracked.
		u0 = ieq_cell%P_HeI + ieq_cell%a_ion_HeI*n_e
		d1 = ieq_cell%rcheiiB*n_e
		if (thereis_HeITR) d1 = d1 + ieq_cell%rcheiTR*n_e
		u1 = ieq_cell%P_HeII + ieq_cell%a_ion_HeII*n_e
		d2 = ieq_cell%rcheiiiB*n_e
		w1 = u0*d2
		w2 = u0*u1
		s  = d1*d2 + w1 + w2
		if (s .gt. 0.0d0) then
			x(2) = w1/s
			x(3) = w2/s
		endif

		if (thereis_HeITR) then
			itr = 4
			if (thereis_mol) itr = 8
			xneu     = max(1.0d0 - x(2) - x(3), 0.0d0)
			n_hi_loc = max(1.0d0 - x(1), 0.0d0)*ieq_cell%nh
			s = ieq_cell%P_HeITR + ieq_cell%A31                          &
			  + n_hi_loc*ieq_cell%Q31                                    &
			  + (ieq_cell%q31a + ieq_cell%q31b                           &
			     + ieq_cell%a_ion_HeITR)*n_e
			if (s .gt. 0.0d0) x(itr) = min(                              &
			      n_e*(x(2)*ieq_cell%rcheiTR + xneu*ieq_cell%q13)/s, xneu)
		endif
	endif

	if (thereis_metals)                                                  &
		call metal_ionization_balance_at_fixed_ne(x,n,mbase,n_e)

	end subroutine ionization_balance_at_fixed_ne

	!----------------------------------!

	real*8 function element_budget_violation(x,n,mbase) result(viol)
	! How far a root of an equilibrium system lies outside the physically
	! allowed states: the largest of the negative stage fractions and of the
	! amounts by which one element's tracked stages exceed its nuclei, both
	! measured as fractions of the element. Zero for an admissible state.
	! Used to pick, among roots that are all inadmissible, the one closest to
	! a state (see clamp_fractions_to_element_budget). Same layout as
	! ionization_fractions_physical.

	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	integer :: im, ix, iy
	real*8  :: s

	viol = 0.0d0
	do ix = 1,n
		viol = max(viol, -x(ix))
	enddo

	s = x(1)
	if (thereis_mol) s = s + x(4) + x(5) + x(6) + x(7)
	! The H nuclei bound in OH and H2O (per H nucleus), as in
	! ionization_fractions_physical.
	if (thereis_oxychem .and. ieq_cell%nh .gt. 0.0d0) then
		iy = oxygen_row_base()
		s = s + (x(iy) + 2.0d0*x(iy+1))*ieq_cell%n_ofam/ieq_cell%nh
	endif
	viol = max(viol, s - 1.0d0)

	if (thereis_He) then
		s = x(2) + x(3)
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			s = s + x(ix)
		endif
		! The He nucleus bound in HeH+ (x(7) is per H nucleus), as in
		! ionization_fractions_physical.
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0)                &
			s = s + x(7)*ieq_cell%nh/ieq_cell%nhe
		viol = max(viol, s - 1.0d0)
	endif

	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			s  = x(ix)
			if (melem_top(im) .ge. 2) s = s + x(ix+1)
			! The oxygen budget also carries OH and H2O.
			if (thereis_oxychem .and. im .eq. iel_O) then
				iy = oxygen_row_base()
				s  = s + x(iy) + x(iy+1)
			endif
			viol = max(viol, s - 1.0d0)
		enddo
	endif

	end function element_budget_violation

	!----------------------------------!

	double precision function normalized_reaction_residual(x,n,mbase,       &
	                                              n_e_ref,row_out) result(res)
	! Dimensionless reaction imbalance of a candidate equilibrium state: the
	! residual vector of the SAME system the cell solve used, evaluated at x,
	! with each balance row measured against its turnover scale -- the rate
	! at which the species that row balances can be produced or destroyed in
	! this cell (the section-93 convention: every species set to the whole of
	! its element). A root has res -> 0 to the accuracy of the solve; element
	! bounds and conservation alone say nothing about it, which is why the
	! acceptance test of ioniz_eq requires BOTH
	! (docs/supersonic_molecular_base.md section 11.5).
	!
	! The molecular systems already reach hybrd1 with their rows divided by
	! the turnover scale (set_mol_turnover_rates, set_mol_metal_turnover_-
	! rates), so their residual is used as returned, widened only by the
	! metal <-> H/He charge-exchange bound that the solver scale leaves out
	! of its rows: one contribution among several for the solver's path, but
	! in the wind, where charge exchange couples a trace metal's rows to the
	! whole H reservoir, the dominant term of the row's turnover -- judging a
	! candidate against a scale without it would reject roots. The atomic
	! systems are unscaled inside the solver (scaling would change the
	! MINPACK path and the byte-identical history for no physical need), so
	! the same style of scale is built here from this cell's own
	! coefficients (ieq_cell, met_*): photo- and collisional ionization,
	! recombination, charge exchange, and the He 2^3S kinetics where the
	! triplet is tracked. n_e_ref is the cell's incoming electron density,
	! as in set_mol_turnover_rates. A row whose scale vanishes carries no
	! reaction: its residual is exactly zero and the scale is left at 1, as
	! are the dimensionless identity rows (absent elements, the pinned X++
	! of a two-stage element, a transported carrier's x - x_fix row).
	!
	! Serves the He-branch acceptance only; the H-only system (System_H) has
	! a single analytic balance with a unique root in [0,1] and no
	! acceptance ladder.
	!
	! row_out (optional) returns the same measure ROW BY ROW,
	! |fvec(i)|/srow(i), so that a caller can name the balance that is out
	! rather than only its largest value: the He 2^3S level row of B1a
	! section 2.5 is one row of this vector, and the certification reports
	! it beside the closure it belongs to.
	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	real*8,  intent(in) :: n_e_ref
	real*8,  intent(out), optional :: row_out(n)
	real*8  :: fv(n), srow(n), cxb(n), par(60), el_tot(12)
	real*8  :: nH, nHe, ne, cx_heh, s
	integer :: iflag, i, e, ix

	par(:)    = 0.0d0                 ! transport argument only, unread
	iflag     = 1
	srow(1:n) = 1.0d0

	if (thereis_metals) then
		el_tot(:) = 0.0d0
		el_tot(1:met_nelem) = met_ntot(1:met_nelem)
		el_tot(11) = ieq_cell%nh      ! cx_H
		el_tot(12) = ieq_cell%nhe     ! cx_He
	endif

	if (.not. thereis_He) then
		! The H-only system: one balance row, and its turnover is the rate
		! at which the cell can make or unmake a proton -- photoionization
		! and collisional ionization of H0, radiative recombination.
		call ion_system_H(n,x,fv,iflag,par)
		srow(1) = (ieq_cell%P_HI                                          &
		           + (ieq_cell%a_ion_HI + ieq_cell%rchiiB)*n_e_ref)       &
		          *ieq_cell%nh
	else if (thereis_mol) then
		if (thereis_metals) then
			cx_metal_base = mbase
			call ion_system_HeH_mol_metals(n,x,fv,iflag,par)
			! Charge-exchange bound of each row (absolute), converted to
			! the dimensionless factor (s_solver + s_cx)/s_solver each
			! already-scaled row is divided by. Rows charge exchange never
			! reaches keep exactly 1.
			cxb(1:n) = 0.0d0
			call cx_add_to_turnover(cxb, el_tot)
			cx_metal_base = 4
			do i = 1,n
				srow(i) = 1.0d0 + cxb(i)*mol_inv_turnover(i)
			enddo
		else
			call ion_system_HeH_mol(n,x,fv,iflag,par)
		endif
	else
		! Atomic layouts: raw balance rows [cm^-3 s^-1] of the system the
		! solve used, then the turnover scale of each row.
		if (thereis_HeITR .and. thereis_metals) then
			cx_metal_base = 5
			call ion_system_HeH_TR_metals(n,x,fv,iflag,par)
			cx_metal_base = 4
		else if (thereis_HeITR) then
			call ion_system_HeH_TR(n,x,fv,iflag,par)
		else if (thereis_metals) then
			call ion_system_HeH_metals(n,x,fv,iflag,par)
		else
			call ion_system_HeH(n,x,fv,iflag,par)
		endif

		nH  = ieq_cell%nh
		nHe = ieq_cell%nhe
		ne  = n_e_ref
		cx_heh = (ieq_cell%kcx_He0_Hp + ieq_cell%kcx_Hep_H0)*nH*nHe
		! (1) H+ : photo- and collisional ionization of H0, radiative
		!     recombination, He <-> H charge exchange; with the triplet
		!     tracked, the Penning ionization source of heh_tr_rows.
		srow(1) = (ieq_cell%P_HI                                          &
		           + (ieq_cell%a_ion_HI + ieq_cell%rchiiB)*ne)*nH + cx_heh
		if (thereis_HeITR) then
			srow(1) = srow(1) + f_penning_HeI23S*ieq_cell%Q31*nHe*nH
			! (2) summed He I balance of heh_tr_rows: both photoionization
			!     channels, collisional ionization of ground and metastable
			!     He I, and both recombination paths back into He I.
			srow(2) = (ieq_cell%P_HeI + ieq_cell%P_HeITR                  &
			           + (ieq_cell%a_ion_HeI + ieq_cell%a_ion_HeITR       &
			              + ieq_cell%rcheiTR + ieq_cell%rcheiiB)*ne)*nHe  &
			        + cx_heh
			! (3) He+ <-> He++.
			srow(3) = (ieq_cell%P_HeII                                    &
			           + (ieq_cell%a_ion_HeII + ieq_cell%rcheiiiB)*ne)*nHe
			! (4) He 2^3S: populated from He+ recombination and 1^1S
			!     excitation, drained by photoionization, A31,
			!     de-excitation, electron-impact and Penning ionization.
			srow(4) = (ieq_cell%P_HeITR + ieq_cell%A31                    &
			           + (ieq_cell%rcheiTR + ieq_cell%q13 + ieq_cell%q31a &
			              + ieq_cell%q31b + ieq_cell%a_ion_HeITR)*ne)*nHe &
			        + ieq_cell%Q31*nHe*nH
		else
			! (2)(3) the standard heh_rows balances.
			srow(2) = (ieq_cell%P_HeI                                     &
			           + (ieq_cell%a_ion_HeI + ieq_cell%rcheiiB)*ne)*nHe  &
			        + cx_heh
			srow(3) = (ieq_cell%P_HeII                                    &
			           + (ieq_cell%a_ion_HeII + ieq_cell%rcheiiiB)*ne)*nHe
		endif
		if (thereis_metals) then
			! X0 <-> X+ and X+ <-> X++ of each element (the metal_rows
			! balances), plus the charge-exchange bound on every row it
			! reaches; identity rows are then forced back to their
			! dimensionless scale of 1.
			do e = 1,met_nelem
				ix = mbase + 2*(e-1)
				if (met_ntot(e) .le. 1.0d-30) cycle
				srow(ix) = met_ntot(e)*(met_g0(e)                         &
				           + (met_b0(e) + met_a1(e))*ne)
				if (met_top(e) .ge. 2)                                    &
					srow(ix+1) = met_ntot(e)*(met_g1(e)                   &
					             + (met_b1(e) + met_a2(e))*ne)
			enddo
			cx_metal_base = mbase
			call cx_add_to_turnover(srow, el_tot)
			cx_metal_base = 4
			do e = 1,met_nelem
				ix = mbase + 2*(e-1)
				if (met_ntot(e) .le. 1.0d-30) then
					srow(ix)   = 1.0d0
					srow(ix+1) = 1.0d0
				else if (met_top(e) .lt. 2) then
					srow(ix+1) = 1.0d0
				endif
			enddo
		endif
	endif

	res = 0.0d0
	do i = 1,n
		s = srow(i)
		if (s .le. 0.0d0) s = 1.0d0
		res = max(res, abs(fv(i))/s)
		if (present(row_out)) row_out(i) = abs(fv(i))/s
	enddo

	end function normalized_reaction_residual

	!----------------------------------!

	subroutine keep_background_of_adopted_state
	! bg_cell_adopted <- bg_cell. Called by the steady solver right after it
	! has evaluated the residual at a state it HOLDS (its current iterate, or
	! a trial the line search has just accepted), so that the background of
	! that state survives the probe sweeps that follow.
	if (.not. transported_rows_exist()) return
	if (.not. allocated(bg_cell_adopted)) allocate(bg_cell_adopted(1-Ng:N+Ng))
	bg_cell_adopted = bg_cell
	end subroutine keep_background_of_adopted_state

	!----------------------------------!

	subroutine keep_background_of_best_iterate
	! bg_cell_best <- bg_cell_adopted. Called where the steady solver records
	! a new best iterate, alongside Ybest and f_sp_best.
	if (.not. transported_rows_exist()) return
	if (.not. allocated(bg_cell_adopted)) return
	if (.not. allocated(bg_cell_best)) allocate(bg_cell_best(1-Ng:N+Ng))
	bg_cell_best = bg_cell_adopted
	end subroutine keep_background_of_best_iterate

	!----------------------------------!

	subroutine adopt_background_of_best_iterate
	! bg_cell_adopted <- bg_cell_best. Called where the steady solver returns
	! the best iterate rather than the last one visited.
	if (.not. transported_rows_exist()) return
	if (.not. allocated(bg_cell_best)) return
	if (.not. allocated(bg_cell_adopted)) allocate(bg_cell_adopted(1-Ng:N+Ng))
	bg_cell_adopted = bg_cell_best
	end subroutine adopt_background_of_best_iterate

	!----------------------------------!

	subroutine install_background_of_adopted_state
	! bg_cell <- bg_cell_adopted. The last thing the steady solver does, so
	! that the background the carrier transport reads describes the state the
	! solver hands back and not the last state it happened to evaluate.
	if (.not. transported_rows_exist()) return
	if (.not. allocated(bg_cell_adopted)) return
	bg_cell = bg_cell_adopted
	end subroutine install_background_of_adopted_state

	!----------------------------------!

	logical function transported_rows_exist() result(any_row)
	! Does the transport-chemistry operator have a row to solve in this run?
	!
	! Two independent sets of rows reach that operator. The molecular
	! carriers -- H2, and OH, H2O, CO with the oxygen cycle -- exist where
	! the molecular network is on and its transport is selected. The three
	! ionization stages x(H II), x(He II) and x(He III) are rows on their
	! own key and in ANY gas: they are stages of an element, transported on
	! that element's own nucleus face flux, and a hydrogen and helium
	! mixture has those stages whether or not it has molecules.
	!
	! Everything that follows from "the operator runs" is conditioned on
	! this and not on the molecular configuration alone: its two entry
	! points, the frozen background the rows are evaluated on, the pass cap
	! and alternation of the stationary outer iteration, and the rows the
	! certification measures.
	any_row = (thereis_mol .and. carrier_transport) .or.                  &
	          ionization_transport
	end function transported_rows_exist

	!----------------------------------!

	logical function finite_real(x)
	! Whether x is an ordinary real number: not a NaN, not an infinity.
	! NaN fails every comparison including with itself, and an infinity
	! exceeds the largest representable finite value, so the two tests
	! together cover both. Written out rather than taken from the
	! ieee_arithmetic intrinsic module so that the build's generated module
	! dependency graph stays over the source tree; the code is compiled
	! without any fast-math option, so the IEEE comparison semantics this
	! relies on hold.
	real*8, intent(in) :: x
	finite_real = (x .eq. x) .and. (abs(x) .le. huge(x))
	end function finite_real

	!----------------------------------!

	logical function mass_projection_of_the_sweep() result(on)
	! Whether the composition a sweep returns is put on the density the
	! sweep was given (item L19). Default ON, because the mass fractions
	! are defined by that density; EXHALE_MASS_PROJECTION=0 is the control
	! experiment that restores the unprojected composition.
	!
	! Read once and kept: ioniz_eq is the hot path of every route in the
	! code and is entered thousands of times per solve, so the environment
	! is not re-read per sweep.
	character(len=8) :: env
	if (.not. mass_projection_known) then
		env = ' '
		call get_environment_variable('EXHALE_MASS_PROJECTION', env)
		mass_projection_on    = (trim(env) .ne. '0')
		mass_projection_known = .true.
	endif
	on = mass_projection_on
	end function mass_projection_of_the_sweep

	!----------------------------------!

	subroutine set_ioniz_eq_sweep_state_kind(kind)
	! Tag the state the following ioniz_eq sweeps are evaluating:
	! ieq_state_marching (the state the run holds and advances),
	! ieq_state_steady_iterate (the steady solver's current iterate, also a
	! state the run holds), or ieq_state_steady_candidate (a
	! finite-difference or Krylov probe, or a line-search trial the solver
	! has not adopted). It selects the ledger the sweep totals are added to
	! and decides whether the non-root streak and the acceptance report see
	! the sweep at all. Called from OUTSIDE any parallel region.
	integer, intent(in) :: kind
	ieq_sweep_state_kind = kind
	end subroutine set_ioniz_eq_sweep_state_kind

	!----------------------------------!

	subroutine add_to_ioniz_eq_ledger(a)
	! Add one sweep's totals to the ledger of the state kind currently
	! tagged. Serial: called once per sweep, outside the cell loop.
	type(ioniz_eq_ledger), intent(in) :: a

	select case (ieq_sweep_state_kind)
	case (ieq_state_steady_iterate)
		call accumulate_ioniz_eq_ledger(ieq_steady_iterate_ledger, a)
	case (ieq_state_steady_candidate)
		call accumulate_ioniz_eq_ledger(ieq_steady_candidate_ledger, a)
	case default
		call accumulate_ioniz_eq_ledger(ieq_marching_ledger(ledger_family), a)
	end select
	end subroutine add_to_ioniz_eq_ledger

	!----------------------------------!

	subroutine accumulate_ioniz_eq_ledger(tot, a)
	! tot <- tot + a, counters summed and the "largest seen" fields maxed.
	type(ioniz_eq_ledger), intent(inout) :: tot
	type(ioniz_eq_ledger), intent(in)    :: a

	tot%n_sweep       = tot%n_sweep       + a%n_sweep
	tot%n_reseed      = tot%n_reseed      + a%n_reseed
	tot%n_retry       = tot%n_retry       + a%n_retry
	tot%n_unphys      = tot%n_unphys      + a%n_unphys
	tot%n_noroot      = tot%n_noroot      + a%n_noroot
	tot%n_mol_clamped = tot%n_mol_clamped + a%n_mol_clamped
	tot%n_mol_info(:) = tot%n_mol_info(:) + a%n_mol_info(:)
	tot%acc_n(:)      = tot%acc_n(:)      + a%acc_n(:)
	tot%acc_resmax(:) = max(tot%acc_resmax(:), a%acc_resmax(:))
	tot%hist_conv(:)  = tot%hist_conv(:)  + a%hist_conv(:)
	tot%hist_uncv(:)  = tot%hist_uncv(:)  + a%hist_uncv(:)
	tot%n_cce_attempt = tot%n_cce_attempt + a%n_cce_attempt
	tot%n_cce_root    = tot%n_cce_root    + a%n_cce_root
	tot%n_cce_solve   = tot%n_cce_solve   + a%n_cce_solve
	tot%cce_seconds   = tot%cce_seconds   + a%cce_seconds
	tot%n_nonfinite   = tot%n_nonfinite   + a%n_nonfinite
	tot%n_offsimplex  = tot%n_offsimplex  + a%n_offsimplex
	tot%viol_worst    = max(tot%viol_worst, a%viol_worst)
	tot%streak_peak   = max(tot%streak_peak, a%streak_peak)
	end subroutine accumulate_ioniz_eq_ledger

	!----------------------------------!

	subroutine write_ioniz_eq_acceptance_report
	! End-of-run acceptance report, one block per state kind: the states the
	! marching loop held, the states the steady solver held as its iterate,
	! and the states the steady solver only probed. Keeping them apart is
	! what makes the first two readable at all -- a JFNK phase evaluates its
	! residual tens of times per outer iteration on states no one adopted,
	! and their acceptance classes swamp the run's own (section 121).

	call write_one_ioniz_eq_ledger(ieq_marching_ledger(ledger_family_init), &
		'marching states')
	call write_one_ioniz_eq_ledger(ieq_marching_ledger(ledger_family_phys), &
		'marching states (accepted physical steps)')
	call write_one_ioniz_eq_ledger(ieq_steady_iterate_ledger,             &
		'steady-solver iterates')
	call write_one_ioniz_eq_ledger(ieq_steady_candidate_ledger,           &
		'steady-solver probes (states the run did not adopt)')
	end subroutine write_ioniz_eq_acceptance_report

	!----------------------------------!

	subroutine write_one_ioniz_eq_ledger(led, label)
	! The acceptance block of one ledger. Silent for a ledger nothing was
	! recorded in, so an ordinary marching run prints exactly the lines it
	! printed before the split.
	type(ioniz_eq_ledger), intent(in) :: led
	character(len=*),      intent(in) :: label

	if (led%n_sweep .le. 0) return

	write(*,'(A,A,A,I0,A)') '     ioniz-eq [', label, ']: ',              &
		led%n_sweep, ' equilibrium sweep(s)'

	if (led%n_retry + led%n_unphys + led%n_reseed + led%n_noroot .gt. 0)  &
		write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                                &
			'     ioniz-eq roots: ', led%n_retry,                         &
			' cell solve(s) restarted, ', led%n_unphys,                   &
			' root(s) outside the physical simplex, ', led%n_reseed,      &
			' stored state(s) rejected, ', led%n_noroot,                  &
			' cell(s) left on the ionization balance'

	if (led%n_mol_clamped .gt. 0)                                         &
		write(*,'(A,I0,A)')                                               &
			'     ioniz-eq molecular: ', led%n_mol_clamped,               &
			' cell(s) clamped onto the element budget'

	if (sum(led%n_mol_info) .gt. 0)                                       &
		write(*,'(A,6(I0,A))')                                            &
			'     ioniz-eq molecular hybrd1 info: 1 ', led%n_mol_info(1), &
			', 2 ', led%n_mol_info(2), ', 3 ', led%n_mol_info(3),         &
			', 4 ', led%n_mol_info(4), ', 5 ', led%n_mol_info(5),         &
			', 0 ', led%n_mol_info(0), ''

	if (sum(led%acc_n) .gt. 0) then
		write(*,'(A,I0,A,ES9.2,A,I0,A,ES9.2,A,I0,A,ES9.2,A)')             &
			'     ioniz-eq acceptance: ', led%acc_n(1),                   &
			' converged root(s) (max res ', led%acc_resmax(1), '), ',     &
			led%acc_n(2), ' root(s) without solver convergence (max res ',&
			led%acc_resmax(2), '), ', led%acc_n(3),                       &
			' projected/handback root(s) (max res ', led%acc_resmax(3),   &
			')'
		if (led%acc_n(4) .gt. 0)                                          &
			write(*,'(A,I0,A,ES9.2,A,I0,A)')                              &
			'     ioniz-eq acceptance WARNING: ', led%acc_n(4),           &
			' NON-ROOT state(s) accepted under the relaxation amnesty '// &
			'(max res ', led%acc_resmax(4), ', longest cell streak ',     &
			led%streak_peak, ' sweep(s))'
		if (led%acc_n(6) .gt. 0)                                          &
			write(*,'(A,I0,A,ES9.2,A,ES9.2,A)')                           &
			'     ioniz-eq acceptance WARNING: ', led%acc_n(6),           &
			' cell state(s) were NON-ROOTS above the amnesty cap (max'//  &
			' res ', led%acc_resmax(6), ', cap ', ieq_nonroot_res_cap,    &
			'); each kept the composition it entered the sweep with'
		if (led%acc_n(5) .gt. 0)                                          &
			write(*,'(A,I0,A,ES9.2,A)')                                   &
			'     ioniz-eq acceptance: ', led%acc_n(5),                   &
			' constrained-continuation root(s) (max res ',                &
			led%acc_resmax(5), ')'
		if (led%n_cce_attempt .gt. 0)                                     &
			write(*,'(A,I0,A,I0,A,I0,A,F0.3,A)')                          &
			'     ioniz-eq constrained solve: ', led%n_cce_attempt,       &
			' cell(s) promoted, ', led%n_cce_root,                        &
			' accepted as root(s), ', led%n_cce_solve,                    &
			' field solve(s), ', led%cce_seconds, ' s'
		write(*,'(A,16(I0,1X))')                                          &
			'     ioniz-eq residual decades (<=1e-16 .. >=1e-1) '//       &
			'converged:   ', led%hist_conv
		if (sum(led%hist_uncv) .gt. 0)                                    &
			write(*,'(A,16(I0,1X))')                                      &
			'     ioniz-eq residual decades (<=1e-16 .. >=1e-1) '//       &
			'unconverged: ', led%hist_uncv
	endif

	if (led%n_nonfinite .gt. 0)                                           &
		write(*,'(A,I0,A)')                                               &
			'     ioniz-eq WARNING: ', led%n_nonfinite,                   &
			' cell state(s) with a non-finite reaction residual'

	if (led%n_offsimplex .gt. 0 .or. led%viol_worst .gt. 0.0d0)           &
		write(*,'(A,I0,A,ES9.2)')                                         &
			'     ioniz-eq admissibility: ', led%n_offsimplex,            &
			' cell(s) with no in-simplex starting point, largest '//      &
			'element-budget violation ', led%viol_worst

	end subroutine write_one_ioniz_eq_ledger

	!----------------------------------!

	integer function residual_decade(res) result(ib)
	! Decade bin of a normalized reaction residual for the run-wide
	! histograms: bin 0 collects everything at or below 1e-16 (exact zeros
	! included), bin 15 everything at or above 1e-1, and bin k in between
	! covers [10^(k-16), 10^(k-15)).
	!
	! The top bin is selected by the bracket itself rather than by clamping
	! log10, so that a residual which is NOT a finite number lands there too.
	! Such a residual does reach this routine: the steady (JFNK) solver
	! evaluates its own residual on probe states whose stage fractions can
	! drive a balance row to NaN or to overflow, and the equilibrium sweep
	! then measures a non-finite reaction residual for that probe. NaN
	! fails every comparison, so it falls through to the last branch. The
	! former form evaluated log10 of it and converted the result with int(),
	! which yields the integer indefinite value -2^31 and turned the bin
	! into a wild array index into hist_conv / hist_uncv.
	real*8, intent(in) :: res
	if (res .le. 1.0d-16) then
		ib = 0
	else if (res .lt. 1.0d-1) then
		ib = 16 + int(floor(log10(res)))
	else
		ib = 15
	endif
	end function residual_decade

	!----------------------------------!

	subroutine report_accepted_cell_state(j, acc_class, acc_res, info_solver,&
	                                      x, nx)
	! The accepted equilibrium state of ONE cell, written by every sweep that
	! solves it. The normalized reaction residual is the quantity the
	! acceptance is taken on (ieq_res_tol), so it says how wide the band of
	! compositions this cell would have accepted is; the hydrogen partition
	! says where in that band this sweep landed.
	integer, intent(in) :: j, acc_class, info_solver, nx
	real*8,  intent(in) :: acc_res
	real*8,  intent(in) :: x(nx)
	write(*,'(A,I5,A,I2,A,ES11.3,A,I2,A,ES22.15,A,ES22.15)')               &
	     ' (ieq_cell_report) cell ', j, ' class ', acc_class,              &
	     ' reaction residual ', acc_res, ' solver info ', info_solver,      &
	     ' x(H+) ', x(1), ' x(H2 nuclei) ', x(min(4,nx))
	end subroutine report_accepted_cell_state

	!----------------------------------!

	subroutine set_molecular_decay_cells(cells, ncell)
	! WHICH CELLS THE CHEMICAL-DECAY DIAGNOSTIC MEASURES AT. It has to be
	! said before the sweep runs: the reaction Jacobian can only be formed
	! where the cell's own rate coefficients are live, which is inside the
	! sweep's own loop, so a caller that wants a cell measured names it
	! first and reads the result afterwards.
	integer, intent(in) :: ncell, cells(ncell)
	integer :: i
	ieq_decay_cells = 0
	ieq_decay_ncell = 0
	ieq_decay_have  = .false.
	do i = 1, min(ncell, ieq_decay_ncell_max)
		ieq_decay_cells(i) = cells(i)
		ieq_decay_ncell    = i
	enddo
	end subroutine set_molecular_decay_cells

	!----------------------------------!

	subroutine molecular_decay_rate_of_the_cell(j, k, name, lam_row,       &
	                                            lam_slow, lam_fast, have)
	! WHAT THE DIAGNOSTIC MEASURED AT CELL j FOR MOLECULE k (1 H2, 2 H2+,
	! 3 H3+, 4 HeH+):
	!   lam_row   the molecule's own decay rate, -dA(k,k) of the local
	!             reaction system, attributable to that one molecule [s^-1];
	!   lam_slow  the smallest |Re lambda| of the whole local spectrum, the
	!             mode that decides whether the closure may eliminate the
	!             system at all [s^-1];
	!   lam_fast  the largest, which is the stiffness of the same system.
	! have is false where the cell was not armed or carries no molecular
	! network.
	integer,          intent(in)  :: j, k
	character(len=*), intent(out) :: name
	real*8,           intent(out) :: lam_row, lam_slow, lam_fast
	logical,          intent(out) :: have
	integer :: i, islot
	character(len=4), parameter :: mname(n_mol_decay_row) =                &
	                 (/ 'H2  ', 'H2+ ', 'H3+ ', 'HeH+' /)
	name    = ' '
	lam_row = 0.0d0;  lam_slow = 0.0d0;  lam_fast = 0.0d0
	have    = .false.
	if (k .lt. 1 .or. k .gt. n_mol_decay_row) return
	name  = mname(k)
	islot = 0
	do i = 1, ieq_decay_ncell
		if (ieq_decay_cells(i) .eq. j) islot = i
	enddo
	if (islot .eq. 0) return
	if (.not. ieq_decay_have(islot)) return
	lam_row  = ieq_decay_row (islot,k)
	lam_slow = ieq_decay_slow(islot)
	lam_fast = ieq_decay_fast(islot)
	have     = .true.
	end subroutine molecular_decay_rate_of_the_cell

	!----------------------------------!

	subroutine measure_molecular_decay_rates(j, x, nx, mbase, n_e_ref)
	! THE CHEMICAL DECAY RATES OF THE LOCAL REACTION SYSTEM AT ONE CELL,
	! with the elemental conservation modes removed.
	!
	! The local system is dot n = P(n) - L(n) for the species of the cell.
	! Written that way it has one null eigenvalue per conserved element, and
	! those modes are not decays: they are the statement that the reactions
	! move nuclei between species and make none. THE CODE'S OWN
	! PARAMETERIZATION HAS ALREADY REMOVED THEM. The unknowns are stage
	! fractions of a fixed nucleus total -- x(1) = n(H+)/n_H, x(4) =
	! 2 n(H2)/n_H, x(2) = n(He+)/n_He, a metal's x = n(X+)/n_X -- and the
	! neutral atomic stage of each element is the closure of that element's
	! budget rather than an unknown (System_HeH_mol: n_HI = (1 - x1 - x4 -
	! x5 - x6 - x7) n_H). So the Jacobian formed here in x IS the reaction
	! Jacobian on the subspace transverse to the conservation modes, and no
	! projection is applied or needed.
	!
	! THE MATRIX. The residual routine returns each molecular row divided by
	! its turnover scale (set_mol_turnover_rates), and that scale is the
	! row's nucleus total times a rate, s_k = N_k r_k. The rate matrix is
	!    A(k,l) = (c_k/N_k) (1/miu_k) dfv_k/dx_l         [s^-1]
	! with miu_k = mol_inv_turnover(k) and c_k the stoichiometric factor of
	! the fraction the row balances (1 for a stage fraction, 2 for x(H2) and
	! x(H2+), 3 for x(H3+)), because fv_k = miu_k dn_k/dt and
	! dot x_k = c_k (dn_k/dt)/N_k. dfv/dx is a central difference of the
	! production routine the solve itself uses.
	!
	! THE IMPOSED ROWS ARE LIFTED FOR THE MEASUREMENT. Where a partition is
	! transported its balance row is replaced by x - x_fix, an identity with
	! no chemistry in it; the question this diagnostic asks is what the
	! chemistry of that partition WOULD do, which is exactly the number that
	! decides whether it has to be transported. The flags are put back
	! before the routine returns.
	!
	! IDENTITY ROWS OF ABSENT ELEMENTS ARE LEFT OUT of the spectrum: their
	! row is x itself and its eigenvalue would be a property of the writing
	! and not of any reaction.
	!
	! Called from inside the sweep's parallel loop for the cells named by
	! set_molecular_decay_cells, under the sweep's report lock.
	integer, intent(in) :: j, nx, mbase
	real*8,  intent(in) :: x(nx), n_e_ref
	real*8  :: xp(nx), fvp(nx), fvm(nx), par(60), h
	real*8  :: A(nx,nx), wr(nx), wi(nx), vdum(1,1)
	real*8  :: work(8*nx+64), conv(nx), lam, lo, hi
	integer :: iflag, i, l, e, ix, islot, info, lwork, nsub, ksub(nx)
	logical :: keep(nx), h2_fix, hp_fix, ox_fix
	if (.not. thereis_mol) return
	! The molecular block is seven rows before any option adds to it and
	! the rate conversion below names all seven; a shorter system is not
	! this network.
	if (nx .lt. 7) return
	islot = 0
	do i = 1, ieq_decay_ncell
		if (ieq_decay_cells(i) .eq. j) islot = i
	enddo
	if (islot .eq. 0) return

	! The conversion of each row into a rate: c_k/N_k times the turnover
	! scale the residual routine already divided the row by.
	conv = 0.0d0
	keep = .false.
	conv(1) = 1.0d0/max(ieq_cell%nh,  1.0d-300)
	conv(2) = 1.0d0/max(ieq_cell%nhe, 1.0d-300)
	conv(3) = conv(2)
	conv(4) = 2.0d0/max(ieq_cell%nh,  1.0d-300)
	conv(5) = conv(4)
	conv(6) = 3.0d0/max(ieq_cell%nh,  1.0d-300)
	conv(7) = 1.0d0/max(ieq_cell%nh,  1.0d-300)
	keep(1:min(7,nx)) = .true.
	if (thereis_HeITR .and. nx .ge. 8) then
		conv(8) = conv(2)
		keep(8) = .true.
	endif
	if (thereis_oxychem) then
		ix = oxygen_row_base()
		if (ix+1 .le. nx) then
			conv(ix)   = 1.0d0/max(ieq_cell%n_ofam, 1.0d-300)
			conv(ix+1) = conv(ix)
			keep(ix)   = .true.
			keep(ix+1) = .true.
		endif
	endif
	if (thereis_metals) then
		do e = 1, met_nelem
			ix = mbase + 2*(e-1)
			if (ix+1 .gt. nx) cycle
			if (met_ntot(e) .le. 1.0d-30) cycle
			conv(ix) = 1.0d0/met_ntot(e)
			keep(ix) = .true.
			if (met_top(e) .ge. 2) then
				conv(ix+1) = conv(ix)
				keep(ix+1) = .true.
			endif
		enddo
	endif
	do i = 1, nx
		if (mol_inv_turnover(i) .gt. 0.0d0)                            &
			conv(i) = conv(i)/mol_inv_turnover(i)
	enddo

	h2_fix = ieq_cell%x_h2_fixed
	hp_fix = ieq_cell%x_hp_fixed
	ox_fix = ieq_cell%x_ox_fixed
	ieq_cell%x_h2_fixed = .false.
	ieq_cell%x_hp_fixed = .false.
	ieq_cell%x_ox_fixed = .false.

	par   = 0.0d0
	iflag = 1
	A     = 0.0d0
	do l = 1, nx
		! A step of the size of the unknown, with a floor so that a stage
		! sitting at none is still sampled.
		h = max(1.0d-7*abs(x(l)), 1.0d-14)
		xp = x;  xp(l) = x(l) + h
		call molecular_network_rows(nx, xp, fvp, iflag, par, mbase)
		xp = x;  xp(l) = x(l) - h
		call molecular_network_rows(nx, xp, fvm, iflag, par, mbase)
		do i = 1, nx
			A(i,l) = conv(i)*(fvp(i) - fvm(i))/(2.0d0*h)
		enddo
	enddo

	ieq_cell%x_h2_fixed = h2_fix
	ieq_cell%x_hp_fixed = hp_fix
	ieq_cell%x_ox_fixed = ox_fix

	! The molecules' own rates: the diagonal of the rate matrix, which is
	! the rate at which that molecule's own population relaxes when nothing
	! else moves.
	do i = 1, n_mol_decay_row
		ieq_decay_row(islot,i) = -A(3+i,3+i)
	enddo

	! The spectrum, over the rows that carry a reaction.
	nsub = 0
	do i = 1, nx
		if (.not. keep(i)) cycle
		nsub = nsub + 1
		ksub(nsub) = i
	enddo
	if (nsub .le. 0) return
	do i = 1, nsub
		do l = 1, nsub
			A(i,l) = A(ksub(i),ksub(l))
		enddo
	enddo
	lwork = 8*nx + 64
	call dgeev('N','N', nsub, A(1:nsub,1:nsub), nsub, wr, wi, vdum, 1,     &
	           vdum, 1, work, lwork, info)
	if (info .ne. 0) return
	lo =  huge(1.0d0);  hi = 0.0d0
	do i = 1, nsub
		lam = abs(wr(i))
		if (lam .le. 0.0d0) cycle
		lo = min(lo, lam);  hi = max(hi, lam)
	enddo
	if (hi .le. 0.0d0) return
	ieq_decay_slow(islot) = lo
	ieq_decay_fast(islot) = hi
	ieq_decay_have(islot) = .true.
	end subroutine measure_molecular_decay_rates

	!----------------------------------!

	subroutine molecular_network_rows(n, x, fv, iflag, par, mbase)
	! The molecular reaction system of this cell, with or without the metal
	! block, in ONE call: the same two routines the solve and
	! normalized_reaction_residual choose between, so a derivative taken
	! here is a derivative of the system the cell was solved on.
	integer, intent(in)    :: n, mbase
	integer, intent(inout) :: iflag
	real*8,  intent(in)    :: x(n)
	real*8,  intent(inout) :: par(60)
	real*8,  intent(out)   :: fv(n)
	if (thereis_metals) then
		cx_metal_base = mbase
		call ion_system_HeH_mol_metals(n, x, fv, iflag, par)
		cx_metal_base = 4
	else
		call ion_system_HeH_mol(n, x, fv, iflag, par)
	endif
	end subroutine molecular_network_rows

	! ------------------------------------------------------!

	pure integer function nonroot_acceptance_class(res) result(c)
	! WHICH OF THE TWO NON-ROOT CLASSES A CANDIDATE FALLS IN, and the single
	! definition of that rule. Class 4 is the relaxation amnesty: the
	! candidate is adopted, reported and counted. Class 6 is the amnesty's
	! severity bound: the candidate is NOT adopted and the cell keeps the
	! composition it entered the sweep with.
	! The test is written as .not.(res <= cap) so that a residual which is
	! not a finite number falls on the class-6 side. A state that has driven
	! a balance row out of the reals is not one to march on either, and NaN
	! fails every ordinary comparison, so the ordinary form would have made
	! it the mildest case instead of the worst.
	real*8, intent(in) :: res
	if (.not. (res .le. ieq_nonroot_res_cap)) then
		c = 6
	else
		c = 4
	endif
	end function nonroot_acceptance_class

	!----------------------------------!

	subroutine nonroot_streak_update(acc_class,j,step,radius,T,dens,n_e,   &
	                                 info,viol,res,x,n)
	! Relaxation amnesty for a NON-ROOT acceptance, and its limit. A cell
	! whose accepted state is a root (classes 1-3, and class 5, the root of
	! the constrained element-conserving continuation solve) resets its
	! streak; a cell whose state is NOT a root -- class 4, the non-root
	! adopted under the amnesty, and class 6, the non-root above
	! ieq_nonroot_res_cap for which the cell kept the composition it
	! entered the sweep with -- is reported loudly, counted, and
	! allowed to continue for at most ieq_nonroot_streak_stop consecutive
	! sweeps -- the measured signature of a cold-start transient that the
	! next sweeps repair (see ieq_res_tol / ieq_nonroot_streak_stop above).
	! One sweep beyond that the same cell resting on a non-root is a
	! solution built on a state that does not solve the reaction network:
	! continuing would produce a physically meaningless wind
	! (docs/supersonic_molecular_base.md section 11.5-A), so the run stops
	! here with the full cell diagnostics. The two classes share the streak
	! because they are the same statement about the cell -- its chemistry
	! was not solved -- and differ only in which state was kept.
	integer, intent(in) :: acc_class, j, step, info, n
	real*8,  intent(in) :: radius, T, dens, n_e, viol, res
	real*8,  intent(in) :: x(n)

	! A steady-solver probe is not a sweep the run rested on: it is one
	! sample of a directional derivative, or a line-search trial that may be
	! thrown away in the next statement. It neither raises the streak nor
	! clears it, and it is not reported cell by cell -- one measured JFNK
	! phase produced 18,339 non-root probe acceptances against none in the
	! 2002 marching steps before it (section 121). What the streak counts is
	! consecutive sweeps of a state the run HOLDS.
	if (ieq_sweep_state_kind .eq. ieq_state_steady_candidate) return

	! Classes 4 and 6 are the non-root outcomes; every other class, class 5
	! included, is a root and clears the streak.
	if (acc_class .ne. 4 .and. acc_class .ne. 6) then
		ieq_nonroot_streak(j) = 0
		return
	endif

	ieq_nonroot_streak(j) = ieq_nonroot_streak(j) + 1
	if (acc_class .eq. 6) then
		call report_acceptance_event(                                     &
			'NON-ROOT above the amnesty cap, entry composition kept',     &
			j,step,radius,T,info,viol,res)
	else
		call report_acceptance_event(                                     &
			'NON-ROOT accepted (relaxation amnesty)',                     &
			j,step,radius,T,info,viol,res)
	endif

	if (ieq_nonroot_streak(j) .ge. ieq_nonroot_streak_stop) then
		!$omp critical (ieq_acc_report)
		write(*,'(A)') ' (ioniz_eq) STOP: a cell has rested on a '//      &
			'NON-ROOT chemical equilibrium beyond the relaxation amnesty'
		write(*,'(A,I0,A,I0,A,I0,A)') '   cell ', j, '  step ', step,     &
			'  consecutive non-root sweeps ', ieq_nonroot_streak(j), ''
		write(*,'(A,F0.6,A,ES11.4,A,ES11.4,A,ES11.4)')                    &
			'   r [R_p] ', radius, '  T [K] ', T,                         &
			'  n_tot [cm^-3] ', dens, '  n_e [cm^-3] ', n_e
		write(*,'(A,I0,A,ES10.3,A,ES10.3,A,ES10.3)')                      &
			'   solver info ', info, '  element violation ', viol,        &
			'  normalized reaction residual ', res,                       &
			'  (tolerance ', ieq_res_tol, ')'
		write(*,'(A)') '   candidate stage fractions:'
		write(*,'(6ES12.4)') x(1:n)
		write(*,'(A)') '   Bounds and element conservation are necessary'//&
			' conditions, not reaction equilibrium; continuing with a'
		write(*,'(A)') '   known non-root produces a physically'//        &
			' meaningless wind. See docs/Update_EXHALE_stage1.md section 113.'
		flush(6)
		!$omp end critical (ieq_acc_report)
		error stop 'ioniz_eq: persistent non-root chemical equilibrium'
	endif

	end subroutine nonroot_streak_update

	!----------------------------------!

	subroutine report_acceptance_event(label,j,step,radius,T,info,viol,res)
	! One line per NON-ROOT acceptance (nonroot_streak_update), with the
	! cell, step, radius, temperature, solver exit code, element-budget
	! violation and normalized reaction residual -- a non-root acceptance
	! is never silent. Capped at ieq_acc_nprint_max lines per run so a
	! pathological run cannot flood the log; the run-wide counters keep
	! the full population.
	character(len=*), intent(in) :: label
	integer, intent(in) :: j, step, info
	real*8,  intent(in) :: radius, T, viol, res

	! Only states the run holds. A probe's acceptance is a property of a
	! state the solver invented to differentiate at, and printing it -- and
	! spending the run's print budget on it -- hides the run's own events.
	if (ieq_sweep_state_kind .eq. ieq_state_steady_candidate) return

	!$omp critical (ieq_acc_report)
	if (ieq_acc_nprint .lt. ieq_acc_nprint_max) then
		ieq_acc_nprint = ieq_acc_nprint + 1
		write(*,'(A,A,A,I0,A,I0,A,F0.6,A,ES10.3,A,I0,A,ES10.3,A,ES10.3)') &
			' (ioniz_eq) ', label, ': cell ', j, ' step ', step,          &
			' r ', radius, ' T[K] ', T, ' info ', info,                   &
			' viol ', viol, ' res ', res
		if (ieq_acc_nprint .eq. ieq_acc_nprint_max)                       &
			write(*,'(A)') ' (ioniz_eq) further acceptance report '//     &
				'lines suppressed (run-wide counters keep counting)'
	endif
	!$omp end critical (ieq_acc_report)
	end subroutine report_acceptance_event

	!----------------------------------!

	subroutine clamp_fractions_to_element_budget(x,n,mbase)
	! Move a root that lies just outside the physically allowed states onto the
	! nearest allowed one: no negative populations, and no element with more
	! nuclei in its tracked stages than it has. Negative fractions are set to
	! zero and, where an element's tracked stages still sum above its nuclei,
	! they are scaled down to sum exactly to them, leaving the neutral stage
	! empty.
	!
	! This is the treatment of a root the solver has resolved to a face of the
	! allowed region -- a fully dissociated H2, a fully neutral or fully
	! ionized element -- and delivered as a small number of either sign. The
	! clamped state is that root to the accuracy the solve reached, so it
	! remains a solution of the equilibrium at THIS cell's state; what it must
	! not be replaced by is a composition carried over from an earlier
	! evaluation.

	integer, intent(in)    :: n, mbase
	real*8,  intent(inout) :: x(n)
	integer :: im, ix, iy
	real*8  :: s, hehp_he, ox_h

	do ix = 1,n
		if (x(ix) .lt. 0.0d0) x(ix) = 0.0d0
	enddo

	! Each metal element separately.  The oxygen element carries its two
	! molecular carriers in the same budget and is scaled with them, so that
	! the clamped state still has a non-negative free atomic oxygen. This
	! block runs FIRST because the oxygen carriers also hold H nuclei, and
	! the hydrogen budget below has to see them already capped. It is
	! disjoint from the hydrogen and helium unknowns, so an oxygen-chemistry
	! -free run is unaffected by the reordering.
	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			s  = x(ix)
			if (melem_top(im) .ge. 2) s = s + x(ix+1)
			iy = 0
			if (thereis_oxychem .and. im .eq. iel_O) then
				iy = oxygen_row_base()
				s  = s + x(iy) + x(iy+1)
			endif
			if (s .gt. 1.0d0) then
				x(ix) = x(ix)/s
				if (melem_top(im) .ge. 2) x(ix+1) = x(ix+1)/s
				if (iy .gt. 0) then
					x(iy)   = x(iy)/s
					x(iy+1) = x(iy+1)/s
				endif
			endif
		enddo
	endif

	! Hydrogen nuclei: H II, the H bound in molecules where tracked, and the
	! H bound in the oxygen carriers. ox_h is that last part, expressed per
	! H nucleus; it has just been capped against the oxygen budget, so it is
	! taken as given here and the H species are scaled into what it leaves,
	! exactly as HeH+ is treated in the helium budget below.
	ox_h = 0.0d0
	if (thereis_oxychem .and. ieq_cell%nh .gt. 0.0d0) then
		iy   = oxygen_row_base()
		ox_h = (x(iy) + 2.0d0*x(iy+1))*ieq_cell%n_ofam/ieq_cell%nh
		if (ox_h .gt. 1.0d0) then
			x(iy)   = x(iy)/ox_h
			x(iy+1) = x(iy+1)/ox_h
			ox_h    = 1.0d0
		endif
	endif
	s = x(1)
	if (thereis_mol) s = s + x(4) + x(5) + x(6) + x(7)
	if (s .gt. 1.0d0 - ox_h) then
		x(1) = x(1)*(1.0d0 - ox_h)/s
		if (thereis_mol) then
			do ix = 4,7
				x(ix) = x(ix)*(1.0d0 - ox_h)/s
			enddo
		endif
	endif

	! Helium nuclei: He II, He III, the 2^3S metastable and, in the molecular
	! layout, the He nucleus bound in HeH+. x(7) is HeH+ per H nucleus and
	! belongs to BOTH element budgets; it has just been capped against the H
	! nuclei, so it is capped against the He nuclei here as well (which can
	! only relax the H budget), and the free-He stages are scaled into what
	! is left. hehp_he is exactly zero without molecules, where the scaling
	! reduces to the plain x/s of every non-molecular run.
	if (thereis_He) then
		hehp_he = 0.0d0
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0) then
			hehp_he = x(7)*ieq_cell%nh/ieq_cell%nhe
			if (hehp_he .gt. 1.0d0) then
				x(7)    = x(7)/hehp_he
				hehp_he = 1.0d0
			endif
		endif
		s = x(2) + x(3)
		ix = 0
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			s = s + x(ix)
		endif
		if (s .gt. 1.0d0 - hehp_he) then
			x(2) = x(2)*(1.0d0 - hehp_he)/s
			x(3) = x(3)*(1.0d0 - hehp_he)/s
			if (ix .gt. 0) x(ix) = x(ix)*(1.0d0 - hehp_he)/s
		endif
	endif

	end subroutine clamp_fractions_to_element_budget

	!----------------------------------!

	subroutine dissociation_ionization_balance_at_fixed_ne(x,n,mbase,n_e,   &
	                                                       p_bar,T_gas)
	! H2 dissociation equilibrium at the local (p, T), with every element in
	! its own ionization balance at a given electron density: the limit of the
	! molecular network in which the couplings between the molecular ions and
	! the ionization of H and He are dropped.
	!
	!   - The H nuclei are partitioned between H2 and atomic H by the
	!     chemical-equilibrium mixing ratio q_H2(p, T) of Koskinen et al.
	!     (2022) Eq. 11, the same fit the molecular base boundary condition
	!     uses; per H nucleus the fraction bound in H2 is
	!     2 q (1 + He/H)/(1 + q) (lower_column::mu_mixture).
	!   - The atomic remainder is ionized in the H ionization balance of
	!     ionization_balance_at_fixed_ne, which also sets He, the 2^3S
	!     metastable and the metal stages.
	!   - H2+, H3+ and HeH+ are left at zero: each is a trace intermediate
	!     whose abundance is set by the very couplings this limit drops, and
	!     each holds far fewer H nuclei than H2 wherever the gas is molecular.
	!
	! Every fraction returned is non-negative and each element's stages sum to
	! at most its nuclei, so the starting point is a physically allowed state,
	! and it depends only on the local (p, T, n_e) and on the rate
	! coefficients of the cell -- not on the cell's previous composition.
	!
	! Molecular layout only: x(4) is the H-nucleus fraction bound in H2 there,
	! and the He 2^3S metastable sits at x(8).

	integer, intent(in)  :: n, mbase
	real*8,  intent(in)  :: n_e, p_bar, T_gas
	real*8,  intent(out) :: x(n)
	integer :: iy, ix
	real*8 :: qh2, x_h2, nh2s, nhis, f_oh, f_h2o, s_free

	call ionization_balance_at_fixed_ne(x,n,mbase,n_e)

	qh2  = q_h2_equilibrium(p_bar, T_gas)
	x_h2 = 2.0d0*qh2*(1.0d0 + HeH)/(1.0d0 + qh2)
	if (x_h2 .gt. 1.0d0) x_h2 = 1.0d0
	x(4) = x_h2
	! Only the H nuclei left atomic are available to ionize.
	x(1) = x(1)*(1.0d0 - x_h2)

	! Oxygen carriers at the chemical equilibrium of the same partition:
	! O + H2 <-> OH + H and OH + H2 <-> H2O + H at this cell's H2/H ratio.
	! Zero is a valid root of the water cycle and hybrd1 is bistable from a
	! zero molecular seed, so the molecular basin has to be seeded with a
	! molecular oxygen partition as well, not only a molecular hydrogen one.
	! The free atomic oxygen then closes the oxygen budget, and the metal
	! stage fractions ionization_balance_at_fixed_ne wrote are scaled into
	! what is left of the element.
	if (thereis_oxychem) then
		iy   = oxygen_row_base()
		nh2s = 0.5d0*x_h2*ieq_cell%nh
		nhis = max(1.0d0 - x_h2, 0.0d0)*(1.0d0 - x(1))*ieq_cell%nh
		call oxygen_chemical_equilibrium_fractions(T_gas, nh2s, nhis,     &
		                                           f_oh, f_h2o)
		x(iy)   = f_oh
		x(iy+1) = f_h2o
		s_free  = max(1.0d0 - f_oh - f_h2o, 0.0d0)
		ix = mbase + 2*(iel_O-1)
		x(ix)   = x(ix)*s_free
		x(ix+1) = x(ix+1)*s_free
	endif

	end subroutine dissociation_ionization_balance_at_fixed_ne

	!----------------------------------!

	subroutine metal_ionization_balance_at_fixed_ne(x,n,mbase,n_e)
	! Metal stage fractions in each element's OWN ionization balance at a
	! given electron density, element by element in canonical order (the
	! metal part of ionization_balance_at_fixed_ne above; see there for the
	! balance itself). Written separately because the molecular retry needs
	! exactly this part: its H/He/molecular unknowns come from the
	! chemical-equilibrium H2 fit, its metal unknowns from here. Only
	! x(mbase..) is touched, so the caller's other seeds are preserved. An
	! element that is absent keeps its stages at zero, as its residual rows do.

	use System_HeH_metals, only: met_nelem, met_ntot, met_g0, met_g1,     &
	                             met_b0, met_b1, met_a1, met_a2, met_top

	integer, intent(in)    :: n, mbase
	real*8,  intent(in)    :: n_e
	real*8,  intent(inout) :: x(n)
	integer :: im, ix, top
	real*8  :: u0,u1,d1,d2,w1,w2,s

	do im = 1,met_nelem
		ix = mbase + 2*(im-1)
		x(ix)   = 0.0d0
		x(ix+1) = 0.0d0
		if (met_ntot(im) .le. 1.0d-30) cycle
		top = met_top(im)
		u0  = met_g0(im) + met_b0(im)*n_e
		d1  = met_a1(im)*n_e
		if (top .ge. 2) then
			u1 = met_g1(im) + met_b1(im)*n_e
			d2 = met_a2(im)*n_e
			w1 = u0*d2
			w2 = u0*u1
			s  = d1*d2 + w1 + w2
			if (s .gt. 0.0d0) then
				x(ix)   = w1/s
				x(ix+1) = w2/s
			endif
		else
			s = d1 + u0
			if (s .gt. 0.0d0) x(ix) = u0/s
		endif
	enddo

	end subroutine metal_ionization_balance_at_fixed_ne

	! End of module

	!----------------------------------!
	! THE ELIMINATED-SPECIES CLOSURE, MEASURED ON A STATE
	!
	! B1a section 2.4: every run eliminates the neutral stage of each
	! element and the electron density by a closure, and until now that
	! closure was evaluated only inside the sweep that also rewrites f_sp
	! and rho. What follows measures it WITHOUT solving. The rate state of
	! each cell is the one the sweep built from this state's own radiation
	! field (ieq_rate_cell); the composition is read from the caller's
	! f_sp; the number returned is normalized_reaction_residual, the same
	! measure the sweep judges a root by. A cell the sweep accepted as a
	! root -- acceptance class 1, 2, 3 or 5 -- therefore reads at or below
	! ieq_res_tol, and a class-4 cell above it.
	!
	! NOTHING IS WRITTEN: not f_sp, not rho, not a ledger, not a streak.
	! The scratch of one cell that these routines install -- ieq_cell, the
	! metal coefficients of System_HeH_metals, the molecular and oxygen
	! coefficients and the row turnover scales of System_HeH_mol -- is held
	! aside and put back, and the round trip is asserted rather than
	! assumed (ieq_closure_scratch_matches).
	!
	! ONE EXCEPTION, STATED. cx_kc, the charge-exchange rate coefficients of
	! a cell, is private to charge_exchange, which exports no accessor for
	! it: it can be re-set from a temperature (cx_set_cell) but not read.
	! It is scratch that every cell of every sweep overwrites from its own
	! T before any residual reads it, so its content between
	! sweeps has no meaning; the loop below leaves it holding the
	! coefficients of the last cell it measured, exactly as the sweep
	! leaves it holding those of the last cell one of its threads ran.

	subroutine save_ieq_closure_scratch(ws)
	type(ieq_closure_scratch), intent(out) :: ws
	ws%cell = ieq_cell
	ws%met_alloc = allocated(met_ntot)
	if (ws%met_alloc) then
		ws%nelem = met_nelem
		ws%ntot(1:met_nelem) = met_ntot(1:met_nelem)
		ws%g0(1:met_nelem)   = met_g0(1:met_nelem)
		ws%g1(1:met_nelem)   = met_g1(1:met_nelem)
		ws%b0(1:met_nelem)   = met_b0(1:met_nelem)
		ws%b1(1:met_nelem)   = met_b1(1:met_nelem)
		ws%a1(1:met_nelem)   = met_a1(1:met_nelem)
		ws%a2(1:met_nelem)   = met_a2(1:met_nelem)
		ws%top(1:met_nelem)  = met_top(1:met_nelem)
	endif
	ws%molc = pack_mol_coefficients()
	ws%minv = mol_inv_turnover
	end subroutine save_ieq_closure_scratch

	!----------------------------------!

	subroutine restore_ieq_closure_scratch(ws)
	type(ieq_closure_scratch), intent(in) :: ws
	ieq_cell = ws%cell
	if (ws%met_alloc) then
		met_nelem = ws%nelem
		met_ntot(1:ws%nelem) = ws%ntot(1:ws%nelem)
		met_g0(1:ws%nelem)   = ws%g0(1:ws%nelem)
		met_g1(1:ws%nelem)   = ws%g1(1:ws%nelem)
		met_b0(1:ws%nelem)   = ws%b0(1:ws%nelem)
		met_b1(1:ws%nelem)   = ws%b1(1:ws%nelem)
		met_a1(1:ws%nelem)   = ws%a1(1:ws%nelem)
		met_a2(1:ws%nelem)   = ws%a2(1:ws%nelem)
		met_top(1:ws%nelem)  = ws%top(1:ws%nelem)
	endif
	call unpack_mol_coefficients(ws%molc)
	mol_inv_turnover = ws%minv
	end subroutine restore_ieq_closure_scratch

	!----------------------------------!

	logical function ieq_closure_scratch_matches(ws) result(same)
	! Bit comparison of the live scratch against the copy taken before a
	! measurement. A measurement that changed the state it measured would
	! move the next cell solve by an amount nothing else would explain.
	type(ieq_closure_scratch), intent(in) :: ws
	type(ieq_closure_scratch) :: now
	call save_ieq_closure_scratch(now)
	same = all(now%molc .eq. ws%molc) .and. all(now%minv .eq. ws%minv)   &
	       .and. (now%met_alloc .eqv. ws%met_alloc)
	if (same .and. ws%met_alloc) then
		same = (now%nelem .eq. ws%nelem)                                  &
		       .and. all(now%ntot .eq. ws%ntot)                           &
		       .and. all(now%g0 .eq. ws%g0) .and. all(now%g1 .eq. ws%g1)  &
		       .and. all(now%b0 .eq. ws%b0) .and. all(now%b1 .eq. ws%b1)  &
		       .and. all(now%a1 .eq. ws%a1) .and. all(now%a2 .eq. ws%a2)  &
		       .and. all(now%top .eq. ws%top)
	endif
	if (same) same = ion_rates_match(now%cell, ws%cell)
	end function ieq_closure_scratch_matches

	!----------------------------------!

	logical function ion_rates_match(a, b) result(same)
	! The cell rate state compared through the quantities the closure
	! evaluator installs; every one of them is a number the residual rows
	! read.
	type(ion_rates), intent(in) :: a, b
	same = (a%P_HI .eq. b%P_HI) .and. (a%P_HeI .eq. b%P_HeI)             &
	  .and. (a%P_HeII .eq. b%P_HeII) .and. (a%P_HeITR .eq. b%P_HeITR)    &
	  .and. (a%rchiiB .eq. b%rchiiB) .and. (a%rcheiiB .eq. b%rcheiiB)    &
	  .and. (a%rcheiiiB .eq. b%rcheiiiB) .and. (a%rcheiTR .eq. b%rcheiTR)&
	  .and. (a%nh .eq. b%nh) .and. (a%nhe .eq. b%nhe)                    &
	  .and. (a%a_ion_HI .eq. b%a_ion_HI)                                 &
	  .and. (a%a_ion_HeI .eq. b%a_ion_HeI)                               &
	  .and. (a%a_ion_HeII .eq. b%a_ion_HeII)                             &
	  .and. (a%a_ion_HeITR .eq. b%a_ion_HeITR)                           &
	  .and. (a%q13 .eq. b%q13) .and. (a%q31a .eq. b%q31a)                &
	  .and. (a%q31b .eq. b%q31b) .and. (a%Q31 .eq. b%Q31)                &
	  .and. (a%A31 .eq. b%A31)                                           &
	  .and. (a%kcx_He0_Hp .eq. b%kcx_He0_Hp)                             &
	  .and. (a%kcx_Hep_H0 .eq. b%kcx_Hep_H0)                             &
	  .and. (a%P_H2 .eq. b%P_H2) .and. (a%k_LW .eq. b%k_LW)              &
	  .and. (a%T_K .eq. b%T_K) .and. (a%ntot .eq. b%ntot)                &
	  .and. (a%n_ofam .eq. b%n_ofam) .and. (a%n_co .eq. b%n_co)
	end function ion_rates_match

	!----------------------------------!

	function pack_mol_coefficients() result(c)
	! The molecular and oxygen rate coefficients of System_HeH_mol in one
	! vector, so that they can be held aside and put back by value. The
	! order is the declaration order of that module.
	real*8 :: c(27)
	c = [ mk5, mk6, mk7, mk8, mk9, mk10, mk11, mk12, mk13, mk14, mk15,   &
	      mk16, mk17, mk18, mk19, mk_h2p_he, mk23, mk_ion_H2,            &
	      ok1, ok1r, ok2, ok2r, ok6, oj3, oj4, oj5, oj7 ]
	end function pack_mol_coefficients

	!----------------------------------!

	subroutine unpack_mol_coefficients(c)
	real*8, intent(in) :: c(27)
	mk5  = c(1);  mk6  = c(2);  mk7  = c(3);  mk8  = c(4)
	mk9  = c(5);  mk10 = c(6);  mk11 = c(7);  mk12 = c(8)
	mk13 = c(9);  mk14 = c(10); mk15 = c(11); mk16 = c(12)
	mk17 = c(13); mk18 = c(14); mk19 = c(15); mk_h2p_he = c(16)
	mk23 = c(17); mk_ion_H2 = c(18)
	ok1  = c(19); ok1r = c(20); ok2  = c(21); ok2r = c(22)
	ok6  = c(23); oj3  = c(24); oj4  = c(25); oj5  = c(26)
	oj7  = c(27)
	end subroutine unpack_mol_coefficients

	!----------------------------------!

	integer function ieq_triplet_row() result(ir)
	! Row of the He 2^3S level balance in the system the run solves
	! (tr_triplet_row of ion_residual_core): x(4) in the atomic layouts,
	! x(8) in the molecular one. Zero where the level is not carried.
	ir = 0
	if (.not. thereis_He) return
	if (.not. thereis_HeITR) return
	ir = 4
	if (thereis_mol) ir = 8
	end function ieq_triplet_row

	!----------------------------------!

	subroutine ieq_install_cell(j)
	! Put cell j back into the scratch the residual rows read: the rate
	! state the sweep built for it, its metal coefficients, and the
	! molecular, oxygen and turnover coefficients rebuilt from the cell's
	! own T, particle count and electron density by the SAME routines the
	! sweep calls, so there is no second copy of any of them.
	integer, intent(in) :: j
	ieq_cell = ieq_rate_cell(j)
	if (thereis_mol) then
		call set_mol_coeffs(ieq_TK_cell(j), ieq_ntot_cell(j))
		if (thereis_oxychem)                                              &
			call set_oxygen_coeffs(ieq_TK_cell(j), j_h2o_fuv(j,:),        &
			                                       j_oh_fuv(j,:))
		call set_mol_turnover_rates(ieq_ne_cell(j), 1.0d0)
		if (thereis_oxychem)                                              &
			call set_oxygen_turnover_rates(ieq_iox_stored, 1.0d0)
	endif
	if (thereis_metals) then
		call cx_set_cell(ieq_TK_cell(j))
		call set_metal_coeffs(n_melem, ieq_met_coef(:,1,j),               &
		                      ieq_met_coef(:,2,j), ieq_met_coef(:,3,j),   &
		                      ieq_met_coef(:,4,j), ieq_met_coef(:,5,j),   &
		                      ieq_met_coef(:,6,j), ieq_met_coef(:,7,j),   &
		                      melem_top)
		if (thereis_mol)                                                  &
			call set_mol_metal_turnover_rates(ieq_ne_cell(j), 1.0d0)
	endif
	end subroutine ieq_install_cell

	!----------------------------------!

	subroutine ieq_state_vector(j, rho, f_sp, x, nx)
	! The composition of cell j in the fraction layout of the system the
	! run solves, built from the caller's state. The denominators are the
	! ELEMENT TOTALS the cell's rate state carries, which is what the rows
	! of that system divide by, so the vector and the rates describe one
	! cell.
	integer, intent(in) :: j
	real*8,  intent(in) :: rho(1-Ng:N+Ng)
	real*8,  intent(in) :: f_sp(1-Ng:N+Ng,n_species)
	real*8,  intent(out) :: x(:)
	integer, intent(out) :: nx
	real*8  :: nd, nHl, nHel, nOf, nXl
	integer :: im, i0, itr
	nx = ieq_neq_stored
	x(1:nx) = 0.0d0
	nd   = rho(j)*n0
	nHl  = max(ieq_rate_cell(j)%nh,  1.0d-99)
	nHel = max(ieq_rate_cell(j)%nhe, 1.0d-99)
	x(1) = f_sp(j,2)*nd/nHl
	if (.not. thereis_He) return
	x(2) = f_sp(j,4)*nd/nHel
	x(3) = f_sp(j,5)*nd/nHel
	itr  = ieq_triplet_row()
	if (thereis_mol) then
		x(4) = 2.0d0*f_sp(j,isp_H2) *nd/nHl
		x(5) = 2.0d0*f_sp(j,isp_H2p)*nd/nHl
		x(6) = 3.0d0*f_sp(j,isp_H3p)*nd/nHl
		x(7) =       f_sp(j,isp_HeHp)*nd/nHl
		if (thereis_oxychem) then
			nOf = max(ieq_rate_cell(j)%n_ofam, 1.0d-99)
			x(ieq_iox_stored)   = f_sp(j,isp_OH) *nd/nOf
			x(ieq_iox_stored+1) = f_sp(j,isp_H2O)*nd/nOf
		endif
	endif
	if (itr .gt. 0) x(itr) = f_sp(j,6)*nd/nHel
	if (thereis_metals) then
		do im = 1, n_melem
			i0  = melem_i0(im)
			nXl = max(ieq_met_coef(im,1,j), 1.0d-99)
			x(ieq_mbase_stored+2*(im-1)) = f_sp(j,mion_fsp(i0+1))*nd/nXl
			if (melem_top(im) .ge. 2)                                     &
				x(ieq_mbase_stored+1+2*(im-1)) =                          &
					f_sp(j,mion_fsp(i0+2))*nd/nXl
		enddo
	endif
	end subroutine ieq_state_vector

	!----------------------------------!

	subroutine ionization_closure_residual_cell(j, x, nx, res, rowres)
	! The normalized reaction residual of ONE cell at a GIVEN composition,
	! with that cell's rate state installed and the scratch put back. It
	! solves nothing and writes no state.
	integer, intent(in)  :: j, nx
	real*8,  intent(in)  :: x(nx)
	real*8,  intent(out) :: res
	real*8,  intent(out) :: rowres(nx)
	type(ieq_closure_scratch) :: ws
	call save_ieq_closure_scratch(ws)
	call ieq_install_cell(j)
	res = normalized_reaction_residual(x, nx, ieq_mbase_stored,          &
	                                   ieq_ne_cell(j), rowres)
	call restore_ieq_closure_scratch(ws)
	if (.not. ieq_closure_scratch_matches(ws))                           &
		write(*,'(A)') ' (closure residual) WARNING: the cell scratch '// &
		     'was NOT reinstated after the isolated evaluation'
	end subroutine ionization_closure_residual_cell

	!----------------------------------!

	subroutine ionization_closure_residual_profile(rho, f_sp, res,       &
	                                               res_triplet, ok, why)
	! The closure residual of every physical cell of a state, and the He
	! 2^3S level row of the same evaluation (B1a section 2.5: the
	! metastable is an unknown of this system, so its balance is one row of
	! it). res_triplet is zero where the level is not carried.
	!
	! ok is false, with the reason in why, when no sweep has left a cell
	! state to measure against, or when the ELEMENT TOTALS of the state
	! handed in are not the ones those rates were built with -- a different
	! gas, whose closure these rates cannot state.
	real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho
	real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
	real*8, dimension(1:N),                 intent(out) :: res, res_triplet
	logical,                                intent(out) :: ok
	character(len=*),                       intent(out) :: why
	type(ieq_closure_scratch) :: ws
	real*8, dimension(1-Ng:N+Ng) :: nhil, nhiil, nheil, nheiil, nheiiil
	real*8, dimension(1-Ng:N+Ng) :: nHl, nHel
	real*8, dimension(1-Ng:N+Ng,4) :: nmoll
	real*8, dimension(1-Ng:N+Ng,3) :: noxl
	real*8  :: xv(10+2*n_melem), rv(10+2*n_melem), dnh, dnhe
	integer :: j, nx, itr

	res = 0.0d0;  res_triplet = 0.0d0;  ok = .false.;  why = ''
	if (.not. ieq_rates_ready) then
		why = 'no equilibrium sweep has left a cell state to measure'
		return
	endif

	! The element totals of the state handed in, counted the way the sweep
	! counts them (one definition, utilities).
	nhil   = f_sp(:,1)*rho*n0
	nhiil  = f_sp(:,2)*rho*n0
	nheil  = f_sp(:,3)*rho*n0
	nheiil = f_sp(:,4)*rho*n0
	nheiiil= f_sp(:,5)*rho*n0
	nmoll  = 0.0d0
	noxl   = 0.0d0
	if (thereis_mol) then
		nmoll(:,1) = f_sp(:,isp_H2)  *rho*n0
		nmoll(:,2) = f_sp(:,isp_H2p) *rho*n0
		nmoll(:,3) = f_sp(:,isp_H3p) *rho*n0
		nmoll(:,4) = f_sp(:,isp_HeHp)*rho*n0
	endif
	if (thereis_oxychem) then
		noxl(:,1) = f_sp(:,isp_OH) *rho*n0
		noxl(:,2) = f_sp(:,isp_H2O)*rho*n0
		noxl(:,3) = f_sp(:,isp_CO) *rho*n0
	endif
	call hydrogen_helium_nuclei_density(nhil,nhiil,nheil,nheiil,nheiiil, &
	                                    nHl,nHel,nmoll,noxl)
	dnh  = 0.0d0
	dnhe = 0.0d0
	do j = 1, N
		dnh  = max(dnh,  abs(nHl(j)  - ieq_rate_cell(j)%nh)               &
		                 /max(ieq_rate_cell(j)%nh, 1.0d-99))
		if (thereis_He) dnhe = max(dnhe,                                  &
		                 abs(nHel(j) - ieq_rate_cell(j)%nhe)              &
		                 /max(ieq_rate_cell(j)%nhe, 1.0d-99))
	enddo
	if (max(dnh, dnhe) .gt. 1.0d-8) then
		why = 'the element totals have moved since the sweep that '//     &
		      'built these rates'
		return
	endif

	itr = ieq_triplet_row()
	call save_ieq_closure_scratch(ws)
	do j = 1, N
		call ieq_state_vector(j, rho, f_sp, xv, nx)
		call ieq_install_cell(j)
		res(j) = normalized_reaction_residual(xv(1:nx), nx,               &
		                    ieq_mbase_stored, ieq_ne_cell(j), rv(1:nx))
		if (itr .gt. 0 .and. itr .le. nx) res_triplet(j) = rv(itr)
	enddo
	call restore_ieq_closure_scratch(ws)
	if (.not. ieq_closure_scratch_matches(ws))                           &
		write(*,'(A)') ' (closure residual) WARNING: the cell scratch '// &
		     'was NOT reinstated after the isolated evaluation'
	ok = .true.
	end subroutine ionization_closure_residual_profile

	!----------------------------------!

	subroutine save_ieq_rate_state(s)
	! Hold aside the rate state of every cell, allocation status included.
	type(ieq_rate_state), intent(out) :: s
	if (allocated(ieq_rate_cell))      s%rate_cell      = ieq_rate_cell
	if (allocated(ieq_ne_cell))        s%ne_cell        = ieq_ne_cell
	if (allocated(ieq_TK_cell))        s%TK_cell        = ieq_TK_cell
	if (allocated(ieq_ntot_cell))      s%ntot_cell      = ieq_ntot_cell
	if (allocated(ieq_met_coef))       s%met_coef       = ieq_met_coef
	if (allocated(ieq_nonroot_streak)) s%nonroot_streak =                 &
	                                              ieq_nonroot_streak
	s%rates_ready  = ieq_rates_ready
	s%neq_stored   = ieq_neq_stored
	s%mbase_stored = ieq_mbase_stored
	s%iox_stored   = ieq_iox_stored
	end subroutine save_ieq_rate_state

	!----------------------------------!

	subroutine restore_ieq_rate_state(s)
	! Put back what save_ieq_rate_state held: the module comes back to the
	! state it was in and not merely to the same numbers.
	type(ieq_rate_state), intent(in) :: s
	if (allocated(s%rate_cell)) then
		ieq_rate_cell = s%rate_cell
	else if (allocated(ieq_rate_cell)) then
		deallocate(ieq_rate_cell)
	endif
	call put_back_ieq_1d(s%ne_cell,   ieq_ne_cell)
	call put_back_ieq_1d(s%TK_cell,   ieq_TK_cell)
	call put_back_ieq_1d(s%ntot_cell, ieq_ntot_cell)
	if (allocated(s%met_coef)) then
		ieq_met_coef = s%met_coef
	else if (allocated(ieq_met_coef)) then
		deallocate(ieq_met_coef)
	endif
	if (allocated(s%nonroot_streak)) then
		ieq_nonroot_streak = s%nonroot_streak
	else if (allocated(ieq_nonroot_streak)) then
		deallocate(ieq_nonroot_streak)
	endif
	ieq_rates_ready  = s%rates_ready
	ieq_neq_stored   = s%neq_stored
	ieq_mbase_stored = s%mbase_stored
	ieq_iox_stored   = s%iox_stored
	end subroutine restore_ieq_rate_state

	!----------------------------------!

	subroutine put_back_ieq_1d(kept, live)
	real*8, allocatable, intent(in)    :: kept(:)
	real*8, allocatable, intent(inout) :: live(:)
	if (allocated(kept)) then
		live = kept
	else if (allocated(live)) then
		deallocate(live)
	endif
	end subroutine put_back_ieq_1d

	!----------------------------------!

	end module ionization_equilibrium

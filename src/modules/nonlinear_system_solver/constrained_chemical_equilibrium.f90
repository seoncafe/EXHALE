	module constrained_chemical_equilibrium
	! Chemical equilibrium of the molecular H/He network -- H0, H+, H2, H2+,
	! H3+, HeH+, He(1^1S), He+, He++, the He 2^3S metastable where it is
	! tracked, the OH/H2O oxygen carriers where the oxygen chemistry is on,
	! and the trace-metal ionization stages -- solved in SPECIES NUMBER
	! DENSITIES with LOGARITHMIC unknowns, with element conservation written
	! as explicit residual rows, the electron density taken from charge
	! neutrality, and the whole system embedded in a continuation in the
	! radiation field from the dense molecular limit up to the full local
	! field (docs/supersonic_molecular_base.md section 11.5-B).
	!
	! WHY THE FORMULATION IS THIS ONE. The fraction systems
	! (System_HeH_mol / System_HeH_mol_metals) solve for stage fractions of
	! each element and close the neutral stage as the remainder, so nothing
	! stops an iterate from leaving the physical composition domain: a
	! negative free H0 or free He, or an element with more nuclei in its
	! tracked stages than it has. When that happens the accepted state has to
	! be projected back onto the element budget afterwards, and a projection
	! is not a solve -- it moves the state off the reaction balances it was
	! supposed to satisfy. The formulation here removes the possibility
	! instead of repairing it:
	!
	!   * the unknowns are u_k = ln(n_k), so n_k = exp(u_k) > 0 for EVERY
	!     iterate, converged or not. Logarithmic composition variables for
	!     exactly this reason are the CEA method of Gordon & McBride (1994),
	!     NASA RP-1311, which iterates on the logarithms of the compositions
	!     precisely so that no iterate can go negative;
	!   * element conservation is not a closure but a RESIDUAL ROW, one per
	!     element, driven to zero along with the reaction rows. HeH+ carries
	!     one H nucleus and one He nucleus, and the oxygen carriers OH and
	!     H2O carry H nuclei as well as O nuclei; each of them therefore
	!     appears in TWO conservation rows at once, and both are satisfied at
	!     the accepted root. That simultaneous satisfaction is the section
	!     11.5-B requirement, and it is what post-solve rescaling of stage
	!     fractions cannot do -- rescaling one element's stages onto its
	!     budget moves a shared species and breaks the other element's;
	!   * the electron density is never an unknown. It is the charge-
	!     neutrality sum of the ion densities, so neutrality holds identically
	!     at every iterate rather than approximately at the root.
	!
	! THE CONTINUATION. The molecular network is bistable in its starting
	! point because the cell is: a molecular basin in the dense, shielded base
	! and an atomic basin above the H2 -> H front. Three unconnected starting
	! points either land in the right basin or do not. A natural-parameter
	! continuation replaces the guess by a PATH: the radiation field is scaled
	! by lambda and the solution is tracked from a WEAK field, six decades
	! below the cell's own, at which the composition is known in closed form
	! (the H2 dissociation equilibrium of the local p, T with every element in
	! its own ionization balance at that field), up to lambda = 1, the cell's
	! actual field. Each rung starts from the previous rung's converged root,
	! so the solver is never asked to cross a basin boundary in one step; a
	! rung that fails is bisected geometrically and retried from the last
	! converged root. The method is the natural-parameter continuation of
	! Allgower & Georg, Introduction to Numerical Continuation Methods (SIAM
	! Classics, 2003), in its simplest form -- no tangent predictor, because
	! the cost of one extra rung is far below the cost of a Jacobian. Why the
	! path stops short of lambda = 0, which is a singular point of it rather
	! than an end of it, is recorded at the ladder below.
	!
	! Each rung is solved by the in-tree MINPACK hybrd (Powell's hybrid
	! trust-region method with a forward-difference Jacobian; More, Garbow &
	! Hillstrom 1980) -- the general driver, not its hybrd1 wrapper, and for
	! one reason. hybrd1 fixes the trust-region parameter at
	! factor = 100, which makes the INITIAL trust radius 100*||x||. For the
	! fraction systems, whose unknowns are O(1), that is the intended
	! generous first step. Here the unknowns are logarithms of densities,
	! ||u|| ~ 130, so the same rule opens the first trust region at ~10^4 in
	! ln-space -- a first step of e^10000 in every density. Measured, the
	! solve then made NO progress at all (info = 4 at every rung, with the
	! residual left exactly proportional to lambda, i.e. equal to the
	! radiation term it never moved to accommodate). The unknowns are
	! therefore taken RELATIVE to the previous rung's composition,
	! u_k = ln(n_k/n_k,ref), and the trust-region parameter is divided by the
	! norm of the starting point so that MINPACK's delta = factor*||x|| opens
	! at trust_radius_ln, an initial trust radius set directly in ln-space.
	! One is the natural value: a first step
	! of at most a factor of e in any density, which is the size of a
	! continuation rung. Positivity is untouched -- n_k = n_k,ref exp(u_k)
	! with n_k,ref > 0 -- and so is the variable scaling, unity in every
	! direction (mode = 2, diag = 1), which is right here because every
	! unknown is a logarithm and they all share one natural step.
	!
	! The algorithm inside hybrd is MINPACK's own. The single thing this
	! module selects in it is the forward-difference step rule of fdjac1: the
	! hybrd call below passes unit_step_floor = .true., which floors the step
	! at unit magnitude instead of letting it shrink with the variable. Every
	! other caller of hybrd reaches it through hybrd1, which passes .false.
	! and so keeps MINPACK's own rule exactly. See fd_eta.
	!
	! WHY THE TRACE SPECIES ARE NOT ALL ITERATED ON, MEASURED. With every
	! species in the Newton system the first rung converged to 1e-14 but
	! ADVANCING the field stalled, for a reason that is a property of the
	! formulation rather than of the step size (initial trust radii of 1 and
	! of 10 were both measured, and both failed the same way). A trace
	! intermediate that the network holds far below its element -- H2+, H3+,
	! HeH+, He++, the upper metal stages in a shielded base -- has a balance
	! row whose production and loss are both far below the row's turnover
	! scale, so the scaled row is zero to round-off and carries no
	! information about that species. Its direction is null: the measured
	! iterates moved such unknowns by up to e^3700 without changing the
	! residual, MINPACK grew its trust region because nothing got worse, and
	! the state that came back could not close the rows that do matter (the
	! Ca I <-> Ca II balance, in the measured case). Those species are
	! therefore held out of the iteration and placed afterwards from their
	! own balance rows -- the trace-species treatment of the CEA method, see
	! eps_hold below.
	!
	! WHAT IS AND IS NOT GUARANTEED. Positive species densities are
	! structural: they hold for every iterate at every lambda. Element
	! conservation and charge neutrality hold AT THE ACCEPTED ROOT --
	! neutrality identically, conservation to the internal tolerance the
	! conservation rows are driven to. Nothing here guarantees that the root
	! reached at lambda = 1 is the root the cell wants, which is why THIS
	! MODULE DOES NOT ACCEPT ANYTHING. It returns a candidate in the fraction
	! layout of the cell's system and says whether the continuation reached
	! the full field; the caller judges it with the same judge every other
	! candidate faces (ionization_fractions_physical and
	! normalized_reaction_residual against ieq_res_tol, section 11.5-A). A
	! successful library status and a small step in the transformed variables
	! are evidence of neither.
	!
	! The rows themselves are NOT re-derived here. The reaction balances are
	! assembled by calling the same routines the fraction systems call --
	! mol_heh_rows, oxygen_carrier_rows, metal_rows, cx_add_to_fvec,
	! he_h_cx_fvec -- with this module's densities and with the photo rates
	! multiplied by lambda, and each row is divided by the same turnover
	! scale (mol_inv_turnover) the fraction systems use. There is one
	! definition of the network in the code, and this module is a different
	! set of unknowns over it, not a second copy of it.

	use global_parameters, only: thereis_HeITR, thereis_metals,           &
	                             thereis_oxychem, HeH
	use species_table,     only: n_melem, iel_O
	use ion_cell_state,    only: ieq_cell
	use ion_residual_core, only: metal_rows, metal_electron_sum
	use System_HeH_mol,    only: mol_heh_rows, oxygen_carrier_rows,        &
	                             mol_inv_turnover,                         &
	                             set_mol_turnover_rates,                   &
	                             set_oxygen_turnover_rates,                &
	                             oj3, oj4, oj5, oj7,                       &
	                             ! the diagnostic row decomposition reads the
	                             ! molecular coefficients from the module
	                             ! that owns them, and rebuilds them through
	                             ! the routine that owns that
	                             set_mol_coeffs,                           &
	                             mk9, mk10, mk13, mk17, mk20, mk23
	use Cooling_Coefficients, only: f_penning_HeI23S
	use System_HeH_mol_metals, only: set_mol_metal_turnover_rates
	use System_HeH_metals, only: met_nelem, met_ntot, met_g0, met_g1,      &
	                             met_b0, met_b1, met_a1, met_a2, met_top,  &
	                             set_metal_coeffs
	use charge_exchange,   only: cx_add_to_fvec, he_h_cx_fvec,             &
	                             cx_metal_base, cx_add_to_turnover,        &
	                             cx_set_cell
	use lower_column,      only: q_h2_equilibrium
	use oxygen_rates,      only: oxygen_chemical_equilibrium_fractions, &
	                             rk_D1_Hep_CO

	implicit none
	private
	public :: equilibrium_from_molecular_limit
	! Opt-in diagnostic (docs/charge_exchange_cancellation_limit.md section 5).
	! cce_probe_from_dump is the entry the standalone driver calls; it is not
	! reached by any production path.
	public :: cce_probe_from_dump
	! The layout of the cell's network, and the species slots and array
	! length that index it. Read by the structural assertions of
	! src/tests/constrained_network_layout/; no production path calls them.
	public :: constrained_network_layout_of_cell
	public :: is_HII, is_H2, n_species_max

	! ---------------------------------------------------------------
	! Species slots. Every species the network can carry has a fixed slot in
	! the density array, present or not; a species the cell does not carry
	! keeps density zero and contributes zero to every row it appears in, so
	! no branch is needed inside the residual.
	integer, parameter :: is_HI     = 1
	integer, parameter :: is_HII    = 2
	integer, parameter :: is_H2     = 3
	integer, parameter :: is_H2p    = 4
	integer, parameter :: is_H3p    = 5
	integer, parameter :: is_HeHp   = 6
	integer, parameter :: is_HeI_SI = 7      ! free neutral He, ground singlet
	integer, parameter :: is_HeII   = 8
	integer, parameter :: is_HeIII  = 9
	integer, parameter :: is_HeITR  = 10     ! He 2^3S metastable
	integer, parameter :: is_OH     = 11
	integer, parameter :: is_H2O    = 12
	! Metal element e occupies is_met0 + 3*(e-1) + (0,1,2) = X0, X+, X++.
	integer, parameter :: is_met0   = 13
	integer, parameter :: n_species_max = 12 + 3*n_melem

	! Largest fraction-layout system this module is ever handed: the
	! molecular block (7) + the He 2^3S metastable (1) + the two oxygen
	! carriers + two rows per metal element. Same bound as n_x_max in
	! ioniz_eq, which builds the N_eq this module is called with.
	integer, parameter :: n_fraction_rows_max = 10 + 2*n_melem

	! Length of MINPACK's packed upper-triangular factor for the largest
	! system (hybrd's lr = n(n+1)/2).
	integer, parameter :: lr_max = (n_species_max*(n_species_max+1))/2

	! Initial trust radius of every rung, in ln-space: a first step of at most
	! a factor of e in any one density. MINPACK opens the region at
	! delta = factor*||diag*x||, so the caller divides by the norm of the
	! starting point to leave the radius at this value. See the module header.
	real*8,  parameter :: trust_radius_ln = 1.0d0

	! Offset carried by the unknowns: v_j = u_shift + ln(n_j/n_ref_j), so a
	! rung starts at v = u_shift rather than at zero.
	!
	! IT NO LONGER CARRIES THE FINITE-DIFFERENCE STEP. The step is now
	! floored inside fdjac1 itself, h_j = eps*max(1,|v_j|), selected by the
	! unit_step_floor argument of the hybrd call below. An offset can only
	! move the point at which MINPACK's own h = eps*|v_j| collapses; it
	! cannot remove it, and after the previous offset-only fix the collapse
	! merely sat at u = -1, i.e. at n = 0.368 n_ref, an ordinary value for a
	! rung to pass through. The floor removes it everywhere.
	!
	! WHAT THE OFFSET DOES CARRY, MEASURED. MINPACK's termination tests are
	! RELATIVE to the norm of the scaled unknowns: hybrd stops on
	! delta <= xtol*xnorm with xnorm = ||diag*x||, and reports info = 3 on
	! p1*max(p1*delta,pnorm) <= epsmch*xnorm. The unknowns here are
	! DEPARTURES, so without an offset they start at zero and a converged
	! rung leaves them near zero: xnorm collapses with them, the relative
	! test cannot fire, and the solve runs on to a slow-progress or
	! evaluation-count exit instead of stopping at its own root. The offset
	! gives the same composition a norm of order sqrt(n_unknown) and puts the
	! test back in reach.
	!
	! Measured on the He/H = 1 hot-Uranus case at 12000 steps, with the floor
	! in place in BOTH: u_shift = 1 left no cell on the 11.4 K temperature
	! floor and none below 50 K, and promoted 481 cells to this solve
	! (3208 field solves); u_shift = 0 pinned 22 cells at the floor and 46
	! below 50 K, and promoted 675 (4331 field solves). The difference is not
	! the step size. Sweeping fd_eta over 1e-4, 1e-5 and 1e-6 -- two decades
	! of h -- left both counts at zero with the offset in place, while
	! dropping the offset changes h by at most a factor of two anywhere the
	! iterates reach.
	real*8,  parameter :: u_shift = 1.0d0

	! Finite-difference size for the constrained solve, chosen by the
	! measured sweep of docs/charge_exchange_cancellation_limit.md phase B
	! rather than left at MINPACK's default. fdjac1 forms eps =
	! sqrt(max(epsfcn, epsmch)) and, under unit_step_floor, uses
	! h_j = eps*max(1,|v_j|), so passing epsfcn = fd_eta^2 sets the step to
	! fd_eta*max(1,|v_j|) -- the absolute-floor rule the review asks for,
	! with no value of the iterate at which it degenerates. Measured on the
	! cell-387 state, against a
	! central-difference reference: MINPACK's own sqrt(epsmch) = 1.5e-8 got
	! the H3+ column wrong by 70% and the HeH+ column by 86%, and the
	! agreement improved monotonically with a larger step -- 8.6e-1 at
	! eta = 1e-7, 6.8e-2 at 1e-6, 1.1e-2 at 1e-5, 4.6e-3 at 1e-4. 1e-4 is
	! the measured optimum and 1e-5 is within a factor 2.5 of it, so the
	! choice sits on a plateau rather than at a point. Larger steps are not
	! explored because truncation error then grows as eta.
	!
	! With the floor in place the value stopped mattering to the solution, as
	! it should: on the He/H = 1 hot-Uranus case above, 1e-4, 1e-5 and 1e-6 all
	! left no cell on the temperature floor, all gave the same log10 Mdot to
	! the printed 8.89, and differed only in how many cells the sweep promoted
	! to this solve at all (481, 524, 483 of 502).
	real*8,  parameter :: fd_eta = 1.0d-4

	! Radiation-field continuation ladder, from a field six decades below the
	! cell's own up to the cell's own. The rungs are logarithmically spaced
	! because the composition responds to the FIRST photons far more strongly
	! than to the last ones: over the first decades a shielded molecular base
	! goes from a purely collisional ionization balance to a photoionized
	! one, while between 0.3 and 1 nothing qualitative happens.
	!
	! WHY THE LADDER DOES NOT START AT lambda = 0, WHICH WAS TRIED AND
	! MEASURED. In the shielded molecular base the trace ionization stages
	! are made by photons and removed by recombination, so their densities
	! are LINEAR in lambda: n_X+ ~ lambda g n_X0/(alpha n_e). Their
	! logarithms therefore run to -infinity as lambda -> 0, and the
	! radiation-free composition puts them not near the lambda > 0 branch but
	! at whatever floor the seed uses -- measured, 22 decades below the
	! lambda = 1e-6 solution. At that point the row those species balance is
	! numerically independent of them (their term is 22 decades below the
	! photo term), the finite-difference Jacobian column is zero to round-off,
	! and MINPACK correctly reports that no step improves anything: every
	! rung returned info = 5 with the residual left exactly proportional to
	! lambda. lambda = 0 is a singular point of the continuation, not a
	! starting point. The first rung is instead a WEAK but nonzero field, at
	! which the H2 partition is still the thermochemical one (six decades of
	! photodissociation below the cell's own changes nothing in a shielded
	! base) while the trace ionization is already on the branch the
	! continuation follows.
	integer, parameter :: n_rung = 6
	real*8,  parameter :: lambda_rung(n_rung) =                            &
		(/ 1.0d-6, 1.0d-4, 1.0d-2, 1.0d-1, 3.0d-1, 1.0d0 /)
	! Smallest field step the bisection will take before giving up: below
	! this the field is so far below the collisional rates that a failure is
	! not a continuation-step failure but a failure of the rung itself.
	real*8,  parameter :: lambda_floor = 1.0d-12
	! Total hybrd1 calls one cell may spend on the continuation.
	integer, parameter :: max_field_solves = 40
	! Internal acceptance of ONE rung: largest absolute scaled residual
	! (reaction rows against their turnover rate, conservation rows already
	! dimensionless) of the vector hybrd1 returns. The MINPACK exit code
	! takes no part -- info = 1 is an xtol statement about the step and
	! info = 4 routinely returns finished roots it cannot certify
	! (docs/Update_EXHALE_stage1.md section 113). This is the tolerance of an
	! INTERMEDIATE rung; the final candidate is judged by the caller.
	real*8,  parameter :: rung_res_tol = 1.0d-6

	! Acceptance of an INTERMEDIATE rung, one at a field below the cell's
	! own. Such a rung is not an answer, it is a place to continue from, and
	! holding it to the final tolerance is stricter than the method needs:
	! measured on the H2 front of the hot-Uranus gate, the continuation
	! advanced to lambda ~ 7e-6 and then spent its whole field budget
	! creeping, because the full residual sat at 1.0-1.6e-6 -- within a
	! factor of 1.6 of the final tolerance -- on one metal row, while
	! element conservation was already at 1e-11. Loosening the intermediate
	! rungs cannot loosen the answer: the rung at lambda = 1 is still held to
	! rung_res_tol over the full row set, and the caller then judges that
	! candidate again with ionization_fractions_physical and
	! normalized_reaction_residual against ieq_res_tol. A path that drifts
	! can therefore only produce a candidate the caller rejects, never one it
	! wrongly accepts.
	real*8,  parameter :: rung_res_tol_intermediate = 1.0d-4

	! Largest log-RATIO, in magnitude, that the exponential is evaluated at.
	! The bound is applied on BOTH sides (bounded_log_ratio): exp(250) times
	! the reference density sits far above any density this network reaches
	! (n_H <= ~1e18 cm^-3 at a microbar base), and exp(-250) times it is
	! still a normal positive number, so neither side binds anywhere near a
	! root. They only keep a trust-region step that has wandered decades
	! outside the domain from overflowing to infinity or underflowing to
	! exactly zero -- either of which would destroy the information the
	! residual carries and, in the zero case, the module's own positivity
	! claim (docs/charge_exchange_cancellation_limit.md section 3.4).
	real*8,  parameter :: ln_ratio_max = 2.5d2

	! Seed floor of a species whose closed-form molecular limit is zero,
	! as a fraction of its element's nuclei. Logarithmic unknowns need a
	! strictly positive starting point -- ln(0) does not exist -- and 1e-25
	! of the element is 15 decades below any density the network resolves,
	! so it is a starting point and not a physical statement.
	real*8,  parameter :: seed_floor_fraction = 1.0d-25

	! Sweeps used to bring the starting point's electron density to its own
	! charge-neutrality fixed point (molecular_limit_composition). The map is
	! damped by a geometric mean, so it converges linearly with ratio 1/2 in
	! the logarithm: eight sweeps close the measured factor of 27 to better
	! than a per cent, which is far inside what a starting point needs.
	integer, parameter :: n_neutrality_sweep = 8

	! ---------------------------------------------------------------
	! TRACE-SPECIES HOLDOUT (the CEA treatment).
	!
	! A species the network holds far below its element has a balance row
	! whose production and loss are both far below the row's turnover scale.
	! The scaled row is then zero to round-off: it carries no information
	! about that species, its Jacobian column vanishes, and the trust-region
	! step is free to move it anywhere. Measured on the hot-Uranus base, such
	! unknowns moved by up to e^3700 with no change in the residual, and the
	! state that came back could not close the rows that do matter.
	!
	! The remedy is the one Gordon & McBride (1994), NASA RP-1311 use for the
	! same reason in the CEA method: species below a size threshold are
	! OMITTED from the iteration equations and restored afterwards from the
	! equilibrium relations. Here the setting is kinetic rather than
	! free-energy minimization, so the relation a held species is restored
	! from is its own balance row, but the treatment and its justification
	! are the same -- a species that carries none of the element budget and
	! none of the charge cannot influence the rest of the system, so solving
	! for it and solving for everything else are separable.
	!
	! A species is held when its density falls below this fraction of its
	! element's nuclei. 1e-10 is where the holdout has to sit: a species at
	! 1e-10 of its element changes the element's conservation row by 1e-10
	! and the electron density by at most that, both four decades below the
	! rung tolerance, so holding it cannot move the solution; while a decade
	! or two higher would start holding species that genuinely carry the
	! chemistry of the shielded base (H3+ runs at ~1e-8 of the H nuclei
	! there). The value is the one the measurement below was taken at.
	real*8,  parameter :: eps_hold = 1.0d-10

	! A held species whose placed density comes back above this multiple of
	! the holdout threshold was held wrongly: it is re-admitted as an unknown
	! and the rung is solved again. The factor keeps a species that lands
	! just at the threshold from oscillating in and out between cycles.
	real*8,  parameter :: readmit_factor = 1.0d1

	! Fixed-point sweeps over the held set when placing it. The held species
	! feed each other (H2+ makes H3+, HeH+ makes H3+), so one pass places
	! each against stale partners; five is well past the point where the
	! placement stops moving for the couplings this network has.
	integer, parameter :: n_placement_sweep = 5

	! Times one rung may re-partition and re-solve after a re-admission.
	integer, parameter :: max_partition_cycle = 3

	! ---------------------------------------------------------------
	! Cell state of the solve in progress. The equilibrium sweep is
	! OpenMP-parallel over cells and the residual is reached through
	! MINPACK, which carries no user context, so this is threadprivate
	! exactly as ieq_cell and mol_inv_turnover are.
	!
	! The layout splits in two. The CELL layout -- which species the cell
	! carries, which balance row each one owns, how many conservation rows
	! there are -- is fixed for the cell. The RUNG partition -- which of
	! those species are held out of the iteration as trace, and hence which
	! unknowns and rows the Newton system actually has -- is rebuilt at every
	! rung and every re-admission cycle.
	integer, save :: n_fraction_row         ! N_eq of the cell's fraction system
	integer, save :: n_conservation_row     ! H, He and each element carried
	integer, save :: metal_base, oxygen_base
	real*8,  save :: field_scale            ! lambda of the rung being solved
	logical, save :: h2_is_fixed, oxygen_carriers_are_fixed, hp_is_fixed
	real*8,  save :: nuclei_H, nuclei_He, nuclei_O
	logical, save :: element_carried(n_melem)
	! Does the cell carry this species at all (an absent element, the pinned
	! X++ of a two-stage element and the oxygen carriers of a run without the
	! oxygen chemistry do not exist).
	logical, save :: species_exists(n_species_max)
	! Is this species' partition owned by the carrier transport rather than
	! by this solve: never an unknown, and its balance row never a row.
	logical, save :: species_fixed(n_species_max)
	! The balance row this species owns in the fraction layout, or 0. HI,
	! He(1^1S) and each metal's neutral stage own NO row -- they are what the
	! element conservation rows determine, which is exactly why the system is
	! square: those species number 2 + (elements carried), the conservation
	! rows number the same.
	integer, save :: row_of_species(n_species_max)
	! Held out of the iteration as trace for the rung in progress.
	logical, save :: species_held(n_species_max)
	integer, save :: n_unknown              ! size of the log-density system
	integer, save :: n_reaction_row         ! balance rows in the iteration
	integer, save :: species_of_unknown(n_species_max)
	integer, save :: fraction_row_of_reaction(n_fraction_rows_max)
	! Composition the current rung measures its unknowns against: the
	! previous rung's converged densities (the seed at the first rung), so
	! that u = ln(n/n_ref) and every rung starts at u = 0, i.e. at
	! v = u_shift in the shifted variable MINPACK sees. Strictly positive
	! for every unknown species, which is what keeps n = n_ref exp(u)
	! positive.
	real*8,  save :: reference_density(n_species_max)
	! Reciprocal turnover rate each reaction row is divided by, in the
	! ACCEPTANCE normalization: the solver scale mol_inv_turnover widened by
	! the metal <-> H/He charge-exchange bound that scale leaves out of its
	! metal rows. See set_turnover_scales_at_field for why the widening is
	! not optional here.
	real*8,  save :: row_scale(n_fraction_rows_max)
	! ---------------------------------------------------------------
	! Opt-in state dump of a continuation that failed to reach the full
	! field (docs/charge_exchange_cancellation_limit.md section 5, phase A).
	! Off unless the environment variable EXHALE_CCE_DUMP names a file, and
	! then written ONCE, for the first failing cell of the run. These three
	! are deliberately NOT threadprivate: the once-only guard has to be
	! shared, and it is taken inside a critical region.
	character(len=512), save :: cce_dump_path  = ' '
	logical, save :: cce_dump_written = .false.
	logical, save :: cce_ok_written   = .false.
	! Seed perturbation, for the probe's reproducibility test only: when the
	! flag is false -- which is every production path, and the initialized
	! value -- the seed is untouched and the array is never read. Shared,
	! not threadprivate, and written only by the single-threaded probe.
	logical, save :: seed_perturb_on = .false.
	real*8,  save :: seed_perturb(n_species_max)
	!$omp threadprivate(n_unknown, n_reaction_row, n_fraction_row,         &
	!$omp               n_conservation_row,                                &
	!$omp               metal_base, oxygen_base, field_scale,              &
	!$omp               h2_is_fixed, oxygen_carriers_are_fixed,            &
	!$omp               hp_is_fixed,                                       &
	!$omp               nuclei_H, nuclei_He, nuclei_O,                     &
	!$omp               species_exists, species_fixed, row_of_species,     &
	!$omp               species_held,                                      &
	!$omp               species_of_unknown, fraction_row_of_reaction,      &
	!$omp               element_carried, reference_density, row_scale)

	contains

	!----------------------------------!

	subroutine equilibrium_from_molecular_limit(x, nx, mbase, iox, n_e_ref, &
	                                            p_bar, reached_full_field,  &
	                                            n_field_solves)
	! Track the chemical equilibrium of this cell from the dense molecular
	! (radiation-free) limit up to its actual radiation field, in positive
	! species densities with element conservation and charge neutrality
	! enforced throughout, and return the composition reached at the full
	! field in the fraction layout of the cell's own system.
	!
	! x(nx)  inout. Ignored on entry. Written only when the continuation
	!        reaches lambda = 1, and then it holds the candidate in the
	!        System_HeH_mol / System_HeH_mol_metals layout: x(1) = n_HII/n_H,
	!        x(2) = n_HeII/n_He, x(3) = n_HeIII/n_He, x(4) = 2 n_H2/n_H,
	!        x(5) = 2 n_H2+/n_H, x(6) = 3 n_H3+/n_H, x(7) = n_HeH+/n_H,
	!        x(8) = n_HeITR/n_He, the oxygen carriers over the free oxygen
	!        family at iox, and the metal stages over their element totals
	!        from mbase upwards. Layout slots the cell does not carry (an
	!        absent element, the pinned X++ of a two-stage element) are 0.
	! mbase, iox  first metal row and first oxygen-carrier row of that layout.
	! n_e_ref     the cell's incoming electron density, the same value the
	!             driver built the turnover scales with.
	! p_bar       the cell's gas pressure [bar], for the H2 dissociation
	!             equilibrium of the lambda = 0 seed.
	! reached_full_field  did the continuation get to lambda = 1.
	! n_field_solves      hybrd1 calls consumed, for the cost ledger.
	!
	! The temperature is ieq_cell%T_K; nothing is passed twice.
	!
	! THE RESULT IS A CANDIDATE, NOT AN ACCEPTANCE. Reaching lambda = 1 says
	! the continuation completed, not that the composition solves the cell's
	! network: the caller applies ionization_fractions_physical and
	! normalized_reaction_residual against ieq_res_tol, the same judge every
	! other candidate faces (docs/supersonic_molecular_base.md section
	! 11.5-A/B).

	integer, intent(in)    :: nx, mbase, iox
	real*8,  intent(in)    :: n_e_ref, p_bar
	real*8,  intent(inout) :: x(nx)
	logical, intent(out)   :: reached_full_field
	integer, intent(out)   :: n_field_solves

	real*8  :: sden(n_species_max), sden_ref(n_species_max)
	real*8  :: sden_work(n_species_max)
	real*8  :: u(n_species_max), fres(n_species_max)
	real*8  :: par(60)
	! MINPACK hybrd work space (the arrays hybrd1 carves out of its wa).
	real*8  :: diag_u(n_species_max)
	real*8  :: fjac(n_species_max,n_species_max)
	real*8  :: rtri(lr_max), qtf(n_species_max)
	real*8  :: wa1(n_species_max), wa2(n_species_max)
	real*8  :: wa3(n_species_max), wa4(n_species_max)
	real*8  :: xtol, epsfcn, dpmpar
	integer :: maxfev, nfev, ml, mu, mode, nprint, lr_loc
	real*8  :: lam_try, lam_goal, lam_ref, lam_next, resmax, rung_tol
	real*8  :: trust_factor
	integer :: info_loc, k, isp, itarget, icycle
	logical :: have_ref, seeded, readmit

	reached_full_field = .false.
	n_field_solves     = 0

	! No hydrogen or no helium is not a cell this network describes: every
	! conservation row of the formulation divides by an element total.
	if (ieq_cell%nh .le. 0.0d0 .or. ieq_cell%nhe .le. 0.0d0) return

	call set_molecular_network_layout(nx, mbase, iox)

	! NO RUNG EXISTS YET FOR THIS CELL. Everything below the cell layout --
	! the holdout, the unknown and row lists, the reference composition the
	! logarithmic unknowns are measured against, the row scales and the
	! field of the rung -- is built per rung by the loop that follows. Until
	! it is, it must not still describe the previous cell this thread
	! solved: a cell's equilibrium is a function of that cell's state, not
	! of the order the sweep handed cells to threads
	! (docs/Update_EXHALE_stage1.md section 121).
	species_held(:)             = .false.
	n_unknown                   = 0
	n_reaction_row              = 0
	species_of_unknown(:)       = 0
	fraction_row_of_reaction(:) = 0
	reference_density(:)        = 0.0d0
	row_scale(:)                = 0.0d0
	field_scale                 = 0.0d0

	par(:)     = 0.0d0         ! MINPACK transport argument only, unread
	xtol       = sqrt(dpmpar(1))
	! fdjac1 forms eps = sqrt(max(epsfcn,epsmch)), so with the floored rule
	! this sets the difference step to fd_eta*max(1,|v|) (see fd_eta above).
	epsfcn     = fd_eta*fd_eta
	mode       = 2             ! caller-supplied variable scaling
	diag_u(:)  = 1.0d0         ! unity: every unknown is a logarithm
	nprint     = 0

	! Natural-parameter continuation. lam_ref is the largest field strength
	! whose root is in hand, lam_goal the ladder rung being worked toward,
	! and lam_try the field the next solve is aimed at. A failed rung is
	! bisected geometrically toward lam_ref and retried FROM lam_ref's root,
	! never from the failed iterate; a successful BISECTED rung then aims at
	! lam_goal again from the closer root rather than at the next ladder
	! value, which would jump straight back over the field level that just
	! failed and loop.
	have_ref = .false.
	lam_ref  = 0.0d0
	itarget  = 1
	lam_goal = lambda_rung(1)
	lam_try  = lam_goal

	do
		if (n_field_solves .ge. max_field_solves) exit

		! Until the first root is in hand there is nothing to continue
		! from, so the closed-form molecular composition is rebuilt at
		! whatever weak field is currently being attempted. Once a root
		! exists the continuation carries it forward instead.
		if (.not. have_ref) then
			call molecular_limit_composition(sden, p_bar, n_e_ref,     &
			                                 lam_try, seeded)
			if (.not. seeded) exit
			sden_ref(:) = sden(:)
		endif

		field_scale = lam_try
		call set_turnover_scales_at_field(n_e_ref, lam_try)

		! The rung, with the trace species held out of the Newton system and
		! placed afterwards. A held species that comes back large was held
		! wrongly, so it is re-admitted and the rung is solved again from
		! the updated composition; the cycle count bounds that, the field
		! budget bounds it again.
		sden_work(:) = sden_ref(:)
		do icycle = 1,max_partition_cycle
			if (n_field_solves .ge. max_field_solves) exit

			call set_rung_partition(sden_work)

			! Each rung measures its unknowns against the composition in
			! hand, so it starts at a departure of zero -- at v = u_shift
			! in the variable MINPACK sees (see the module header).
			reference_density(:) = sden_work(:)
			u(1:n_unknown) = u_shift

			ml     = n_unknown - 1   ! dense Jacobian, as hybrd1 sets it
			mu     = n_unknown - 1
			maxfev = 200*(n_unknown + 1)
			lr_loc = (n_unknown*(n_unknown + 1))/2

			! MINPACK opens the trust region at factor*||diag*x||, and
			! the unknowns start at u_shift, so the factor is divided
			! by that norm to leave the initial radius at
			! trust_radius_ln -- the value the method was measured with.
			trust_factor = trust_radius_ln                             &
			             /(u_shift*sqrt(dble(n_unknown)))

			! The trailing .true. selects the floored difference step
			! h_j = eps*max(1,|v_j|) in fdjac1 (see fd_eta above).
			call hybrd(constrained_equilibrium_residual, n_unknown, u,  &
			           fres, xtol, maxfev, ml, mu, epsfcn, diag_u,      &
			           mode, trust_factor, nprint, info_loc, nfev,      &
			           fjac, n_species_max, rtri, lr_loc, qtf, wa1,     &
			           wa2, wa3, wa4, par, .true.)
			n_field_solves = n_field_solves + 1

			do k = 1,n_unknown
				isp = species_of_unknown(k)
				sden_work(isp) = reference_density(isp)                &
				               * exp(bounded_log_ratio(u(k)))
			enddo

			call place_held_species(sden_work)

			! Re-admission: a held species whose placed density is well
			! above the holdout threshold belongs in the iteration.
			readmit = .false.
			do isp = 1,n_species_max
				if (.not. species_held(isp)) cycle
				if (sden_work(isp) .gt. readmit_factor*eps_hold        &
				    *element_nuclei_of_species(isp)) readmit = .true.
			enddo
			if (.not. readmit) exit
		enddo

		! The rung is accepted on the FULL system -- the held rows at their
		! placed densities included -- not on the reduced one the Newton
		! step saw.
		resmax = full_network_residual(sden_work)

		! The rung at the cell's own field is the answer and is held to the
		! final tolerance; the rungs below it are only places to continue
		! from (see rung_res_tol_intermediate).
		if (lam_try .ge. 1.0d0) then
			rung_tol = rung_res_tol
		else
			rung_tol = rung_res_tol_intermediate
		endif

		if (resmax .le. rung_tol) then
			sden_ref(:) = sden_work(:)
			have_ref = .true.
			lam_ref  = lam_try
			if (lam_try .ge. 1.0d0) then
				reached_full_field = .true.
				exit
			endif
			if (lam_try .ge. lam_goal) then
				! The goal itself is reached; aim at the next rung.
				itarget  = itarget + 1
				lam_goal = lambda_rung(min(itarget, n_rung))
			endif
			lam_try = lam_goal
		else
			! The rung failed. With a root in hand the step is bisected
			! geometrically toward it; without one the whole path is moved
			! to a weaker field, where the closed-form start above is
			! rebuilt -- the composition of a more strongly shielded layer,
			! which is the one the continuation can begin from.
			if (lam_ref .le. 0.0d0) then
				lam_next = lam_try/1.0d1
			else
				lam_next = sqrt(lam_ref*lam_try)
			endif
			if (lam_next .le. lambda_floor) exit
			if (lam_next .ge. lam_try)      exit
			lam_try = lam_next
		endif
	enddo

	! RESTORE THE SOLVER-PATH TURNOVER SCALES, ON EVERY EXIT PATH. The
	! scales are threadprivate module state of System_HeH_mol shared with
	! the cell solve and with normalized_reaction_residual, and this routine
	! has been rewriting them at each rung's lambda. Leaving them at the
	! last rung's field would make the caller judge this candidate -- and
	! every later candidate of this cell -- against the wrong normalization,
	! which is the acceptance judge itself going wrong silently. Nothing
	! below this line may return before it.
	call set_turnover_scales_at_field(n_e_ref, 1.0d0)

	if (.not. reached_full_field) then
		! Opt-in: preserve the first such state so that it can be
		! re-evaluated without the hydrodynamics (phase A of the review
		! plan). Costs one environment-variable read per run when off.
		call write_continuation_dump('EXHALE_CCE_DUMP', cce_dump_written,&
		                             n_e_ref, p_bar, nx, mbase, iox,    &
		                             lam_ref, lam_try, sden_work,       &
		                             info_loc)
		return
	endif

	! Opt-in: preserve the first state that DID reach the full field, so
	! that the seed-robustness test has an actual root to reproduce.
	call write_continuation_dump('EXHALE_CCE_DUMP_OK', cce_ok_written,   &
	                             n_e_ref, p_bar, nx, mbase, iox,         &
	                             lam_ref, lam_try, sden_ref, info_loc)

	call fraction_layout_from_densities(sden_ref, x, nx)

	end subroutine equilibrium_from_molecular_limit

	!----------------------------------!

	subroutine set_turnover_scales_at_field(n_e_ref, photo_scale)
	! Turnover scale of every row of the cell's fraction-layout system at a
	! radiation field scaled by photo_scale, in the order the driver builds
	! them: set_mol_turnover_rates resets the whole array, the oxygen rows
	! are appended to it, and the metal rows after those. photo_scale = 1
	! reproduces the driver's own scales exactly (each rate enters as
	! 1.0d0*rate, which is exact in IEEE arithmetic, with the term order
	! unchanged).
	!
	! It then builds row_scale, the reciprocal turnover this module's own
	! residual divides by, which is NOT the same as mol_inv_turnover: it is
	! mol_inv_turnover widened by the metal <-> H/He charge-exchange bound
	! (cx_add_to_turnover), exactly as normalized_reaction_residual widens it
	! at acceptance time.
	!
	! The widening is not a refinement here, it is the difference between a
	! solvable and an unsolvable rung. set_mol_metal_turnover_rates leaves
	! charge exchange out of the metal rows deliberately -- for the fraction
	! systems it is one contribution among several to the solver's path. But
	! in a shielded molecular base the O I <-> O II row is DOMINATED by the
	! resonant O + H+ <-> O+ + H pair, so dividing that row by a scale built
	! from photoionization and recombination alone leaves a residual that no
	! composition can push below the rung tolerance: measured, the rung
	! stalled at 3e-2 on precisely that row, at every field down to 1e-12.
	! Widening it also makes this module drive down exactly the quantity the
	! caller's acceptance judge measures, which is the property that lets a
	! reached rung mean anything. A positive diagonal scaling has the same
	! zeros, so it changes the path and not the root.
	real*8, intent(in) :: n_e_ref, photo_scale
	real*8  :: cxb(n_fraction_rows_max), el_tot(12)
	integer :: i

	call set_mol_turnover_rates(n_e_ref, photo_scale)
	if (thereis_oxychem) call set_oxygen_turnover_rates(oxygen_base,       &
	                                                    photo_scale)
	if (thereis_metals)  call set_mol_metal_turnover_rates(n_e_ref,        &
	                                                       photo_scale)

	do i = 1,n_fraction_row
		row_scale(i) = mol_inv_turnover(i)
	enddo

	if (thereis_metals) then
		el_tot(:) = 0.0d0
		el_tot(1:met_nelem) = met_ntot(1:met_nelem)
		el_tot(11) = nuclei_H       ! cx_H
		el_tot(12) = nuclei_He      ! cx_He
		cxb(1:n_fraction_row) = 0.0d0
		cx_metal_base = metal_base
		call cx_add_to_turnover(cxb, el_tot)
		cx_metal_base = 4
		do i = 1,n_fraction_row
			row_scale(i) = mol_inv_turnover(i)                     &
			             /(1.0d0 + cxb(i)*mol_inv_turnover(i))
		enddo
	endif

	end subroutine set_turnover_scales_at_field

	!----------------------------------!

	subroutine set_molecular_network_layout(nx, mbase, iox)
	! Which species this cell carries, which balance row of the fraction
	! layout each of them owns, and how many element conservation rows there
	! are. Fixed for the cell; the rung partition below chooses the Newton
	! system out of it.
	!
	! A carrier whose partition is imposed is marked fixed rather than
	! absent: something outside this cell owns it and the system is handed
	! the answer. Such a species still appears in every reaction row and
	! every conservation row it belongs to, but it is never an unknown and
	! its own balance row is never a row -- one fewer unknown and one fewer
	! row, so the system stays square. H2, H+ and the oxygen carriers are
	! marked independently, because H2 can be imposed on its own (the
	! lower-boundary reservoir composition) while H+ is imposed by the
	! transported ionization state ("Ionization transport") and the oxygen
	! carriers only by carrier transport.

	integer, intent(in) :: nx, mbase, iox
	integer :: e, i

	n_fraction_row = nx
	metal_base     = mbase
	oxygen_base    = iox

	h2_is_fixed               = ieq_cell%x_h2_fixed
	oxygen_carriers_are_fixed = ieq_cell%x_ox_fixed
	hp_is_fixed               = ieq_cell%x_hp_fixed
	nuclei_H  = ieq_cell%nh
	nuclei_He = ieq_cell%nhe
	nuclei_O  = ieq_cell%n_ofam

	element_carried(:) = .false.
	if (thereis_metals) then
		do e = 1,met_nelem
			element_carried(e) = (met_ntot(e) .gt. 1.0d-30)
		enddo
	endif

	species_exists(:)  = .false.
	species_fixed(:)   = .false.
	row_of_species(:)  = 0

	! The H/He block. HI and He(1^1S) own no balance row: they are what the
	! hydrogen and helium conservation rows determine.
	species_exists(is_HI)     = .true.
	species_exists(is_HII)    = .true.
	species_exists(is_H2)     = .true.
	species_exists(is_H2p)    = .true.
	species_exists(is_H3p)    = .true.
	species_exists(is_HeHp)   = .true.
	species_exists(is_HeI_SI) = .true.
	species_exists(is_HeII)   = .true.
	species_exists(is_HeIII)  = .true.
	row_of_species(is_HII)    = 1
	row_of_species(is_HeII)   = 2
	row_of_species(is_HeIII)  = 3
	row_of_species(is_H2)     = 4
	row_of_species(is_H2p)    = 5
	row_of_species(is_H3p)    = 6
	row_of_species(is_HeHp)   = 7
	if (thereis_HeITR) then
		species_exists(is_HeITR) = .true.
		row_of_species(is_HeITR) = 8
	endif
	if (thereis_oxychem) then
		species_exists(is_OH)  = .true.
		species_exists(is_H2O) = .true.
		row_of_species(is_OH)  = oxygen_base
		row_of_species(is_H2O) = oxygen_base + 1
	endif
	if (h2_is_fixed) species_fixed(is_H2) = .true.
	! The transported ionization state owns the H+ partition of the cell:
	! the fraction system this continuation is judged against replaces the
	! H+ balance row by the constraint x(1) = x_hp_fix
	! (System_HeH_mol_metals.f90 282-290), so solving that row here as a
	! reaction row would return the LOCAL photoionization/recombination root
	! and be measured against a row it was never asked to satisfy.
	if (hp_is_fixed) species_fixed(is_HII) = .true.
	if (oxygen_carriers_are_fixed .and. thereis_oxychem) then
		species_fixed(is_OH)  = .true.
		species_fixed(is_H2O) = .true.
	endif

	! The metal block. Each carried element contributes its neutral stage
	! (no row -- the element conservation row determines it) and the ionized
	! stages it tracks.
	n_conservation_row = 2
	if (thereis_metals) then
		do e = 1,met_nelem
			if (.not. element_carried(e)) cycle
			n_conservation_row = n_conservation_row + 1
			i = is_met0 + 3*(e-1)
			species_exists(i)   = .true.
			species_exists(i+1) = .true.
			row_of_species(i+1) = metal_base + 2*(e-1)
			if (met_top(e) .ge. 2) then
				species_exists(i+2) = .true.
				row_of_species(i+2) = metal_base + 2*(e-1) + 1
			endif
		enddo
	endif

	end subroutine set_molecular_network_layout

	!----------------------------------!

	logical function seed_species_of_cell(isp)
	! Whether the SEED has to start this species positive: a species the
	! cell carries and whose partition this solve owns. It is the rung
	! partition's unknown set with the holdout test removed -- no rung
	! exists when the seed is built -- and it is a function of the cell's
	! own layout alone, which is what the seed has to be.
	integer, intent(in) :: isp
	seed_species_of_cell = species_exists(isp) .and. .not. species_fixed(isp)
	end function seed_species_of_cell

	!----------------------------------!

	subroutine set_rung_partition(sden)
	! Split the cell's species into the ones the Newton system iterates on
	! and the ones held out of it as trace, at the composition sden, and
	! build the unknown and row lists that follow.
	!
	! A species is held when it sits below eps_hold of its element AND owns
	! a balance row. The second condition is not a convenience: HI, He(1^1S)
	! and a metal's neutral stage own no row, so holding one of them would
	! remove an unknown without removing a row and the system would not be
	! square -- and there is nothing to restore them from afterwards either,
	! since what determines them is element conservation, which stays in the
	! system. In a shielded base they are also the largest species there
	! are, so the condition never binds on them in practice; it is stated
	! here so that it cannot start binding on a state nobody anticipated.
	!
	! Squareness is structural. Every existing species that is neither fixed
	! nor held contributes one unknown, and one row if it owns one; the ones
	! owning no row number 2 + (elements carried), which is exactly
	! n_conservation_row. The check below is therefore an assertion about
	! this routine, not about the cell, and a failure is a coding error.

	real*8, intent(in) :: sden(n_species_max)
	integer :: i

	species_held(:) = .false.
	do i = 1,n_species_max
		if (.not. species_exists(i)) cycle
		if (species_fixed(i))        cycle
		if (row_of_species(i) .le. 0) cycle
		if (sden(i) .lt. eps_hold*element_nuclei_of_species(i))        &
			species_held(i) = .true.
	enddo

	n_unknown      = 0
	n_reaction_row = 0
	do i = 1,n_species_max
		if (.not. species_exists(i)) cycle
		if (species_fixed(i))        cycle
		if (species_held(i))         cycle
		n_unknown = n_unknown + 1
		species_of_unknown(n_unknown) = i
		if (row_of_species(i) .gt. 0) then
			n_reaction_row = n_reaction_row + 1
			fraction_row_of_reaction(n_reaction_row) =             &
				row_of_species(i)
		endif
	enddo

	if (n_reaction_row + n_conservation_row .ne. n_unknown) then
		write(*,'(A,I0,A,I0,A,I0,A)')                                  &
			' (constrained_chemical_equilibrium) partition is not square:', &
			n_reaction_row, ' reaction + ', n_conservation_row,        &
			' conservation rows against ', n_unknown, ' unknowns'
		flush(6)
		error stop 'constrained_chemical_equilibrium: non-square system'
	endif

	end subroutine set_rung_partition

	!----------------------------------!

	subroutine constrained_network_layout_of_cell(nx, mbase, iox, sden,    &
	                              n_unknown_out, n_reaction_row_out,       &
	                              n_conservation_row_out, is_unknown)
	! The layout the continuation would give the cell now standing in
	! ieq_cell, at the composition sden: how many unknowns it has, how many
	! reaction and conservation rows, and which species slots are unknowns.
	!
	! It exists so that the two structural statements of the layout can be
	! asserted without relaxing a wind: that reaction rows plus conservation
	! rows equal unknowns (squareness), and that a species whose partition
	! is imposed from outside the cell -- H2 by the lower-boundary reservoir,
	! H+ by the transported ionization state, OH and H2O by carrier
	! transport -- is not among the unknowns. It writes only the module state
	! that set_molecular_network_layout and set_rung_partition write, which
	! every solve rebuilds for its own cell before using.

	integer, intent(in)  :: nx, mbase, iox
	real*8,  intent(in)  :: sden(n_species_max)
	integer, intent(out) :: n_unknown_out, n_reaction_row_out
	integer, intent(out) :: n_conservation_row_out
	logical, intent(out) :: is_unknown(n_species_max)
	integer :: i

	call set_molecular_network_layout(nx, mbase, iox)
	call set_rung_partition(sden)

	n_unknown_out          = n_unknown
	n_reaction_row_out     = n_reaction_row
	n_conservation_row_out = n_conservation_row
	is_unknown(:) = .false.
	do i = 1,n_unknown
		is_unknown(species_of_unknown(i)) = .true.
	enddo

	end subroutine constrained_network_layout_of_cell

	!----------------------------------!

	subroutine constrained_equilibrium_residual(nu, u, fres, iflag, params)
	! Residual of the constrained system at log-densities u, at the current
	! continuation field field_scale. The first n_reaction_row entries are
	! the selected balance rows of the cell's own fraction-layout system,
	! assembled by the SAME routines that system calls and divided by the
	! SAME turnover scales; the remaining entries are the dimensionless
	! element conservation statements. Charge neutrality is not a row: the
	! electron density is built from the ion densities, so it holds
	! identically.
	!
	! params is the MINPACK transport argument and is never read, matching
	! the fraction systems.

	integer :: nu, iflag
	real*8  :: u(nu), fres(nu)
	real*8  :: params(60)

	real*8  :: sden(n_species_max)
	real*8  :: fvec(n_fraction_rows_max)
	integer :: i, k, ir

	! Densities. Positive by construction -- the reference density is
	! positive and the exponential cannot be -- which is the whole point of
	! the logarithmic unknowns. Species that are not unknowns of this rung
	! (a transported carrier, a species held out as trace) keep the density
	! the reference composition carries, so they still enter every row.
	sden(:) = reference_density(:)
	do k = 1,nu
		i = species_of_unknown(k)
		sden(i) = reference_density(i)                                 &
		        * exp(bounded_log_ratio(u(k)))
	enddo

	call network_balance_rows(sden, fvec)

	do ir = 1,n_reaction_row
		fres(ir) = fvec(fraction_row_of_reaction(ir))
	enddo

	call element_conservation_rows(sden, fres(n_reaction_row+1))

	return
	end subroutine constrained_equilibrium_residual

	!----------------------------------!

	double precision function bounded_log_ratio(v) result(w)
	! The exponent actually used to form a density, n = n_ref exp(w), from
	! an unknown v = u_shift + ln(n/n_ref).
	!
	! BOUNDED ON BOTH SIDES, and the lower bound is the point. The module
	! claims that no iterate can produce a nonpositive density; with only
	! the upper side capped that claim was false, because a sufficiently
	! negative exponent underflows exp to exactly zero
	! (docs/charge_exchange_cancellation_limit.md section 3.4). Clamping the
	! low side restores the claim exactly: exp(-250) = 2.6e-109 times any
	! reference density this network carries is still a normal, strictly
	! positive number.
	!
	! Neither bound is reachable near a solution -- they sit 108 decades from
	! any density the chemistry resolves -- so neither changes a root. What
	! the lower bound does change is the failure mode of a trust-region step
	! that has wandered far outside the domain: instead of a zero density
	! with a structurally dead Jacobian column, it produces a tiny positive
	! one whose column is merely uninformative, and the solve can come back.
	real*8, intent(in) :: v
	w = min(max(v - u_shift, -ln_ratio_max), ln_ratio_max)
	end function bounded_log_ratio

	!----------------------------------!

	subroutine network_balance_rows(sden, fvec)
	! Every balance row of the cell's fraction-layout system, evaluated at
	! the species densities sden and at the continuation's current field, and
	! divided by the acceptance turnover scale. One definition, used by the
	! Newton residual, by the placement of the held species and by the rung
	! acceptance -- so that all three judge the same rows.
	!
	! The rows themselves are not written here: they are the same routines
	! the fraction systems call, in the same order, with the photo rates
	! multiplied by the continuation parameter.

	real*8, intent(in)  :: sden(n_species_max)
	real*8, intent(out) :: fvec(n_fraction_rows_max)

	real*8  :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
	real*8  :: gm0(n_melem), gm1(n_melem)
	real*8  :: xzero(n_fraction_rows_max)
	real*8  :: n_e, n_hei, lam
	integer :: i, e

	lam = field_scale

	! Metal stage densities in canonical element order; an element the cell
	! does not carry stays at zero, as its fraction-system rows do.
	nm0(:) = 0.0d0
	nm1(:) = 0.0d0
	nm2(:) = 0.0d0
	if (thereis_metals) then
		do e = 1,met_nelem
			if (.not. element_carried(e)) cycle
			i = is_met0 + 3*(e-1)
			nm0(e) = sden(i)
			nm1(e) = sden(i+1)
			nm2(e) = sden(i+2)
		enddo
	endif

	! Charge neutrality: every molecular ion carries +1, He++ and X++ carry
	! +2. Same sum, in the same order, as the fraction systems build.
	n_e = sden(is_HII) + sden(is_H2p) + sden(is_H3p) + sden(is_HeHp)      &
	    + sden(is_HeII) + 2.0d0*sden(is_HeIII)
	if (thereis_metals) call metal_electron_sum(n_e, met_nelem, nm1, nm2)

	! Free neutral helium, singlet plus metastable: the He I reservoir the
	! metal charge exchange reacts with (the metastable is a level inside
	! the neutral stage). The He <-> H pair below takes the ground singlet
	! alone, exactly as in the fraction systems.
	n_hei = sden(is_HeI_SI) + sden(is_HeITR)

	! --- reaction rows, in the fraction layout ---
	fvec(1:n_fraction_row) = 0.0d0

	! Only the radiation-driven rates carry lambda; every collisional,
	! recombination, charge-exchange and three-body rate is the cell's own.
	call mol_heh_rows(fvec, sden(is_HI), sden(is_HII), sden(is_H2),       &
	                  sden(is_H2p), sden(is_H3p), sden(is_HeHp),          &
	                  sden(is_HeI_SI), sden(is_HeITR), sden(is_HeII),     &
	                  sden(is_HeIII), n_e, ieq_cell%ntot,                 &
	                  lam*ieq_cell%P_HI, lam*ieq_cell%P_HeI,              &
	                  lam*ieq_cell%P_HeII, lam*ieq_cell%P_HeITR,          &
	                  lam*ieq_cell%P_H2, lam*ieq_cell%P_H2_di,             &
	                  lam*ieq_cell%P_H2_dd, lam*ieq_cell%P_H2_nd,          &
	                  lam*ieq_cell%k_LW,                                  &
	                  ieq_cell%rchiiB, ieq_cell%rcheiiB,                  &
	                  ieq_cell%rcheiiiB, ieq_cell%rcheiTR,                &
	                  ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,              &
	                  ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR,          &
	                  ieq_cell%q13, ieq_cell%q31a, ieq_cell%q31b,         &
	                  ieq_cell%Q31, ieq_cell%A31)

	! Oxygen carriers, and the oxygen cycle's exchange with the H2 row.
	! sden(is_HI) is the free atomic H and nm0(iel_O) the free atomic O by
	! construction here -- both are unknowns of their own, so neither has
	! to be formed as a difference.
	if (thereis_oxychem)                                                  &
		call oxygen_carrier_rows(fvec, oxygen_base, sden(is_HI),          &
		                         sden(is_H2), sden(is_OH), sden(is_H2O),  &
		                         nm0(iel_O),                              &
		                         lam*oj3, lam*oj4, lam*oj5, lam*oj7)

	if (thereis_metals) then
		do e = 1,met_nelem
			gm0(e) = lam*met_g0(e)
			gm1(e) = lam*met_g1(e)
		enddo
		! metal_rows writes the identity rows fvec = x for an absent
		! element and for the pinned X++ of a two-stage element; those
		! rows are not selected below, so the argument only has to exist.
		xzero(1:n_fraction_row) = 0.0d0
		call metal_rows(fvec, xzero, metal_base, met_nelem, met_ntot,     &
		                gm0, gm1, met_b0, met_b1, met_a1, met_a2,         &
		                met_top, nm0, nm1, nm2, n_e)
		cx_metal_base = metal_base
		call cx_add_to_fvec(n_fraction_row, fvec, nm0, nm1, nm2,          &
		                    sden(is_HI), sden(is_HII), n_hei,             &
		                    sden(is_HeII), sden(is_HeIII))
		cx_metal_base = 4
	endif

	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,     &
	                  sden(is_HI), sden(is_HII), sden(is_HeI_SI),         &
	                  sden(is_HeII), 1.0d0)

	! Each row against its own turnover rate, in the normalization the
	! caller's acceptance judge uses (set_turnover_scales_at_field).
	do i = 1,n_fraction_row
		fvec(i) = fvec(i)*row_scale(i)
	enddo

	end subroutine network_balance_rows

	!----------------------------------!

	subroutine element_conservation_rows(sden, cons)
	! The dimensionless element conservation statements, in the canonical
	! order: hydrogen, helium, then each element the cell carries. HeH+
	! appears in the hydrogen row AND in the helium row, and the oxygen
	! carriers appear in the hydrogen row AND in the oxygen row; both budgets
	! are satisfied at the root, which is the requirement that motivated this
	! formulation and the thing a rescaling of stage fractions cannot do.

	real*8, intent(in)  :: sden(n_species_max)
	real*8, intent(out) :: cons(n_conservation_row)
	real*8  :: s
	integer :: i, e, ir

	ir = 0

	! Hydrogen nuclei: free atomic H and H+, two per H2 and per H2+, three
	! per H3+, one in HeH+, and where the oxygen chemistry is on one in OH
	! and two in H2O.
	s = sden(is_HI) + sden(is_HII) + 2.0d0*sden(is_H2)                    &
	  + 2.0d0*sden(is_H2p) + 3.0d0*sden(is_H3p) + sden(is_HeHp)
	if (thereis_oxychem) s = s + sden(is_OH) + 2.0d0*sden(is_H2O)
	ir = ir + 1
	cons(ir) = s/nuclei_H - 1.0d0

	! Helium nuclei: the ground singlet, the 2^3S metastable, He+, He++, and
	! the He nucleus bound in HeH+.
	s = sden(is_HeI_SI) + sden(is_HeITR) + sden(is_HeII)                  &
	  + sden(is_HeIII) + sden(is_HeHp)
	ir = ir + 1
	cons(ir) = s/nuclei_He - 1.0d0

	! Each element the cell carries: its ionization stages, and for oxygen
	! with the carrier chemistry on the O nuclei bound in OH and H2O.
	if (thereis_metals) then
		do e = 1,met_nelem
			if (.not. element_carried(e)) cycle
			i = is_met0 + 3*(e-1)
			s = sden(i) + sden(i+1)
			if (met_top(e) .ge. 2) s = s + sden(i+2)
			if (thereis_oxychem .and. e .eq. iel_O)                &
				s = s + sden(is_OH) + sden(is_H2O)
			ir = ir + 1
			cons(ir) = s/met_ntot(e) - 1.0d0
		enddo
	endif

	end subroutine element_conservation_rows

	!----------------------------------!

	subroutine place_held_species(sden)
	! Restore every held species from its OWN balance row at the composition
	! the Newton solve returned (Gordon & McBride 1994, NASA RP-1311: the
	! species omitted from the iteration equations are recovered afterwards
	! from the equilibrium relations).
	!
	! No species of this network reacts with itself, so a species' balance
	! row is affine in its own density: row = a - b n, with a the production
	! and b the loss rate coefficient. Two evaluations of the row therefore
	! determine it exactly -- at n = 0 the row gives a, at a probe density it
	! gives a - b n_probe -- and the balancing density is a/b. (The only
	! coupling that is not affine runs through the electron density, and a
	! held species carries less than eps_hold of its element, so its own
	! charge changes n_e by less than that; the extraction is exact to the
	! same order that justifies holding it at all.)
	!
	! A row with no production or no loss has no balancing density to offer:
	! the species keeps the value it came in with, and the rung acceptance
	! over the full row set decides whether that is admissible.
	!
	! The held species feed each other -- H2+ makes H3+, HeH+ makes H3+ --
	! so the placement is swept over the set with the others held at their
	! current values.

	real*8, intent(inout) :: sden(n_species_max)
	real*8  :: fvec(n_fraction_rows_max)
	real*8  :: n_keep, n_probe, prod, loss, f_probe, n_el
	integer :: sweep, i, r

	do sweep = 1,n_placement_sweep
		do i = 1,n_species_max
			if (.not. species_held(i)) cycle
			r    = row_of_species(i)
			n_el = element_nuclei_of_species(i)
			if (n_el .le. 0.0d0) cycle
			n_keep  = sden(i)
			n_probe = max(n_keep, eps_hold*n_el)

			sden(i) = 0.0d0
			call network_balance_rows(sden, fvec)
			prod = fvec(r)

			sden(i) = n_probe
			call network_balance_rows(sden, fvec)
			f_probe = fvec(r)

			loss = (prod - f_probe)/n_probe
			if (prod .gt. 0.0d0 .and. loss .gt. 0.0d0) then
				sden(i) = prod/loss
			else
				sden(i) = n_keep
			endif
		enddo
	enddo

	end subroutine place_held_species

	!----------------------------------!

	double precision function full_network_residual(sden) result(res)
	! Largest scaled residual of the COMPLETE system at a composition: every
	! balance row the cell carries, the held species' rows included at their
	! placed densities, and every element conservation row. This is what a
	! rung has to satisfy. A Newton solve that converged over the reduced
	! system says nothing about the rows that were taken out of it, so the
	! reduced residual is not a rung acceptance and is not used as one.

	real*8, intent(in) :: sden(n_species_max)
	real*8  :: fvec(n_fraction_rows_max)
	real*8  :: cons(n_species_max)
	integer :: i, r

	call network_balance_rows(sden, fvec)

	res = 0.0d0
	do i = 1,n_species_max
		if (.not. species_exists(i)) cycle
		if (species_fixed(i))        cycle
		r = row_of_species(i)
		if (r .le. 0) cycle
		res = max(res, abs(fvec(r)))
	enddo

	call element_conservation_rows(sden, cons)
	do i = 1,n_conservation_row
		res = max(res, abs(cons(i)))
	enddo

	end function full_network_residual

	!----------------------------------!

	subroutine molecular_limit_composition(sden, p_bar, n_e_ref,          &
	                                       photo_scale, ok)
	! The composition of the cell in the DENSE MOLECULAR LIMIT at a
	! radiation field scaled by photo_scale, in closed form: the starting
	! point of the continuation.
	!
	!   * The H nuclei are partitioned between H2 and atomic H by the
	!     chemical-equilibrium mixing ratio q_H2(p, T) of Koskinen et al.
	!     (2022) Eq. 11 -- the same fit the molecular base boundary condition
	!     and the driver's molecular-basin retry use. Per H nucleus the
	!     fraction bound in H2 is 2 q (1 + He/H)/(1 + q). This part is the
	!     THERMOCHEMICAL limit and does not see photo_scale at all: at the
	!     weak fields this routine is called at, photodissociation is decades
	!     below the collisional dissociation equilibrium of a shielded base.
	!   * Every element is then placed in its own ionization balance AT THAT
	!     FIELD: the rate out of a stage is photo_scale*gamma + beta n_e and
	!     the rate back into it alpha n_e, which is
	!     ionization_equilibrium::ionization_balance_at_fixed_ne with its
	!     photoionization scaled. It is written here rather than imported
	!     because the imported one carries the unscaled field, and because
	!     importing it would make this module depend on the driver that
	!     calls it.
	!
	!     Keeping the photo term is what makes the starting point usable.
	!     Dropped, the trace ionization stages of a shielded base fall to the
	!     seed floor -- 22 decades below where even 1e-6 of the local field
	!     puts them -- and in logarithmic unknowns that is a singular
	!     starting point, not a nearby one (see the ladder above).
	!   * The oxygen carriers come from the chemical equilibrium of the same
	!     H2/H partition (oxygen_chemical_equilibrium_fractions), with free
	!     atomic oxygen and the oxygen ions sharing what the carriers leave.
	!   * The molecular ions H2+, H3+ and HeH+ are trace intermediates whose
	!     abundance is set by the very couplings this closed form drops, so
	!     it has nothing to say about them; they, and any other species it
	!     leaves empty, start at the seed floor.
	!
	! ok is false when the limit cannot be built at all, which is a cell
	! this network does not describe rather than a solver failure.

	real*8,  intent(out) :: sden(n_species_max)
	real*8,  intent(in)  :: p_bar, n_e_ref, photo_scale
	logical, intent(out) :: ok

	real*8  :: n_e, n_hi_seed, n_h2_seed
	real*8  :: f_oh, f_h2o, s_free, sden_floor
	integer :: e, i, isp, k

	ok      = .false.
	sden(:) = 0.0d0

	! CHARGE NEUTRALITY OF THE STARTING POINT. The ionization balances below
	! need an electron density, and the obvious one -- the cell's incoming
	! n_e, which is what the turnover scales use -- is the wrong one HERE. In
	! the shielded molecular base the trace metals are their own dominant
	! electron donors: measured on the hot-Uranus base, the balances at the
	! incoming n_e = 9.3e3 produce a composition whose own charge neutrality
	! gives 2.5e5, a factor of 27 apart, and every metal is then seeded
	! over-ionized by the same factor. The residual enforces neutrality
	! identically, so such a start is far from the root in a direction that
	! moves every metal at once, and the solve stalls on it.
	!
	! The electron density of the starting point is therefore brought to its
	! own fixed point: balances at n_e, neutrality of the result, repeat.
	! The map is a contraction here (more electrons means more recombination
	! means fewer ions means fewer electrons), and the geometric mean damping
	! keeps it one even where recombination is weak. This is the starting
	! point only; the solve owns the answer.
	n_e = n_e_ref
	do k = 1, n_neutrality_sweep
		call molecular_limit_stages(sden, p_bar, photo_scale, n_e)
		s_free = sden(is_HII) + sden(is_H2p) + sden(is_H3p)            &
		       + sden(is_HeHp) + sden(is_HeII) + 2.0d0*sden(is_HeIII)
		if (thereis_metals) then
			do e = 1,met_nelem
				if (.not. element_carried(e)) cycle
				i = is_met0 + 3*(e-1)
				s_free = s_free + sden(i+1) + 2.0d0*sden(i+2)
			enddo
		endif
		if (s_free .le. 0.0d0) exit
		n_e = sqrt(n_e*s_free)
	enddo
	call molecular_limit_stages(sden, p_bar, photo_scale, n_e)

	! --- oxygen carriers, and the oxygen budget they leave ---
	if (thereis_oxychem) then
		if (oxygen_carriers_are_fixed) then
			f_oh  = ieq_cell%x_oh_fix
			f_h2o = ieq_cell%x_h2o_fix
		else
			n_h2_seed = sden(is_H2)
			n_hi_seed = sden(is_HI)
			call oxygen_chemical_equilibrium_fractions(ieq_cell%T_K,  &
			                        n_h2_seed, n_hi_seed, f_oh, f_h2o)
		endif
		sden(is_OH)  = f_oh*nuclei_O
		sden(is_H2O) = f_h2o*nuclei_O
		! The carriers and the atomic/ionized oxygen share one budget: the
		! stages just written are scaled into what the carriers leave.
		s_free = max(1.0d0 - f_oh - f_h2o, 0.0d0)
		i = is_met0 + 3*(iel_O-1)
		sden(i)   = sden(i)*s_free
		sden(i+1) = sden(i+1)*s_free
		sden(i+2) = sden(i+2)*s_free
	endif

	! Probe-only: displace the closed-form start, to test that the solve
	! reaches the same root from a nearby seed rather than from this one.
	if (seed_perturb_on) then
		do isp = 1,n_species_max
			if (.not. seed_species_of_cell(isp)) cycle
			sden(isp) = sden(isp)*seed_perturb(isp)
		enddo
	endif

	! --- seed floors: a logarithmic unknown needs a positive start ---
	!
	! The set floored here is the CELL's, not a rung's. This routine builds
	! the seed, so it runs before set_rung_partition has chosen any rung for
	! this cell: n_unknown and species_of_unknown still describe whichever
	! cell this thread solved last. Reading them here made a cell's seed --
	! which species are floored, and, through the return below, whether the
	! continuation is seeded at all -- a function of the order the cells
	! happened to be handed to threads, and with it the accepted state of
	! every class-5 cell (docs/Update_EXHALE_stage1.md section 121).
	do isp = 1,n_species_max
		if (.not. seed_species_of_cell(isp)) cycle
		sden_floor = seed_floor_fraction*element_nuclei_of_species(isp)
		if (sden(isp) .lt. sden_floor) sden(isp) = sden_floor
		if (sden(isp) .le. 0.0d0) return
	enddo

	ok = .true.

	end subroutine molecular_limit_composition

	!----------------------------------!

	subroutine molecular_limit_stages(sden, p_bar, photo_scale, n_e)
	! The H2 partition and the ionization stages of every element at a GIVEN
	! electron density, at a radiation field scaled by photo_scale. Split out
	! of molecular_limit_composition so the starting point's own charge
	! neutrality can be closed by calling it repeatedly; see the fixed point
	! there for why that matters. The oxygen carriers and the seed floors are
	! not here: they sit outside the neutrality loop, the carriers because
	! they hold no charge and the floors because they must be applied once,
	! at the end.

	real*8,  intent(inout) :: sden(n_species_max)
	real*8,  intent(in)    :: p_bar, photo_scale, n_e

	real*8  :: qh2, x_h2, x_atomic
	real*8  :: f0, f1, f2, f_tr, drain_tr
	real*8  :: u0, u1, d1, d2
	integer :: e, i

	! --- hydrogen: H2 partition, then the atomic remainder ionized ---
	if (h2_is_fixed) then
		x_h2 = ieq_cell%x_h2_fix
	else
		qh2  = q_h2_equilibrium(p_bar, ieq_cell%T_K)
		x_h2 = 2.0d0*qh2*(1.0d0 + HeH)/(1.0d0 + qh2)
	endif
	if (x_h2 .gt. 1.0d0) x_h2 = 1.0d0
	if (x_h2 .lt. 0.0d0) x_h2 = 0.0d0
	x_atomic = 1.0d0 - x_h2

	sden(is_H2) = 0.5d0*x_h2*nuclei_H

	! WHERE THE IONIZATION STATE IS TRANSPORTED THE SEED USES IT. This
	! routine builds the starting point, and the row that will be solved for
	! x(1) is then not a balance but the constraint x(1) = x_hp_fix, so a
	! seed placed on the local balance instead starts the solve off its own
	! constraint surface and closes the seed's charge neutrality at an
	! electron density belonging to a different ionization state. Both are
	! avoided by seeding at the transported fraction; the H2 partition above
	! is handled the same way for the same reason.
	if (hp_is_fixed) then
		f1 = min(max(ieq_cell%x_hp_fix, 0.0d0), x_atomic)
		sden(is_HII) = f1*nuclei_H
		sden(is_HI)  = max(x_atomic*nuclei_H - sden(is_HII), 0.0d0)
	else
	u0 = photo_scale*ieq_cell%P_HI + ieq_cell%a_ion_HI*n_e
	d1 = ieq_cell%rchiiB*n_e
	call ionization_stage_fractions(u0, 0.0d0, d1, 0.0d0, 1, f0, f1, f2)
	sden(is_HI)  = f0*x_atomic*nuclei_H
	sden(is_HII) = f1*x_atomic*nuclei_H
	endif

	! --- helium: the three stages, then the metastable inside the neutral ---
	u0 = photo_scale*ieq_cell%P_HeI + ieq_cell%a_ion_HeI*n_e
	d1 = ieq_cell%rcheiiB*n_e
	! The metastable capture is a branch of the He+ recombination, so it
	! adds to the rate back into neutral He when the triplet is tracked.
	if (thereis_HeITR) d1 = d1 + ieq_cell%rcheiTR*n_e
	u1 = photo_scale*ieq_cell%P_HeII + ieq_cell%a_ion_HeII*n_e
	d2 = ieq_cell%rcheiiiB*n_e
	call ionization_stage_fractions(u0, u1, d1, d2, 2, f0, f1, f2)
	sden(is_HeII)  = f1*nuclei_He
	sden(is_HeIII) = f2*nuclei_He

	f_tr = 0.0d0
	if (thereis_HeITR) then
		! Steady state of the 2^3S level: fed by He+ recombination into the
		! triplet and by 1^1S collisional excitation, drained by
		! photoionization, A31, collisional de-excitation, electron impact
		! ionization and Penning ionization on H0.
		drain_tr = photo_scale*ieq_cell%P_HeITR + ieq_cell%A31         &
		         + sden(is_HI)*ieq_cell%Q31                            &
		         + (ieq_cell%q31a + ieq_cell%q31b                      &
		            + ieq_cell%a_ion_HeITR)*n_e
		if (drain_tr .gt. 0.0d0) f_tr = min(                           &
			n_e*(f1*ieq_cell%rcheiTR + f0*ieq_cell%q13)/drain_tr, f0)
	endif
	sden(is_HeITR)  = f_tr*nuclei_He
	sden(is_HeI_SI) = (f0 - f_tr)*nuclei_He

	! --- metals: each element in its own ionization balance ---
	if (thereis_metals) then
		do e = 1,met_nelem
			if (.not. element_carried(e)) cycle
			i  = is_met0 + 3*(e-1)
			u0 = photo_scale*met_g0(e) + met_b0(e)*n_e
			d1 = met_a1(e)*n_e
			u1 = photo_scale*met_g1(e) + met_b1(e)*n_e
			d2 = met_a2(e)*n_e
			call ionization_stage_fractions(u0, u1, d1, d2,        &
			                                met_top(e), f0, f1, f2)
			sden(i)   = f0*met_ntot(e)
			sden(i+1) = f1*met_ntot(e)
			if (met_top(e) .ge. 2) sden(i+2) = f2*met_ntot(e)
		enddo
	endif

	end subroutine molecular_limit_stages

	!----------------------------------!

	double precision function element_nuclei_of_species(isp) result(n_el)
	! Nuclei density of the element whose budget a species belongs to. The
	! species shared between two elements (HeH+, and the oxygen carriers)
	! are referred to the hydrogen reservoir, which is the larger of their
	! two budgets and therefore the weaker floor.
	integer, intent(in) :: isp
	integer :: e

	if (isp .le. is_H3p .or. isp .eq. is_HeHp) then
		n_el = nuclei_H
	else if (isp .le. is_HeITR) then
		n_el = nuclei_He
	else if (isp .le. is_H2O) then
		n_el = nuclei_H
	else
		e    = (isp - is_met0)/3 + 1
		n_el = met_ntot(e)
	endif

	end function element_nuclei_of_species

	!----------------------------------!

	subroutine ionization_stage_fractions(u0, u1, d1, d2, top, f0, f1, f2)
	! Stage fractions of an element in its OWN ionization balance at a fixed
	! electron density,
	!
	!    n_k (gamma_k + beta_k n_e) = alpha_{k+1} n_e n_{k+1},
	!
	! with the couplings between elements dropped. With u_k the rate out of
	! stage k (photoionization plus electron impact) and d_k the rate back
	! into stage k-1 (radiative recombination), the three-stage solution is
	!
	!    (n_0, n_1, n_2) proportional to (d_1 d_2, u_0 d_2, u_0 u_1),
	!
	! non-negative and normalized to one, so it always lies inside the
	! element's simplex. top = 1 keeps two stages, top >= 2 all three. A
	! vanishing denominator means the element has no ionization channel at
	! all at this field, and the whole element is neutral.
	real*8,  intent(in)  :: u0, u1, d1, d2
	integer, intent(in)  :: top
	real*8,  intent(out) :: f0, f1, f2
	real*8 :: w0, w1, w2, s

	f0 = 1.0d0
	f1 = 0.0d0
	f2 = 0.0d0

	if (top .ge. 2) then
		w0 = d1*d2
		w1 = u0*d2
		w2 = u0*u1
		s  = w0 + w1 + w2
		if (s .gt. 0.0d0) then
			f0 = w0/s
			f1 = w1/s
			f2 = w2/s
		endif
	else
		s = d1 + u0
		if (s .gt. 0.0d0) then
			f0 = d1/s
			f1 = u0/s
		endif
	endif

	end subroutine ionization_stage_fractions

	!----------------------------------!

	subroutine fraction_layout_from_densities(sden, x, nx)
	! Species densities back to the stage fractions of the cell's own
	! fraction-layout system (System_HeH_mol / System_HeH_mol_metals), which
	! is the layout every downstream consumer -- the acceptance judge, the
	! extraction, the output -- reads. Layout slots the cell does not carry
	! (an absent element, the pinned X++ of a two-stage element) are left at
	! zero, matching the identity rows those slots have in the residual.

	real*8,  intent(in)  :: sden(n_species_max)
	integer, intent(in)  :: nx
	real*8,  intent(out) :: x(nx)
	integer :: e, i, ix

	x(1:nx) = 0.0d0

	x(1) = sden(is_HII)/nuclei_H
	x(2) = sden(is_HeII)/nuclei_He
	x(3) = sden(is_HeIII)/nuclei_He
	x(4) = 2.0d0*sden(is_H2)/nuclei_H
	x(5) = 2.0d0*sden(is_H2p)/nuclei_H
	x(6) = 3.0d0*sden(is_H3p)/nuclei_H
	x(7) = sden(is_HeHp)/nuclei_H
	if (thereis_HeITR) x(8) = sden(is_HeITR)/nuclei_He

	if (thereis_oxychem .and. nuclei_O .gt. 0.0d0) then
		x(oxygen_base)   = sden(is_OH)/nuclei_O
		x(oxygen_base+1) = sden(is_H2O)/nuclei_O
	endif

	if (thereis_metals) then
		do e = 1,met_nelem
			if (.not. element_carried(e)) cycle
			i  = is_met0 + 3*(e-1)
			ix = metal_base + 2*(e-1)
			x(ix) = sden(i+1)/met_ntot(e)
			if (met_top(e) .ge. 2) x(ix+1) = sden(i+2)/met_ntot(e)
		enddo
	endif

	end subroutine fraction_layout_from_densities

	! ===============================================================
	!  Opt-in diagnostic of a continuation that did not reach the full
	!  field (docs/charge_exchange_cancellation_limit.md section 5).
	!  Nothing below runs unless EXHALE_CCE_DUMP is set, or unless the
	!  standalone driver calls cce_probe_from_dump. No production path
	!  reaches any of it.
	! ===============================================================

	subroutine write_continuation_dump(envname, already, n_e_ref, p_bar,   &
	                                   nx, mbase, iox, lam_ref, lam_try,   &
	                                   sden, info_loc)
	! Save everything a standalone evaluation of this cell needs. The rate
	! COEFFICIENTS are not written: they are deterministic functions of the
	! temperature and the total particle density (set_mol_coeffs), of the
	! metal rate arrays (set_metal_coeffs) and of the temperature again
	! (cx_set_cell), so the driver regenerates them from the same inputs by
	! calling the same routines. Writing the inputs rather than the derived
	! values is what makes the standalone reproduction a check rather than a
	! transcription.

	character(len=*), intent(in) :: envname
	logical, intent(inout) :: already
	real*8,  intent(in) :: n_e_ref, p_bar, lam_ref, lam_try
	integer, intent(in) :: nx, mbase, iox, info_loc
	real*8,  intent(in) :: sden(n_species_max)
	integer :: iu, ios, e

	if (already) return
	call get_environment_variable(envname, cce_dump_path)
	if (len_trim(cce_dump_path) .eq. 0) return

	!$omp critical (cce_dump_guard)
	if (.not. already) then
		already = .true.
		open(newunit=iu, file=trim(cce_dump_path), status='replace',   &
		     action='write', iostat=ios)
		if (ios .eq. 0) then
			write(iu,*) 1                      ! dump format version
			write(iu,*) thereis_HeITR, thereis_metals, thereis_oxychem
			write(iu,*) HeH
			write(iu,*) nx, mbase, iox, info_loc
			write(iu,*) n_e_ref, p_bar, lam_ref, lam_try
			write(iu,*) ieq_cell%P_HI, ieq_cell%P_HeI,             &
			            ieq_cell%P_HeII, ieq_cell%rchiiB,          &
			            ieq_cell%rcheiiB, ieq_cell%rcheiiiB
			write(iu,*) ieq_cell%nh, ieq_cell%nhe
			write(iu,*) ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,     &
			            ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR
			write(iu,*) ieq_cell%rcheiTR, ieq_cell%A31,            &
			            ieq_cell%P_HeITR
			write(iu,*) ieq_cell%q13, ieq_cell%q31a,               &
			            ieq_cell%q31b, ieq_cell%Q31
			write(iu,*) ieq_cell%P_H2, ieq_cell%P_H2_di,           &
			            ieq_cell%P_H2_dd, ieq_cell%P_H2_nd,         &
			            ieq_cell%k_LW,                              &
			            ieq_cell%T_K, ieq_cell%ntot
			write(iu,*) ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0
			write(iu,*) ieq_cell%n_ofam, ieq_cell%n_co
			write(iu,*) ieq_cell%x_h2_fixed, ieq_cell%x_ox_fixed,   &
			            ieq_cell%x_h2_fix,                          &
			            ieq_cell%x_oh_fix, ieq_cell%x_h2o_fix
			write(iu,*) ieq_cell%x_hp_fixed, ieq_cell%x_hp_fix
			write(iu,*) met_nelem
			do e = 1,met_nelem
				write(iu,*) met_ntot(e), met_g0(e), met_g1(e),     &
				            met_b0(e), met_b1(e), met_a1(e),       &
				            met_a2(e), met_top(e)
			enddo
			write(iu,*) sden(1:n_species_max)
			close(iu)
			write(*,'(A,A,A,A)') ' (constrained_chemical_'//   &
				'equilibrium) cell state (', envname,               &
				') written to ', trim(cce_dump_path)
			flush(6)
		endif
	endif
	!$omp end critical (cce_dump_guard)

	end subroutine write_continuation_dump

	!----------------------------------!

	subroutine hydrogen_helium_row_terms(sden, lbl, val, nt, irow)
	! Every signed term of the H+ balance (irow = 1) or the He+ balance
	! (irow = 2) BEFORE the turnover scaling, as the review plan asks for.
	!
	! This MIRRORS mol_heh_rows term by term and is a diagnostic only: it
	! defines no rate coefficient of its own -- every coefficient is read
	! from the one place that owns it (set_mol_coeffs in System_HeH_mol, and
	! ieq_cell) -- but it does restate the row's STRUCTURE, so a reaction
	! added to mol_heh_rows has to be added here too or this decomposition
	! stops summing to the row. The probe checks exactly that: it compares
	! the sum of these terms against the row the shared assembly produces,
	! and says so when they differ.

	real*8, intent(in)  :: sden(n_species_max)
	character(len=*), intent(out) :: lbl(*)
	real*8,  intent(out) :: val(*)
	integer, intent(out) :: nt
	integer, intent(in)  :: irow

	real*8 :: n_e, n_hei, lam, r1, r2
	real*8 :: nm0(n_melem), nm1(n_melem), nm2(n_melem)
	integer :: e, i

	lam = field_scale

	nm0(:) = 0.0d0
	nm1(:) = 0.0d0
	nm2(:) = 0.0d0
	if (thereis_metals) then
		do e = 1,met_nelem
			if (.not. element_carried(e)) cycle
			i = is_met0 + 3*(e-1)
			nm0(e) = sden(i)
			nm1(e) = sden(i+1)
			nm2(e) = sden(i+2)
		enddo
	endif
	n_e = sden(is_HII) + sden(is_H2p) + sden(is_H3p) + sden(is_HeHp)      &
	    + sden(is_HeII) + 2.0d0*sden(is_HeIII)
	if (thereis_metals) call metal_electron_sum(n_e, met_nelem, nm1, nm2)
	n_hei = sden(is_HeI_SI) + sden(is_HeITR)

	! The H <-> He pair, the two gross rates the review is about.
	r1 = ieq_cell%kcx_He0_Hp*sden(is_HeI_SI)*sden(is_HII)
	r2 = ieq_cell%kcx_Hep_H0*sden(is_HeII)*sden(is_HI)

	nt = 0
	if (irow .eq. 1) then
		call add_term(lbl, val, nt, 'photoion H0',                     &
			lam*ieq_cell%P_HI*sden(is_HI))
		call add_term(lbl, val, nt, 'collion H0',                      &
			ieq_cell%a_ion_HI*n_e*sden(is_HI))
		call add_term(lbl, val, nt, 'dissoc photoion H2',              &
			lam*ieq_cell%P_H2_di*sden(is_H2))
		! Two protons per event, which is why this term carries a 2 and
		! the H2 row does not.
		call add_term(lbl, val, nt, 'double photoion H2',             &
			2.0d0*lam*ieq_cell%P_H2_dd*sden(is_H2))
		call add_term(lbl, val, nt, 'R9  H2+ + H0',                    &
			mk9*sden(is_H2p)*sden(is_HI))
		call add_term(lbl, val, nt, 'R17 He+ + H2',                    &
			mk17*sden(is_HeII)*sden(is_H2))
		call add_term(lbl, val, nt, 'Penning He(2^3S)+H0',             &
			f_penning_HeI23S*ieq_cell%Q31*sden(is_HeITR)*sden(is_HI))
		call add_term(lbl, val, nt, 'recomb H+',                       &
			-ieq_cell%rchiiB*n_e*sden(is_HII))
		call add_term(lbl, val, nt, 'R10+R13 H+ + H2',                 &
			-(mk10 + mk13)*sden(is_HII)*sden(is_H2))
		call add_term(lbl, val, nt, 'CX gross R2 (He+ +H0)', r2)
		call add_term(lbl, val, nt, 'CX gross R1 (He0 +H+)', -r1)
	else
		call add_term(lbl, val, nt, 'photoion He(1^1S)',               &
			lam*ieq_cell%P_HeI*sden(is_HeI_SI))
		call add_term(lbl, val, nt, 'collion He(1^1S)',                &
			ieq_cell%a_ion_HeI*n_e*sden(is_HeI_SI))
		call add_term(lbl, val, nt, 'photoion He(2^3S)',               &
			lam*ieq_cell%P_HeITR*sden(is_HeITR))
		call add_term(lbl, val, nt, 'collion He(2^3S)',                &
			ieq_cell%a_ion_HeITR*n_e*sden(is_HeITR))
		call add_term(lbl, val, nt, 'recomb He++',                     &
			ieq_cell%rcheiiiB*n_e*sden(is_HeIII))
		call add_term(lbl, val, nt, 'recomb He+ (B+TR)',               &
			-(ieq_cell%rcheiiB + ieq_cell%rcheiTR)*n_e*sden(is_HeII))
		call add_term(lbl, val, nt, 'photoion He+',                    &
			-lam*ieq_cell%P_HeII*sden(is_HeII))
		call add_term(lbl, val, nt, 'collion He+',                     &
			-ieq_cell%a_ion_HeII*n_e*sden(is_HeII))
		call add_term(lbl, val, nt, 'R17+R20+R23 He+ + H2',            &
			-(mk17 + mk20 + mk23)*sden(is_HeII)*sden(is_H2))
		! D1 He+ + CO -> C+ + O + He (RATE22 4068), the He+ loss the row
		! of mol_heh_rows carries since B3b-CO2; n_co is the cell's CO
		! background density (zero unless the oxygen chemistry is on).
		call add_term(lbl, val, nt, 'D1  He+ + CO',                     &
			-rk_D1_Hep_CO()*ieq_cell%n_co*sden(is_HeII))
		call add_term(lbl, val, nt, 'CX gross R1 (He0 +H+)', r1)
		call add_term(lbl, val, nt, 'CX gross R2 (He+ +H0)', -r2)
	endif

	end subroutine hydrogen_helium_row_terms

	!----------------------------------!

	subroutine add_term(lbl, val, nt, name, v)
	character(len=*), intent(inout) :: lbl(*)
	real*8,  intent(inout) :: val(*)
	integer, intent(inout) :: nt
	character(len=*), intent(in) :: name
	real*8,  intent(in) :: v
	nt = nt + 1
	lbl(nt) = name
	val(nt) = v
	end subroutine add_term

	!----------------------------------!

	subroutine jacobian_by_differences(u, fbase, fjac, rule, eta)
	! Numerical Jacobian of the reduced constrained system at u, by one of
	! three difference rules:
	!   rule = 0  MINPACK's unfloored rule, h = sqrt(eps)*|u_j|, h = sqrt(eps)
	!             when u_j is EXACTLY zero (fdjac1.f90 with
	!             unit_step_floor = .false., which is what every fraction
	!             system still gets through hybrd1);
	!   rule = 1  scale-aware forward, h = eta*max(1,|u_j|) -- the rule this
	!             module's own solve now uses (unit_step_floor = .true.);
	!   rule = 2  scale-aware central, same h.
	! The unknowns are logarithmic DEPARTURES, so rule 0 collapses as soon
	! as a component leaves zero; that is the hypothesis the two rules were
	! written to separate.

	real*8,  intent(in)  :: u(n_species_max), fbase(n_species_max)
	real*8,  intent(out) :: fjac(n_species_max,n_species_max)
	integer, intent(in)  :: rule
	real*8,  intent(in)  :: eta

	real*8 :: uw(n_species_max), fp(n_species_max), fm(n_species_max)
	real*8 :: par(60), h, epsm, dpmpar
	integer :: j, i, iflag

	par(:) = 0.0d0
	iflag  = 1
	epsm   = sqrt(dpmpar(1))

	do j = 1,n_unknown
		if (rule .eq. 0) then
			h = epsm*abs(u(j))
			if (h .eq. 0.0d0) h = epsm
		else
			h = eta*max(1.0d0, abs(u(j)))
		endif

		uw(1:n_unknown) = u(1:n_unknown)
		uw(j) = u(j) + h
		call constrained_equilibrium_residual(n_unknown, uw, fp,       &
		                                      iflag, par)
		if (rule .eq. 2) then
			uw(1:n_unknown) = u(1:n_unknown)
			uw(j) = u(j) - h
			call constrained_equilibrium_residual(n_unknown, uw, fm,  &
			                                      iflag, par)
			do i = 1,n_unknown
				fjac(i,j) = (fp(i) - fm(i))/(2.0d0*h)
			enddo
		else
			do i = 1,n_unknown
				fjac(i,j) = (fp(i) - fbase(i))/h
			enddo
		endif
	enddo

	end subroutine jacobian_by_differences

	!----------------------------------!

	subroutine cce_probe_from_dump(fname)
	! Standalone re-evaluation of a saved failed continuation: phase B of
	! docs/charge_exchange_cancellation_limit.md. Reads the dump, rebuilds
	! the cell exactly as the driver would, and reports
	!   1 the residual, and the H+/He+ rows summed in binary64 against the
	!     exact sum of the same binary64 terms (real*16 accumulation), which
	!     isolates SUMMATION error -- the thing section 3.1 says must be
	!     measured;
	!   2 Jacobian columns under MINPACK's rule, scale-aware forward and
	!     central differences over a sweep of eta, and the column-by-column
	!     disagreement between them;
	!   3 singular values of the scale-aware Jacobian;
	!   4 the same with the H <-> He pair removed, as a diagnosis only.

	character(len=*), intent(in) :: fname

	real*8  :: n_e_ref, p_bar, lam_ref, lam_try, HeH_in
	integer :: nx, mbase, iox, info_in, ver, nel
	real*8  :: sden(n_species_max), u(n_species_max), fres(n_species_max)
	real*8  :: mg_ntot(n_melem), mg_g0(n_melem), mg_g1(n_melem)
	real*8  :: mg_b0(n_melem), mg_b1(n_melem), mg_a1(n_melem)
	real*8  :: mg_a2(n_melem)
	integer :: mg_top(n_melem)
	real*8  :: fj0(n_species_max,n_species_max)
	real*8  :: fjf(n_species_max,n_species_max)
	real*8  :: fjc(n_species_max,n_species_max)
	real*8  :: fjref(n_species_max,n_species_max)
	real*8  :: fvec(n_fraction_rows_max)
	real*8  :: tval(24), etas(4), svals(n_species_max)
	real*8  :: wrk(8*n_species_max)
	character(len=24) :: tlbl(24)
	real*8  :: dsum, csum, cn0, cnf, cnc, dmax0, dmaxc, par(60)
	real*8  :: xfrac(n_fraction_rows_max), xref(n_fraction_rows_max)
	logical :: okfull
	integer :: nsolv
	real(kind=16) :: qsum
	integer :: iu, ios, e, i, j, k, nt, ie, iflag, lwrk, ierr
	real*8  :: kcx1_save, kcx2_save

	open(newunit=iu, file=trim(fname), status='old', action='read',       &
	     iostat=ios)
	if (ios .ne. 0) then
		write(*,'(A,A)') ' cce_probe: cannot open ', trim(fname)
		return
	endif
	read(iu,*) ver
	if (ver .ne. 1) then
		write(*,'(A,I0)') ' cce_probe: unsupported dump version ', ver
		close(iu)
		return
	endif
	read(iu,*) thereis_HeITR, thereis_metals, thereis_oxychem
	read(iu,*) HeH_in
	HeH = HeH_in
	read(iu,*) nx, mbase, iox, info_in
	read(iu,*) n_e_ref, p_bar, lam_ref, lam_try
	read(iu,*) ieq_cell%P_HI, ieq_cell%P_HeI, ieq_cell%P_HeII,            &
	           ieq_cell%rchiiB, ieq_cell%rcheiiB, ieq_cell%rcheiiiB
	read(iu,*) ieq_cell%nh, ieq_cell%nhe
	read(iu,*) ieq_cell%a_ion_HI, ieq_cell%a_ion_HeI,                     &
	           ieq_cell%a_ion_HeII, ieq_cell%a_ion_HeITR
	read(iu,*) ieq_cell%rcheiTR, ieq_cell%A31, ieq_cell%P_HeITR
	read(iu,*) ieq_cell%q13, ieq_cell%q31a, ieq_cell%q31b, ieq_cell%Q31
	read(iu,*) ieq_cell%P_H2, ieq_cell%P_H2_di,                           &
	           ieq_cell%P_H2_dd, ieq_cell%P_H2_nd, ieq_cell%k_LW,         &
	           ieq_cell%T_K, ieq_cell%ntot
	read(iu,*) ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0
	read(iu,*) ieq_cell%n_ofam, ieq_cell%n_co
	read(iu,*) ieq_cell%x_h2_fixed, ieq_cell%x_ox_fixed,                  &
	           ieq_cell%x_h2_fix,                                         &
	           ieq_cell%x_oh_fix, ieq_cell%x_h2o_fix
	! The transported-proton constraint, on its own record so that a dump
	! written before the proton option is still read by the lines above and
	! only this one fails -- the dump is a debug artifact of one cell, not a
	! restart format, so it is versioned by being appended to.
	read(iu,*) ieq_cell%x_hp_fixed, ieq_cell%x_hp_fix
	read(iu,*) nel
	do e = 1,nel
		read(iu,*) mg_ntot(e), mg_g0(e), mg_g1(e), mg_b0(e), mg_b1(e), &
		           mg_a1(e), mg_a2(e), mg_top(e)
	enddo
	read(iu,*) sden(1:n_species_max)
	close(iu)

	if (thereis_oxychem) then
		write(*,'(A)') ' cce_probe: the oxygen carrier chemistry is '// &
			'on in this dump. Its photolysis channels are cell '//     &
			'state that this dump does not carry, so the probe '//     &
			'refuses rather than evaluating a different cell.'
		return
	endif

	! Rebuild the cell exactly as ioniz_eq does: the same routines, from
	! the same inputs, so the coefficients are regenerated and not copied.
	call set_mol_coeffs(ieq_cell%T_K, ieq_cell%ntot)
	if (thereis_metals) then
		call set_metal_coeffs(nel, mg_ntot, mg_g0, mg_g1, mg_b0,       &
		                      mg_b1, mg_a1, mg_a2, mg_top)
		call cx_set_cell(ieq_cell%T_K)
	endif

	call set_molecular_network_layout(nx, mbase, iox)
	field_scale = lam_try
	call set_turnover_scales_at_field(n_e_ref, lam_try)
	call set_rung_partition(sden)
	reference_density(:) = sden(:)
	u(1:n_unknown) = u_shift

	write(*,'(A)') '=============================================='
	write(*,'(A)') ' constrained-equilibrium probe (phase B)'
	write(*,'(A)') '=============================================='
	write(*,'(A,ES13.6,A,ES13.6)') ' T [K] ', ieq_cell%T_K,               &
		'   n_tot [cm^-3] ', ieq_cell%ntot
	write(*,'(A,ES13.6,A,ES13.6)') ' n_H   ', ieq_cell%nh,                &
		'   n_He           ', ieq_cell%nhe
	write(*,'(A,ES13.6,A,ES13.6)') ' n_e_ref ', n_e_ref,                  &
		' p [bar]        ', p_bar
	write(*,'(A,ES13.6,A,ES13.6)') ' lambda_ref ', lam_ref,               &
		' lambda_try ', lam_try
	write(*,'(A,I0,A,I0,A,I0)') ' unknowns ', n_unknown, '  reaction rows ',&
		n_reaction_row, '  conservation rows ', n_conservation_row
	write(*,'(A,I0)') ' held species ', count(species_held)
	write(*,'(A,ES13.6)') ' kcx He0+H+ ', ieq_cell%kcx_He0_Hp
	write(*,'(A,ES13.6)') ' kcx He+ +H0 ', ieq_cell%kcx_Hep_H0

	! ---- 1. residual, and the summation test on the two H/He rows ----
	par(:) = 0.0d0
	iflag  = 1
	call constrained_equilibrium_residual(n_unknown, u, fres, iflag, par)
	call network_balance_rows(sden, fvec)

	write(*,'(A)') ' --- H+ and He+ balances, signed terms before scaling'
	do ie = 1,2
		call hydrogen_helium_row_terms(sden, tlbl, tval, nt, ie)
		dsum = 0.0d0
		qsum = 0.0_16
		do k = 1,nt
			dsum = dsum + tval(k)
			qsum = qsum + real(tval(k), kind=16)
		enddo
		if (ie .eq. 1) then
			write(*,'(A)') '   row 1 (H+):'
		else
			write(*,'(A)') '   row 2 (He+):'
		endif
		do k = 1,nt
			write(*,'(A,A24,ES23.15)') '     ', tlbl(k), tval(k)
		enddo
		csum = real(qsum, kind=8)
		write(*,'(A,ES23.15)') '     binary64 sum of terms   ', dsum
		write(*,'(A,ES23.15)') '     exact sum of same terms ', csum
		if (abs(csum) .gt. 0.0d0) then
			write(*,'(A,ES12.4)') '     summation error rel.  ',    &
				abs(dsum - csum)/abs(csum)
		endif
		write(*,'(A,ES23.15)') '     shared assembly (scaled)',      &
			fvec(ie)
		write(*,'(A,ES23.15)') '     terms scaled            ',      &
			dsum*row_scale(ie)
	enddo

	write(*,'(A,ES12.4)') ' full network residual ',                      &
		full_network_residual(sden)

	! ---- 2. Jacobian rules ----
	etas = (/ 1.0d-4, 1.0d-5, 1.0d-6, 1.0d-7 /)
	call jacobian_by_differences(u, fres, fj0, 0, 0.0d0)
	call jacobian_by_differences(u, fres, fjref, 2, 1.0d-5)

	write(*,'(A)') ' --- Jacobian column comparison (reference: central,'//&
		' eta = 1e-5)'
	write(*,'(A)') '   col  species   ||MINPACK||   ||central||   '//     &
		'rel.diff'
	dmax0 = 0.0d0
	do j = 1,n_unknown
		cn0 = 0.0d0
		cnc = 0.0d0
		do i = 1,n_unknown
			cn0 = cn0 + fj0(i,j)**2
			cnc = cnc + fjref(i,j)**2
		enddo
		cn0 = sqrt(cn0)
		cnc = sqrt(cnc)
		dsum = 0.0d0
		do i = 1,n_unknown
			dsum = dsum + (fj0(i,j) - fjref(i,j))**2
		enddo
		dsum = sqrt(dsum)
		if (cnc .gt. 0.0d0) dsum = dsum/cnc
		dmax0 = max(dmax0, dsum)
		write(*,'(I6,I9,3ES14.5)') j, species_of_unknown(j), cn0, cnc,  &
			dsum
	enddo
	write(*,'(A,ES12.4)') '   worst column disagreement (MINPACK rule) ', &
		dmax0

	write(*,'(A)') ' --- scale-aware step sweep (forward and central)'
	do k = 1,4
		call jacobian_by_differences(u, fres, fjf, 1, etas(k))
		call jacobian_by_differences(u, fres, fjc, 2, etas(k))
		dmaxc = 0.0d0
		dsum  = 0.0d0
		do j = 1,n_unknown
			cnf = 0.0d0
			cnc = 0.0d0
			do i = 1,n_unknown
				cnf = cnf + (fjf(i,j) - fjref(i,j))**2
				cnc = cnc + (fjc(i,j) - fjref(i,j))**2
			enddo
			cn0 = 0.0d0
			do i = 1,n_unknown
				cn0 = cn0 + fjref(i,j)**2
			enddo
			cn0 = sqrt(cn0)
			if (cn0 .gt. 0.0d0) then
				dsum  = max(dsum,  sqrt(cnf)/cn0)
				dmaxc = max(dmaxc, sqrt(cnc)/cn0)
			endif
		enddo
		write(*,'(A,ES9.2,A,ES12.4,A,ES12.4)') '   eta ', etas(k),      &
			'  forward vs ref ', dsum, '  central vs ref ', dmaxc
	enddo

	! ---- 3. singular values of the scale-aware Jacobian ----
	lwrk = 8*n_species_max
	call dgesvd('N','N', n_unknown, n_unknown, fjref, n_species_max,      &
	            svals, fjf, n_species_max, fjc, n_species_max, wrk,       &
	            lwrk, ierr)
	if (ierr .eq. 0) then
		write(*,'(A)') ' --- singular values of the scale-aware Jacobian'
		write(*,'(A,ES12.4,A,ES12.4)') '   largest ', svals(1),        &
			'   smallest ', svals(n_unknown)
		if (svals(n_unknown) .gt. 0.0d0)                               &
			write(*,'(A,ES12.4)') '   condition number ',          &
				svals(1)/svals(n_unknown)
		write(*,'(A)') '   spectrum:'
		write(*,'(6ES13.4)') svals(1:n_unknown)
	else
		write(*,'(A,I0)') ' dgesvd failed, info = ', ierr
	endif

	! ---- 4. the same with the H <-> He pair removed (diagnosis only) ----
	kcx1_save = ieq_cell%kcx_He0_Hp
	kcx2_save = ieq_cell%kcx_Hep_H0
	ieq_cell%kcx_He0_Hp = 0.0d0
	ieq_cell%kcx_Hep_H0 = 0.0d0
	call set_turnover_scales_at_field(n_e_ref, lam_try)
	call constrained_equilibrium_residual(n_unknown, u, fres, iflag, par)
	call jacobian_by_differences(u, fres, fjref, 2, 1.0d-5)
	call dgesvd('N','N', n_unknown, n_unknown, fjref, n_species_max,      &
	            svals, fjf, n_species_max, fjc, n_species_max, wrk,       &
	            lwrk, ierr)
	write(*,'(A)') ' --- H <-> He pair removed (diagnosis only, not a fix)'
	write(*,'(A,ES12.4)') '   full network residual ',                    &
		full_network_residual(sden)
	if (ierr .eq. 0) then
		write(*,'(A,ES12.4,A,ES12.4)') '   largest ', svals(1),        &
			'   smallest ', svals(n_unknown)
		if (svals(n_unknown) .gt. 0.0d0)                               &
			write(*,'(A,ES12.4)') '   condition number ',          &
				svals(1)/svals(n_unknown)
	endif
	ieq_cell%kcx_He0_Hp = kcx1_save
	ieq_cell%kcx_Hep_H0 = kcx2_save
	call set_turnover_scales_at_field(n_e_ref, lam_try)

	! ---- 5. same root from a displaced seed (validation criterion 4) ----
	! Run the WHOLE continuation, from the closed-form start and from that
	! start displaced by +-10% and by an alternating pattern, and compare
	! the compositions it lands on. A root reached only from one particular
	! seed is not a root of the cell, it is an artifact of the path.
	write(*,'(A)') ' --- seed robustness: full continuation from '//      &
		'displaced starts'
	do k = 0,3
		seed_perturb_on = (k .gt. 0)
		do i = 1,n_species_max
			select case (k)
			case (1) ; seed_perturb(i) = 1.1d0
			case (2) ; seed_perturb(i) = 0.9d0
			case default
				if (mod(i,2) .eq. 0) then
					seed_perturb(i) = 1.1d0
				else
					seed_perturb(i) = 0.9d0
				endif
			end select
		enddo
		call equilibrium_from_molecular_limit(xfrac, nx, mbase, iox,    &
		                                      n_e_ref, p_bar, okfull,   &
		                                      nsolv)
		if (k .eq. 0) then
			xref(1:nx) = xfrac(1:nx)
			write(*,'(A,L1,A,I0)') '   unperturbed   reached ',     &
				okfull, '  field solves ', nsolv
		else
			dsum = 0.0d0
			do i = 1,nx
				dsum = max(dsum, abs(xfrac(i) - xref(i))            &
				       /max(abs(xref(i)), 1.0d-30))
			enddo
			write(*,'(A,F5.2,A,L1,A,I0,A,ES11.4)') '   seed x ',    &
				merge(1.1d0, merge(0.9d0, -1.0d0, k .eq. 2),        &
				      k .eq. 1),                                    &
				'  reached ', okfull, '  solves ', nsolv,           &
				'  max rel. difference vs unperturbed ', dsum
		endif
	enddo
	seed_perturb_on = .false.

	write(*,'(A)') '=============================================='

	end subroutine cce_probe_from_dump

	! End of module
	end module constrained_chemical_equilibrium

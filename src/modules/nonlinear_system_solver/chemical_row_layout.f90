	module chemical_row_layout
	! THE ROWS OF THE CELL'S IONIZATION AND CHEMICAL SYSTEMS, WHAT EACH ROW
	! STATES, AND THE ELEMENT GROUP EACH KINETIC ROW IS MEASURED AGAINST.
	!
	! Two records, kept apart on purpose (md/AM_interface_proposal_20261006.md,
	! plan md/PLAN_20261006_to_be_done_AI_AM_rev3.md section 2):
	!
	!   1. the NUCLEAR STOICHIOMETRY of every species slot, the number of
	!      nuclei of each element a species carries (nuclei_of_species):
	!      H2 two H, H3+ three H, HeH+ one H and one He, OH one H and one O,
	!      H2O two H and one O, CO one C and one O, a metal stage one nucleus
	!      of its element. It is physics: element conservation.
	!   2. the ACCEPTANCE GROUP of every row (acceptance_group): the element
	!      whose chemical turnover a kinetic row's imbalance is measured
	!      against when a candidate composition is judged as a root. It is a
	!      normalization choice, not the stoichiometry: HeH+ (H and He) is
	!      measured against hydrogen and OH, H2O (H and O) against oxygen.
	!      The only tie between the two records is the invariant checked in
	!      ionization_row_layout: a kinetic row is measured against an
	!      element its species does carry.
	!
	! Row kinds (inventory md/AM_row_layout_inventory_20261006.md section 1):
	!   kinetic, in cm^-3 s^-1, a production-loss balance:
	!     row_species_balance       P - L of one species
	!     row_stage_boundary_flow   the net flow across one stage boundary,
	!                               ionization positive
	!     row_summed_he_i_balance   the summed He I (singlet + 2^3S) balance
	!                               of the atomic He 2^3S systems
	!   constraint, dimensionless, no chemical turnover of its own:
	!     row_identity              x = 0 of an absent element, or of the
	!                               unused X++ of a two-stage element
	!     row_pinned_fraction       x - x_fix of a fraction the flow carries
	!     row_imposed_reservoir     x(H2) - x2*(1-x_ion) + x_O of a lower
	!                               ghost; x_O is zero without OH/H2O
	!   row_element_conservation is the kind of the element rows of the
	!   constrained continuation (constrained_chemical_equilibrium); no
	!   fraction layout carries one.
	!
	! A constraint row replaces the kinetic row of the same index after every
	! term of that row is written (ion_residual_core,
	! impose_transported_ionization_fractions; System_HeH_mol*). It has no
	! chemical turnover: it is judged in its own units against
	! constraint_row_tol, and it never enters the turnover of a group
	! (element_group_turnover_reciprocal with exclude_constraint_rows).
	!
	! Layout of the unknowns (the systems' own headers carry the derivation):
	!   H only (System_H)             1: H I <-> H II
	!   atomic H/He (System_HeH)      1: H I <-> H II, 2: He I <-> He II,
	!                                 3: He II <-> He III; metals from 4
	!   atomic + He 2^3S
	!     (System_HeH_TR)             1: H I <-> H II, 2: summed He I,
	!                                 3: He II <-> He III, 4: He 2^3S;
	!                                 metals from 5
	!   molecular (System_HeH_mol)    1 H+, 2 He+, 3 He++, 4 H2, 5 H2+,
	!                                 6 H3+, 7 HeH+, 8 He 2^3S (triplet);
	!                                 with the oxygen chemistry OH, H2O at
	!                                 8 or 9 and 9 or 10; metals after them
	!   metals                        X0 <-> X+ and X+ <-> X++ of each
	!                                 element in canonical order
	!                                 (species_table), two rows each
	!
	! This module holds no state; every routine is safe inside the OpenMP
	! cell sweep.
	! Reservoir partition: ../EXHALE_v1.00/md/Update_EXHALE_stage3.md section E2.

	use species_table,  only: n_melem, iel_C, iel_O
	use ion_cell_state, only: ion_rates

	implicit none

	! Largest row count of any fraction layout: the molecular block (7), the
	! He 2^3S level, the two oxygen carriers and two rows for each metal.
	integer, parameter :: n_layout_rows_max = 10 + 2*n_melem
	! Acceptance groups: 1 hydrogen, 2 helium, 2 + e metal element e.
	integer, parameter :: group_hydrogen = 1
	integer, parameter :: group_helium   = 2
	integer, parameter :: n_layout_groups = 2 + n_melem

	! Species slots. The first 12 + 3 n_melem are the slots of the
	! constrained continuation (constrained_chemical_equilibrium, is_*;
	! the acceptance_scale test asserts that the two numberings agree);
	! CO is a background species of the cell (ieq_cell%n_co), a row of no
	! layout, and has the last slot so that its stoichiometry is stated.
	integer, parameter :: slot_HI     = 1
	integer, parameter :: slot_HII    = 2
	integer, parameter :: slot_H2     = 3
	integer, parameter :: slot_H2p    = 4
	integer, parameter :: slot_H3p    = 5
	integer, parameter :: slot_HeHp   = 6
	integer, parameter :: slot_HeI    = 7     ! neutral He, ground singlet
	integer, parameter :: slot_HeII   = 8
	integer, parameter :: slot_HeIII  = 9
	integer, parameter :: slot_HeITR  = 10    ! He 2^3S
	integer, parameter :: slot_OH     = 11
	integer, parameter :: slot_H2O    = 12
	! Metal element e: slot_met0 + 3 (e - 1) + (0, 1, 2) = X0, X+, X++.
	integer, parameter :: slot_met0   = 13
	integer, parameter :: slot_CO     = 13 + 3*n_melem
	integer, parameter :: n_layout_species = 13 + 3*n_melem

	! Row kinds.
	integer, parameter :: row_species_balance      = 1
	integer, parameter :: row_stage_boundary_flow  = 2
	integer, parameter :: row_summed_he_i_balance  = 3
	integer, parameter :: row_identity             = 11
	integer, parameter :: row_pinned_fraction      = 12
	integer, parameter :: row_imposed_reservoir    = 13
	integer, parameter :: row_element_conservation = 14
	! Units of a row before any scaling.
	integer, parameter :: row_units_volume_rate = 1    ! cm^-3 s^-1
	integer, parameter :: row_units_fraction    = 2    ! dimensionless

	! THE TOLERANCE OF A CONSTRAINT ROW, in its own (fraction) units: a
	! pinned, imposed or identity row is met when |x - x_fix| <= 1e-6 of
	! the nuclei of its element. The value is the one the molecular
	! acceptance has always applied to these rows (they reached it
	! dimensionless with scale 1); the sweep sets the pinned and reservoir
	! unknowns exactly before it judges a candidate
	! (ionization_equilibrium, impose_transported_fractions), so the value
	! matters in the stored-state measure, where the composition judged is
	! one the transport step may have moved.
	real*8, parameter :: constraint_row_tol = 1.0d-6

	! An element whose nuclei density is at or below this is absent from
	! the cell: its two rows are the identity rows x = 0 (ion_residual_core,
	! metal_rows) and carry no reaction.
	real*8, parameter :: absent_element_density = 1.0d-30

	! THE RESOLUTION OF A FRACTION OF AN ELEMENT. Every unknown of the cell
	! systems is a fraction of the nuclei of an element (a share of a
	! partition whose parts sum to one), and every row reads it through
	! sums over that element (the neutral 1 - sum of the ionized and bound
	! fractions, the electrons it gives). A change of a fraction smaller
	! than the spacing of the floating-point numbers at the element total,
	! epsilon(1) = 2.2e-16, changes none of those sums, so a row whose
	! Newton correction along its own unknown is below it is at its root to
	! the resolution of the composition (normalized_reaction_residual,
	! the resolution floor).
	real*8, parameter :: fraction_resolution = epsilon(1.0d0)

	! THE H3+ ROW AGAINST ITS OWN GROSS RATE. The group norm admits a
	! relative error of a species up to about ieq_res_tol over its share of
	! the group turnover; H3+ is the one species measured to carry a share
	! small enough (down to 2e-11) for that to reach 0.36 in cells where it
	! carries 1 to 100 per cent of the cooling (oxygen_chemistry, r =
	! 1.017-1.040 R_p, and mol_lyman_werner; md/Update_EXHALE_stage3.md
	! section 112.6). Its row is therefore also judged against its own
	! gross rate, |f|/(P + L) <= this (the 0.1 per cent output rule).
	real*8, parameter :: h3p_row_balance_tol = 1.0d-3

	! Status of a row in the turnover aggregation.
	integer, parameter :: turnover_of_group               = 0
	integer, parameter :: turnover_row_without_reaction   = 1  ! gross(i) = 0, G > 0
	integer, parameter :: turnover_group_without_reaction = 2  ! G = 0
	integer, parameter :: turnover_not_finite             = 3  ! gross(i) or G
	integer, parameter :: turnover_constraint_row         = 4

	type :: chemical_row_layout_of_cell
		integer :: n_rows = 0
		! Row of the first metal unknown (n_rows + 1 where none), row of
		! the OH unknown (0 where the oxygen chemistry is off), number of
		! metal elements carried.
		integer :: metal_row_base  = 0
		integer :: oxygen_row_base = 0
		integer :: n_metal = 0
		logical :: has_molecules = .false.
		integer :: kind(n_layout_rows_max)             = 0
		integer :: units(n_layout_rows_max)            = 0
		integer :: species(n_layout_rows_max)          = 0
		integer :: acceptance_group(n_layout_rows_max) = 0
		real*8  :: constraint_tol(n_layout_rows_max)   = 0.0d0
	end type chemical_row_layout_of_cell

	contains

	pure real*8 function helium_closure_allowance(has_triplet, has_molecule) result(allowance)
	! The neutral ground singlet closes He II, He III, the metastable level
	! when present, and HeH+ only when the helium budget is positive.
	! Each summed fraction contributes
	! one fraction_resolution to the representable face allowance.
	logical, intent(in) :: has_triplet, has_molecule
	integer :: n_terms
	n_terms = 2
	if (has_triplet) n_terms = n_terms + 1
	if (has_molecule) n_terms = n_terms + 1
	allowance = dble(n_terms)*fraction_resolution
	end function helium_closure_allowance


	pure real*8 function imposed_stage_sum_rounding_bound(n_terms) result(bound)
	! THE ROUNDING BOUND OF AN IMPOSED STAGE PARTITION, in units of the
	! element's nuclei (u = eps/2 the unit roundoff), for n_terms imposed
	! fractions read from a stored or transported composition whose closing
	! stage c (H I; the He I ground singlet n(He I) - n(He 2^3S) beside a
	! carried 2^3S level; the He I total where the level is not carried) is
	! at the rounding of the element total T.  It bounds two quantities,
	! the excess s - 1 of the imposed fraction sum s and the deficit -c/N of
	! the closing density against the computed total N (closing_stage_
	! density_at_rounding), for the expressions this code evaluates:
	!   * the closing stage was formed by the writer as T minus a weighted
	!     sum of the n_terms other stages, weights 1, 2 or 3 (2 H2, 2 H2+,
	!     3 H3+, 2 H2O): n_terms - 1 additions, at most P = 1 product by a
	!     weight that is not a power of two (3 H3+; a factor 2 is exact) and
	!     one subtraction, so |c_stored - c| <= (n_terms + P) u T';
	!   * the He I ground singlet is a further difference n(He I) - n(He
	!     2^3S) of two stored densities, S = 1 more rounding (and the He I
	!     total already contains the level, so the element total is not
	!     summed over the singlet separately);
	!   * here: n_terms + 1 conversions n_i = f_i n, the total N of the
	!     n_terms + 1 densities with its P weighted product (n_terms + P),
	!     n_terms quotients n_i/N and n_terms - 1 additions of the fraction
	!     sum.
	! T' is the sum of the magnitudes, which exceeds T only by the admitted
	! negative fixed terms, each >= -fraction_resolution, so T' <= T (1 + 6
	! n_terms eps) and the excess enters at second order.  Adding the
	! counts, |s - 1 + c/N| <= (5 n_terms + 2 P + S) u (1 + O(eps)) <=
	! (5 n_terms + 3) u, and the closing density alone carries at most
	! (n_terms + P + S + 2) u.  Both are below (5 n_terms + 1) eps =
	! (10 n_terms + 2) u for n_terms >= 1, which is the bound used: the
	! factor of about two over the count is the convention of
	! ionization_stage_sum_rounding_bound (one eps per rounding).
	! MEASURED: the upper ghost N+1 of the LHS 1140 b He/H 0.1 state written
	! by v1.10 has n(He I) - n(He 2^3S) = -3.8 eps of its helium nuclei,
	! within the 21 eps of n_terms = 4.
	integer, intent(in) :: n_terms
	bound = dble(5*n_terms + 1)*fraction_resolution
	end function imposed_stage_sum_rounding_bound


	pure subroutine closing_stage_density_at_rounding(closure, total,     &
	                                                  n_terms, repaired,  &
	                                                  within, deficit)
	! THE SIGNED CLOSING-STAGE DENSITY OF AN IMPOSED PARTITION, BOUNDED ON
	! THE POSITIVE ELEMENT TOTAL BEFORE IT IS ZEROED.  closure is the
	! closing density as the cell holds it (it may be below zero by
	! rounding), total the positive, finite element total N of the same
	! cell, n_terms the number of imposed terms the closure completes.  A
	! closure >= 0 is returned as it is.  A closure whose relative deficit
	! -closure/total is at most imposed_stage_sum_rounding_bound(n_terms)
	! is a rounding-sized deficit and is returned as zero.  A larger
	! deficit is not rounding: within = .false., repaired = closure, and
	! deficit = -closure/total - that bound (the excess in units of the
	! element nuclei).  The bound is enforced on the density itself,
	! because the fraction-sum check is a different floating-point
	! evaluation and does not by itself bound the density that is zeroed
	! (counterexample: a He sum 21 eps above one beside a ground singlet
	! 21.625 eps below zero).  THE COMPARISON IS DIMENSIONLESS: the
	! relative deficit is formed by one division and compared with the
	! bound, so a subnormal total cannot round the admitted deficit upward
	! (total 6e-310, closure -5e-324: 37 eps, refused); a quotient above
	! the representable range is not finite and refuses.
	real*8,  intent(in)  :: closure, total
	integer, intent(in)  :: n_terms
	real*8,  intent(out) :: repaired, deficit
	logical, intent(out) :: within
	real*8 :: rel, bnd
	bnd      = imposed_stage_sum_rounding_bound(n_terms)
	repaired = closure
	deficit  = 0.0d0
	within   = .true.
	if (closure .ge. 0.0d0) return
	rel = -closure/total
	if (rel .le. bnd) then
		repaired = 0.0d0
		return
	endif
	! rel above the bound, or not finite (overflow of the quotient).
	within  = .false.
	deficit = rel - bnd
	end subroutine closing_stage_density_at_rounding


	pure integer function nuclei_of_species(slot, element_group) result(nu)
	! Nuclei of the element of element_group (1 H, 2 He, 2 + e metal e)
	! carried by one particle of the species in slot. The nuclear
	! stoichiometry matrix, species x element.
	integer, intent(in) :: slot, element_group
	integer :: e, k
	nu = 0
	select case (slot)
	case (slot_HI, slot_HII)
		if (element_group .eq. group_hydrogen) nu = 1
	case (slot_H2, slot_H2p)
		if (element_group .eq. group_hydrogen) nu = 2
	case (slot_H3p)
		if (element_group .eq. group_hydrogen) nu = 3
	case (slot_HeHp)
		if (element_group .eq. group_hydrogen) nu = 1
		if (element_group .eq. group_helium)   nu = 1
	case (slot_HeI, slot_HeII, slot_HeIII, slot_HeITR)
		if (element_group .eq. group_helium) nu = 1
	case (slot_OH)
		if (element_group .eq. group_hydrogen)  nu = 1
		if (element_group .eq. 2 + iel_O)       nu = 1
	case (slot_H2O)
		if (element_group .eq. group_hydrogen)  nu = 2
		if (element_group .eq. 2 + iel_O)       nu = 1
	case (slot_CO)
		if (element_group .eq. 2 + iel_C) nu = 1
		if (element_group .eq. 2 + iel_O) nu = 1
	case default
		if (slot .ge. slot_met0 .and. slot .lt. slot_CO) then
			k = slot - slot_met0
			e = k/3 + 1
			if (element_group .eq. 2 + e) nu = 1
		endif
	end select
	end function nuclei_of_species


	subroutine ionization_row_layout(layout, n_rows, has_he, has_triplet,  &
	                                 has_molecules, has_oxygen_chemistry,  &
	                                 n_metal)
	! THE STATIC LAYOUT of the system a run solves, from its switches: the
	! kind, units, species and acceptance group of every row, all rows
	! kinetic. The constraint rows of one cell are marked afterwards by
	! mark_constraint_rows_of_cell. n_rows is the row count of the caller's
	! system; a count that does not match the switches is a programming
	! defect and stops the run.
	type(chemical_row_layout_of_cell), intent(out) :: layout
	integer, intent(in) :: n_rows, n_metal
	logical, intent(in) :: has_he, has_triplet, has_molecules
	logical, intent(in) :: has_oxygen_chemistry
	integer :: i, e, ix, mbase, n_expected

	layout%n_metal       = n_metal
	layout%has_molecules = has_molecules
	layout%oxygen_row_base = 0

	if (.not. has_he) then
		! System_H: the H I <-> H II flow alone.
		mbase = 2
		call set_row(1, row_stage_boundary_flow, slot_HII, group_hydrogen)
	else if (has_molecules) then
		call set_row(1, row_species_balance, slot_HII,   group_hydrogen)
		call set_row(2, row_species_balance, slot_HeII,  group_helium)
		! He++ is the top stage: its species balance is the flow across
		! He II <-> He III.
		call set_row(3, row_stage_boundary_flow, slot_HeIII, group_helium)
		call set_row(4, row_species_balance, slot_H2,    group_hydrogen)
		call set_row(5, row_species_balance, slot_H2p,   group_hydrogen)
		call set_row(6, row_species_balance, slot_H3p,   group_hydrogen)
		call set_row(7, row_species_balance, slot_HeHp,  group_hydrogen)
		mbase = 8
		if (has_triplet) then
			call set_row(8, row_species_balance, slot_HeITR, group_helium)
			mbase = 9
		endif
		if (has_oxygen_chemistry) then
			layout%oxygen_row_base = mbase
			call set_row(mbase,   row_species_balance, slot_OH,  2 + iel_O)
			call set_row(mbase+1, row_species_balance, slot_H2O, 2 + iel_O)
			mbase = mbase + 2
		endif
	else
		call set_row(1, row_stage_boundary_flow, slot_HII, group_hydrogen)
		if (has_triplet) then
			call set_row(2, row_summed_he_i_balance, slot_HeI, group_helium)
			call set_row(3, row_stage_boundary_flow, slot_HeIII, group_helium)
			call set_row(4, row_species_balance, slot_HeITR, group_helium)
			mbase = 5
		else
			call set_row(2, row_stage_boundary_flow, slot_HeII, group_helium)
			call set_row(3, row_stage_boundary_flow, slot_HeIII, group_helium)
			mbase = 4
		endif
	endif

	n_expected = mbase - 1 + 2*n_metal
	if (n_rows .ne. n_expected .or. n_rows .gt. n_layout_rows_max .or.     &
	    (.not. has_he .and. n_metal .gt. 0)) then
		write(*,'(A,I0,A,I0,A)') ' (chemical_row_layout) the caller''s '// &
		     'system has ', n_rows, ' rows; its switches give ',          &
		     n_expected, '. Aborting.'
		error stop 1
	endif

	do e = 1, n_metal
		ix = mbase + 2*(e-1)
		call set_row(ix,   row_stage_boundary_flow, slot_met0 + 3*(e-1) + 1, &
		             2 + e)
		call set_row(ix+1, row_stage_boundary_flow, slot_met0 + 3*(e-1) + 2, &
		             2 + e)
	enddo
	layout%metal_row_base = mbase
	layout%n_rows         = n_rows

	! The one tie between the stoichiometry and the acceptance grouping: a
	! kinetic row is measured against an element its species carries.
	do i = 1, n_rows
		if (nuclei_of_species(layout%species(i),                            &
		                      layout%acceptance_group(i)) .le. 0) then
			write(*,'(A,I0,A)') ' (chemical_row_layout) row ', i,          &
			     ' is grouped with an element its species does not '//   &
			     'carry. Aborting.'
			error stop 1
		endif
	enddo

	contains

	subroutine set_row(irow, row_kind, slot, grp)
	integer, intent(in) :: irow, row_kind, slot, grp
	layout%kind(irow)             = row_kind
	layout%units(irow)            = row_units_volume_rate
	layout%species(irow)          = slot
	layout%acceptance_group(irow) = grp
	layout%constraint_tol(irow)   = 0.0d0
	end subroutine set_row

	end subroutine ionization_row_layout


	subroutine mark_constraint_rows_of_cell(layout, cell, metal_ntot,      &
	                                        metal_top)
	! THE CONSTRAINT ROWS OF ONE CELL, from the same flags and element
	! densities the systems read when they replace a balance row:
	!   absent metal element (n_X <= absent_element_density): both rows
	!     identity; two-stage element (metal_top < 2): its X++ row identity
	!     (ion_residual_core, metal_rows);
	!   lower-ghost reservoir partition: the H2 row imposed; transported H2
	!     carrier: the H2 row pinned (System_HeH_mol*; the reservoir is
	!     written last there and wins);
	!   transported ionization stages and He 2^3S level: rows 1, 2, 3 and
	!     the level's row pinned (impose_transported_ionization_fractions);
	!   transported oxygen carriers: the OH and H2O rows pinned
	!     (System_HeH_mol_metals).
	! metal_ntot and metal_top are the cell's (System_HeH_metals met_ntot,
	! met_top), read only where the layout carries metals.
	type(chemical_row_layout_of_cell), intent(inout) :: layout
	type(ion_rates), intent(in) :: cell
	real*8,  intent(in) :: metal_ntot(*)
	integer, intent(in) :: metal_top(*)
	integer :: e, ix, iox

	do e = 1, layout%n_metal
		ix = layout%metal_row_base + 2*(e-1)
		if (metal_ntot(e) .le. absent_element_density) then
			call set_constraint(ix,   row_identity)
			call set_constraint(ix+1, row_identity)
		else if (metal_top(e) .lt. 2) then
			call set_constraint(ix+1, row_identity)
		endif
	enddo

	if (layout%has_molecules .and. layout%n_rows .ge. 4) then
		if (cell%x_h2_neutral_partition_fixed) then
			call set_constraint(4, row_imposed_reservoir)
		else if (cell%x_h2_fixed) then
			call set_constraint(4, row_pinned_fraction)
		endif
	endif

	if (cell%x_hp_fixed .and. layout%n_rows .ge. 1)                         &
		call set_constraint(1, row_pinned_fraction)
	if (cell%x_heii_fixed .and. layout%n_rows .ge. 2)                       &
		call set_constraint(2, row_pinned_fraction)
	if (cell%x_heiii_fixed .and. layout%n_rows .ge. 3)                      &
		call set_constraint(3, row_pinned_fraction)
	if (cell%x_hetr_fixed .and. cell%x_hetr_row .gt. 0 .and.                &
	    cell%x_hetr_row .le. layout%n_rows)                                 &
		call set_constraint(cell%x_hetr_row, row_pinned_fraction)

	iox = layout%oxygen_row_base
	if (cell%x_ox_fixed .and. iox .gt. 0) then
		call set_constraint(iox,   row_pinned_fraction)
		call set_constraint(iox+1, row_pinned_fraction)
	endif

	contains

	subroutine set_constraint(irow, row_kind)
	integer, intent(in) :: irow, row_kind
	layout%kind(irow)           = row_kind
	layout%units(irow)          = row_units_fraction
	layout%constraint_tol(irow) = constraint_row_tol
	end subroutine set_constraint

	end subroutine mark_constraint_rows_of_cell


	pure logical function is_constraint_row(row_kind)
	integer, intent(in) :: row_kind
	is_constraint_row = (row_kind .eq. row_identity .or.                    &
	                     row_kind .eq. row_pinned_fraction .or.             &
	                     row_kind .eq. row_imposed_reservoir .or.           &
	                     row_kind .eq. row_element_conservation)
	end function is_constraint_row


	pure logical function finite_number(v)
	! Whether v is an ordinary real number: neither a NaN nor an infinity,
	! read off the exponent field (all ones for both) so that the test
	! cannot be folded away by an optimizer that assumes finite arithmetic.
	real*8, intent(in) :: v
	integer*8 :: bits
	bits = transfer(v, bits)
	finite_number = (iand(ishft(bits, -52), int(z'7FF', 8))                 &
	                 .ne. int(z'7FF', 8))
	end function finite_number


	subroutine element_group_turnover_reciprocal(layout, gross, inv,       &
	                                             status,                   &
	                                             exclude_constraint_rows,  &
	                                             turnover)
	! THE CHEMICAL TURNOVER OF EACH ROW'S ACCEPTANCE GROUP, G_E = sum of
	! the gross rates P + L of the rows of group E, summed in row order,
	! and its reciprocal for every row:
	!   status turnover_of_group               inv = 1/G_E (gross(i) > 0)
	!          turnover_row_without_reaction   inv = 1/G_E; the row itself
	!                                          carries no reaction
	!          turnover_group_without_reaction G_E = 0, inv = 0
	!          turnover_not_finite             gross(i) or G_E not a finite
	!                                          number, inv = 0
	!          turnover_constraint_row         a constraint row, which has
	!                                          no turnover, inv = 0 (only
	!                                          with exclude_constraint_rows)
	! With exclude_constraint_rows false every row counts in its group,
	! constraint rows included. With it true a constraint row neither
	! contributes nor is given a turnover. Both the molecular solve scale
	! (System_HeH_mol, mol_row_scale) and the acceptance
	! (ionization_equilibrium, normalized_reaction_residual) exclude them:
	! the kinetic balance a constraint replaced is not an equation of the
	! cell, and its gross rate is not part of the chemistry the other rows
	! are measured against. What a caller does where status is not
	! turnover_of_group is the caller's decision. turnover (optional)
	! returns G_E itself (0 for a constraint row when excluded).
	type(chemical_row_layout_of_cell), intent(in) :: layout
	real*8,  intent(in)  :: gross(*)
	real*8,  intent(out) :: inv(*)
	integer, intent(out) :: status(*)
	logical, intent(in)  :: exclude_constraint_rows
	real*8,  intent(out), optional :: turnover(*)
	real*8  :: g_group(n_layout_groups), g
	integer :: i, n

	n = layout%n_rows
	g_group(:) = 0.0d0
	do i = 1, n
		if (exclude_constraint_rows .and.                                   &
		    is_constraint_row(layout%kind(i))) cycle
		g_group(layout%acceptance_group(i)) =                               &
			g_group(layout%acceptance_group(i)) + gross(i)
	enddo
	do i = 1, n
		if (exclude_constraint_rows .and.                                   &
		    is_constraint_row(layout%kind(i))) then
			status(i) = turnover_constraint_row
			inv(i)    = 0.0d0
			if (present(turnover)) turnover(i) = 0.0d0
			cycle
		endif
		g = g_group(layout%acceptance_group(i))
		if (present(turnover)) turnover(i) = g
		if (.not. (finite_number(gross(i)) .and. finite_number(g))) then
			status(i) = turnover_not_finite
			inv(i)    = 0.0d0
		else if (.not. (g .gt. 0.0d0)) then
			status(i) = turnover_group_without_reaction
			inv(i)    = 0.0d0
		else if (.not. (gross(i) .gt. 0.0d0)) then
			status(i) = turnover_row_without_reaction
			inv(i)    = 1.0d0/g
		else
			status(i) = turnover_of_group
			inv(i)    = 1.0d0/g
		endif
	enddo
	end subroutine element_group_turnover_reciprocal

	end module chemical_row_layout

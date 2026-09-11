      module element_inventory
      ! ONE stoichiometric map of the hydrogen, helium, oxygen and carbon
      ! nuclei of a cell, and the three constraints a composition state is
      ! judged by.  Every multiplicity is read from the species table
      ! (bsp_nH, bsp_nHe, bsp_nO, bsp_nC and the metal stage tables of
      ! species_table), so this module carries no species list of its own
      ! and adding a carrier is adding a row there.
      !
      ! WHY ONE MAP AND NOT ONE BOUND FOR EACH CARRIER.  A carrier holds
      ! nuclei of an element that other species of the same cell also hold,
      ! so the admissible set is a simplex and not a product of intervals.
      ! With the available hydrogen normalized to one, n(H2) <= 0.5 and
      ! n(H+) <= 1 are each satisfied by a state holding two hydrogen
      ! nuclei where the cell has one.  The constraint that state breaks is
      !
      !     2 n(H2) + n(H+) + n(OH) + 2 n(H2O)  <=  available hydrogen,
      !
      ! with n(OH) + n(H2O) + n(CO) <= free oxygen and n(CO) <= free carbon
      ! beside it, and the helium nucleus of HeH+ inside the helium total.
      ! Those are the constraints this module evaluates, and it reports the
      ! MAGNITUDE by which a state breaks them: a count of offending cells
      ! bounds neither the size of the breach nor which cells carry it.
      !
      ! THE THREE CONTEXTS OF THE SAME CONSTRAINT.
      !
      !  (a) THE FULL INVENTORY of a state: every nucleus of every element,
      !      whatever species carries it.  n_total below, and the identity
      !      n_total = n_carried + n_locked + n_reserved + n_closure holds
      !      for each element by construction of the groups.
      !
      !  (b) THE CONDITIONAL BUDGET of a calculation that holds some
      !      densities fixed.  A carrier transport step freezes the
      !      molecular ions and, unless the proton is transported, H II,
      !      and it rescales them rather than re-solving them, so the
      !      nuclei sitting in them are not the carriers' to take
      !      (carrier_budget).  The helium hydride is not even rescaled,
      !      because rescaling it would move a helium nucleus as well, so
      !      its hydrogen leaves the hydrogen the step can move at all
      !      (element_free).  Oxygen and carbon have a second, LOOSER
      !      budget in the code today: the box sides of their carriers
      !      spend the element total, while the chemistry rows close free
      !      atomic oxygen and carbon against the total less the stages
      !      above neutral (budget_less_reserved).  The two are reported
      !      separately rather than averaged into one number.
      !
      !  (c) AFTER CONSISTENT REDISTRIBUTION.  The remainder each element
      !      leaves for its closure species, remainder = element_free -
      !      n_carried.  A negative remainder is a state whose carriers
      !      hold more nuclei than the cell has; a redistribution that
      !      clamps it at zero keeps the carriers and empties the closure
      !      species, so the element total of the written state EXCEEDS the
      !      one it was formed from by the size of that remainder.
      !
      ! EVALUATED AT A STATE THAT CAME FROM ITS OWN ELEMENT CLOSURE the
      ! hydrogen constraint reduces to n(H I) >= 0 and is an identity, up to
      ! the rounding of two different summation orders.  It has content for
      ! a CANDIDATE state: a trial of the coupled solve, whose budget is
      ! frozen at the outer iterate while its carriers move, and the state a
      ! transport step writes, whose budget was formed at the state the step
      ! entered.  inventory_feasible_slack is a representability bound and
      ! not a tolerance: eight units in the last place cover the rounding of
      ! the two sums.

      use global_parameters, only: N, Ng, n_species, n0, he_diffusion
      use species_table, only: n_bsp, bsp_fsp, bsp_nH, bsp_nHe, bsp_nO,   &
                               bsp_nC, bsp_is_excited_level,              &
                               isp_HII, isp_H2, isp_H2p, isp_H3p,         &
                               isp_HeHp, isp_OH, isp_H2O, isp_CO,         &
                               n_mion, mion_fsp, mion_elem, mion_stage,   &
                               melem_i0, melem_top, iel_O, iel_C
      ! THE GATE ON AN ELEMENT RATIO, one definition for the whole code
      ! (element_census, where the same number bounds the census's own
      ! ratio check).  The carrier breach of this module compares two
      ! expressions of the SAME nuclei of ONE cell, so its bound is
      ! representability, inventory_feasible_slack below; a ratio of two
      ! elements is an exact bookkeeping identity carried across every
      ! operator of a run and is bounded by the accumulated round-off of the
      ! repartitions instead.
      use element_census, only: element_ratio_gate

      implicit none
      private

      integer, parameter :: dp = kind(1.0d0)

      ! The elements a molecular carrier can hold: hydrogen, helium, oxygen
      ! and carbon, the elements of H2, H2+, H3+, HeH+, OH, H2O and CO.
      ! Every other element of species_table sits in one stage family that
      ! no carrier touches, and its inventory is the element census.
      integer, parameter, public :: n_inv_element = 4
      integer, parameter, public :: ien_H  = 1
      integer, parameter, public :: ien_He = 2
      integer, parameter, public :: ien_O  = 3
      integer, parameter, public :: ien_C  = 4
      character(len=2), parameter, public ::                              &
         inv_element_name(n_inv_element) = [ 'H ', 'He', 'O ', 'C ' ]

      ! Eight units in the last place of a double, which is what the
      ! rounding of the inventory sum and of the budget sum can put between
      ! two expressions of the same nuclei.  A state on the face of the
      ! simplex is on the face; anything above this is a breach with a size.
      real(dp), parameter, public :: inventory_feasible_slack              &
                                   = 8.0d0*epsilon(1.0d0)

      type, public :: element_inventory_cell
         ! (a) Every nucleus of the element in the cell [cm^-3].
         real(dp) :: n_total(n_inv_element)    = 0.0d0
         ! Nuclei held by the transported carriers of this calculation.
         real(dp) :: n_carried(n_inv_element)  = 0.0d0
         ! Nuclei a carrier step cannot move at all: the hydrogen and the
         ! helium of HeH+, one molecule carrying a nucleus of each.
         real(dp) :: n_locked(n_inv_element)   = 0.0d0
         ! Nuclei the conditional budget withholds from the carriers: the
         ! stages a step holds fixed and rescales rather than re-solves.
         real(dp) :: n_reserved(n_inv_element) = 0.0d0
         ! Nuclei in the species that take the remainder of the element.
         real(dp) :: n_closure(n_inv_element)  = 0.0d0
         ! (b) The element a step can move at all, and the part of it the
         ! carriers may spend.  element_free is the carrier operator's
         ! nH_free, nO_free and nC_free; carrier_budget is what its box
         ! sides spend, hydrogen less the frozen stages and oxygen and
         ! carbon whole.
         real(dp) :: element_free(n_inv_element)   = 0.0d0
         real(dp) :: carrier_budget(n_inv_element) = 0.0d0
         ! The same budget with the frozen stages of every element removed,
         ! which is what the chemistry rows close free atomic oxygen and
         ! carbon against.
         real(dp) :: budget_less_reserved(n_inv_element) = 0.0d0
         ! (c) What the closure species are left, element_free - n_carried.
         real(dp) :: remainder(n_inv_element) = 0.0d0
         ! The breach of the shared constraint of each element, relative to
         ! the budget it is measured against, and zero where there is none.
         real(dp) :: violation(n_inv_element) = 0.0d0
         real(dp) :: violation_less_reserved(n_inv_element) = 0.0d0
         ! How much of the budget the carriers occupy, n_carried divided by
         ! carrier_budget.  At or below one is inside the constraint, and
         ! the number says how tight the state runs.
         real(dp) :: occupancy(n_inv_element) = 0.0d0
         ! The density [cm^-3] of every species this map counts: the base
         ! species except the excited levels, whose population is already
         ! inside their parent, and the oxygen and carbon stages.  It is
         ! here so that a caller asking about one carrier and a caller
         ! asking about a shared element read one state.  Zero for a
         ! species the map does not count.
         real(dp) :: n_of_species(n_species) = 0.0d0
      end type element_inventory_cell

      public :: element_inventory_of_cell, element_inventory_of_candidate
      public :: carrier_species_that_hold, nuclei_per_particle
      public :: element_box_side, carrier_box_breach
      public :: carrier_element_totals, carrier_hydrogen_budget
      public :: inventory_is_feasible, inventory_worst_violation
      public :: element_inventory_report, element_inventory_report_on
      public :: inventory_judges_the_shared_constraint
      public :: helium_to_hydrogen, ratio_departure_from_the_base
      public :: element_ratio_may_depart
      public :: element_constraint_derivatives

      contains

      ! ------------------------------------------------------------- !

      ! NUCLEI OF ELEMENT ie CARRIED BY ONE PARTICLE OF f_sp SPECIES isp.
      ! The one place a stoichiometric coefficient of this map is read, and
      ! it is read from species_table: two hydrogen nuclei in H2 and H2+,
      ! three in H3+, one hydrogen and one helium in HeH+, one hydrogen and
      ! one oxygen in OH, two and one in H2O, one carbon and one oxygen in
      ! CO, and one nucleus of its own element in each metal ion stage.
      integer function nuclei_per_particle(isp, ie) result(nu)
      integer, intent(in) :: isp, ie
      integer :: ib, im
      nu = 0
      do ib = 1, n_bsp
         if (bsp_fsp(ib) .ne. isp) cycle
         if (bsp_is_excited_level(ib)) cycle
         select case (ie)
         case (ien_H);  nu = bsp_nH(ib)
         case (ien_He); nu = bsp_nHe(ib)
         case (ien_O);  nu = bsp_nO(ib)
         case (ien_C);  nu = bsp_nC(ib)
         end select
         return
      enddo
      do im = 1, n_mion
         if (mion_fsp(im) .ne. isp) cycle
         if (ie .eq. ien_O .and. mion_elem(im) .eq. iel_O) nu = 1
         if (ie .eq. ien_C .and. mion_elem(im) .eq. iel_C) nu = 1
         return
      enddo
      end function nuclei_per_particle

      ! ------------------------------------------------------------- !

      ! THE LARGEST DENSITY SPECIES isp CAN REACH ON A BUDGET OF ELEMENT ie
      ! [cm^-3]: the budget divided by the nuclei of that element one
      ! particle holds.  It is one side of the box a step's unknowns live
      ! in, and it is stoichiometry alone.  Which element's budget a carrier
      ! is bounded against is the caller's statement, and a side formed from
      ! one element does NOT bound the nuclei of any other element the same
      ! species holds: a state can satisfy every side and still break a
      ! shared constraint, which is the magnitude
      ! element_inventory_of_cell returns.
      double precision function element_box_side(isp, ie, budget)         &
                               result(nmax)
      integer,  intent(in) :: isp, ie
      real(dp), intent(in) :: budget
      integer :: nu
      nu = nuclei_per_particle(isp, ie)
      if (nu .le. 0) then
         write(*,'(a,i0,a,a)') ' (element_box_side) f_sp species ', isp,  &
            ' holds no nuclei of element ', trim(inv_element_name(ie))
         error stop 1
      endif
      nmax = budget/dble(nu)
      end function element_box_side

      ! ------------------------------------------------------------- !

      ! THE SPECIES A CARRIER OPERATOR CAN TRANSPORT, as a mask over the
      ! f_sp columns.  The set is stoichiometric: H2, OH, H2O and CO hold
      ! every hydrogen, oxygen and carbon nucleus that is neither atomic nor
      ! in a molecular ion, and the proton joins them when the ionization
      ! state is transported.  A calculation that solves fewer of them
      ! passes its own mask; this is the widest one.
      subroutine carrier_species_that_hold(proton_carried, carried)
      logical, intent(in)  :: proton_carried
      logical, intent(out) :: carried(n_species)
      carried = .false.
      carried(isp_H2)  = .true.
      carried(isp_OH)  = .true.
      carried(isp_H2O) = .true.
      carried(isp_CO)  = .true.
      if (proton_carried) carried(isp_HII) = .true.
      end subroutine carrier_species_that_hold

      ! ------------------------------------------------------------- !

      ! THE INVENTORY AND THE CONSTRAINTS OF ONE CELL.
      !
      ! f_cell holds the species fractions of the cell and nd the density
      ! they are fractions of [cm^-3 of m_H], so a species density is
      ! f_cell(isp)*nd, the same nd = rho*n0 every composition operator
      ! uses.  A caller holding DENSITIES passes them as f_cell with
      ! nd = 1, and the products are then exact.
      !
      ! carried marks the species this calculation transports.  Which stages
      ! the conditional budget withholds follows from it: the proton is
      ! either a carrier or a frozen stage, never both.
      subroutine element_inventory_of_cell(nd, f_cell, carried, inv)
      real(dp), intent(in)  :: nd
      real(dp), intent(in)  :: f_cell(n_species)
      logical,  intent(in)  :: carried(n_species)
      type(element_inventory_cell), intent(out) :: inv
      real(dp) :: c, free(n_inv_element)
      integer  :: ib, im, ie, isp

      inv = element_inventory_cell()

      ! ---- (a) the full inventory, by the group each species is in -----
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         isp = bsp_fsp(ib)
         c   = f_cell(isp)*nd
         call add_nuclei(inv, isp, c, carried,                            &
                         [ bsp_nH(ib), bsp_nHe(ib), bsp_nO(ib),           &
                           bsp_nC(ib) ])
      enddo
      ! The metal ion stages hold one nucleus of their element each; only
      ! oxygen and carbon are elements of this map.
      do im = 1, n_mion
         if (mion_elem(im) .ne. iel_O .and. mion_elem(im) .ne. iel_C) cycle
         isp = mion_fsp(im)
         c   = f_cell(isp)*nd
         if (mion_elem(im) .eq. iel_O) then
            call add_nuclei(inv, isp, c, carried, [ 0, 0, 1, 0 ])
         else
            call add_nuclei(inv, isp, c, carried, [ 0, 0, 0, 1 ])
         endif
      enddo

      ! ---- (b) the budgets and (c) the breach --------------------------
      ! The element a step can move at all.  Hydrogen loses the nucleus
      ! locked in HeH+; oxygen and carbon lose nothing, because every stage
      ! of them is rescaled onto the remainder; helium loses nothing here
      ! because no carrier holds it.
      do ie = 1, n_inv_element
         free(ie) = max(inv%n_total(ie) - inv%n_locked(ie), 0.0d0)
      enddo
      free(ien_He) = inv%n_total(ien_He)
      call budgets_and_breach(inv, free)

      end subroutine element_inventory_of_cell

      ! ------------------------------------------------------------- !

      ! THE CONSTRAINTS OF A CANDIDATE STATE AGAINST A FROZEN BUDGET.
      !
      ! This is context (b) of the header written down: a solve holds some
      ! densities fixed, so the element densities its constraint is measured
      ! against belong to the state the budget was frozen at, while the
      ! carriers belong to the candidate.  element_free_frozen carries the
      ! frozen densities [cm^-3] in the order of the element indices, and it
      ! is the carrier operator's nH_free, nO_free and nC_free.
      !
      ! The stages the budget withholds are read from the candidate, which
      ! is the same number the frozen state carries whenever a candidate
      ! moves the carriers alone, and it is the number a caller that moved
      ! them as well has to be measured by.
      subroutine element_inventory_of_candidate(nd, f_cell, carried,      &
                                        element_free_frozen, inv)
      real(dp), intent(in) :: nd
      real(dp), intent(in) :: f_cell(n_species)
      logical,  intent(in) :: carried(n_species)
      real(dp), intent(in) :: element_free_frozen(n_inv_element)
      type(element_inventory_cell), intent(out) :: inv
      call element_inventory_of_cell(nd, f_cell, carried, inv)
      call budgets_and_breach(inv, element_free_frozen)
      end subroutine element_inventory_of_candidate

      ! ------------------------------------------------------------- !

      ! The budgets of context (b) and the breach of context (c), from an
      ! inventory whose groups are already counted and an element density
      ! the caller states.  ONE place, so that the constraint a candidate is
      ! measured by and the constraint a state states about itself cannot
      ! drift apart: they differ only in which element density is handed in.
      subroutine budgets_and_breach(inv, free)
      type(element_inventory_cell), intent(inout) :: inv
      real(dp), intent(in) :: free(n_inv_element)
      real(dp) :: occ_less_reserved
      integer  :: ie
      inv%element_free = free
      ! What the box sides spend: for hydrogen the element a step can move
      ! less the stages it freezes, for oxygen and carbon the element total.
      inv%carrier_budget         = free
      inv%carrier_budget(ien_H)  = max(free(ien_H)                        &
                                       - inv%n_reserved(ien_H), 0.0d0)
      inv%carrier_budget(ien_He) = 0.0d0
      do ie = 1, n_inv_element
         inv%budget_less_reserved(ie) = max(free(ie)                      &
                                            - inv%n_reserved(ie), 0.0d0)
      enddo
      inv%budget_less_reserved(ien_He) = 0.0d0
      do ie = 1, n_inv_element
         inv%remainder(ie) = inv%element_free(ie) - inv%n_carried(ie)
         call relative_breach(inv%n_carried(ie), inv%carrier_budget(ie),  &
                              inv%violation(ie), inv%occupancy(ie))
         call relative_breach(inv%n_carried(ie),                          &
                              inv%budget_less_reserved(ie),              &
                              inv%violation_less_reserved(ie),           &
                              occ_less_reserved)
      enddo
      end subroutine budgets_and_breach

      ! ------------------------------------------------------------- !

      ! One species' nuclei added to the totals and to the group it belongs
      ! to.  THE GROUPS ARE THE CONDITIONAL BUDGET'S: what decides
      ! feasibility is which nuclei a carrier may take.
      !   * locked: HeH+, whose hydrogen cannot be rescaled without moving
      !     its helium;
      !   * carried: the species this calculation transports;
      !   * reserved: the stages held fixed and rescaled, which are H II
      !     when it is not transported, H2+ and H3+ always, and the oxygen
      !     and carbon stages above neutral that the chemistry rows subtract
      !     before closing free atomic O and C;
      !   * closure: what takes the remainder, which is H I, the helium
      !     stages, and the neutral oxygen and carbon stages.
      subroutine add_nuclei(inv, isp, c, carried, nu)
      type(element_inventory_cell), intent(inout) :: inv
      integer,  intent(in) :: isp
      real(dp), intent(in) :: c
      logical,  intent(in) :: carried(n_species)
      integer,  intent(in) :: nu(n_inv_element)
      real(dp) :: q(n_inv_element)
      integer  :: ie
      do ie = 1, n_inv_element
         q(ie) = dble(nu(ie))*c
      enddo
      inv%n_of_species(isp) = c
      inv%n_total = inv%n_total + q
      if (isp .eq. isp_HeHp) then
         inv%n_locked = inv%n_locked + q
      else if (carried(isp)) then
         inv%n_carried = inv%n_carried + q
      else if (isp .eq. isp_HII .or. isp .eq. isp_H2p .or.                &
               isp .eq. isp_H3p .or. stage_above_neutral(isp)) then
         inv%n_reserved = inv%n_reserved + q
      else
         inv%n_closure = inv%n_closure + q
      endif
      end subroutine add_nuclei

      ! ------------------------------------------------------------- !

      ! .true. for an oxygen or carbon stage above the neutral one.
      logical function stage_above_neutral(isp) result(above)
      integer, intent(in) :: isp
      integer :: im
      above = .false.
      do im = 1, n_mion
         if (mion_fsp(im) .ne. isp) cycle
         if (mion_elem(im) .ne. iel_O .and. mion_elem(im) .ne. iel_C) exit
         above = (mion_stage(im) .gt. 0)
         exit
      enddo
      end function stage_above_neutral

      ! ------------------------------------------------------------- !

      ! The breach of one constraint as a fraction of the budget it is
      ! measured against, and the occupancy of that budget.  With no nuclei
      ! of the element in the cell, any demand for them is infeasible and
      ! has no finite fraction, so the breach is huge(1) and says so.
      subroutine relative_breach(demand, budget, breach, occupancy)
      real(dp), intent(in)  :: demand, budget
      real(dp), intent(out) :: breach, occupancy
      if (budget .gt. 0.0d0) then
         breach    = max(demand - budget, 0.0d0)/budget
         occupancy = demand/budget
      else if (demand .gt. 0.0d0) then
         breach    = huge(1.0d0)
         occupancy = huge(1.0d0)
      else
         breach    = 0.0d0
         occupancy = 0.0d0
      endif
      end subroutine relative_breach

      ! ------------------------------------------------------------- !

      ! THE ELEMENT DENSITIES THE CARRIER OPERATOR SPENDS, cell by cell
      ! [cm^-3]: the hydrogen a step can move, which is the element total
      ! less the nucleus locked in HeH+, and the oxygen and carbon totals of
      ! the cell whatever species carries them.  The multiplicities are the
      ! species table's, so these are the map's element_free evaluated over
      ! the grid and the two entry points return the same numbers.
      subroutine carrier_element_totals(nd, f_sp, nH_free, nO_free,       &
                                        nC_free)
      real(dp), dimension(1-Ng:N+Ng),           intent(in)  :: nd
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nH_free
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nO_free
      real(dp), dimension(1-Ng:N+Ng),           intent(out) :: nC_free
      integer :: ib, i0, k
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
      end subroutine carrier_element_totals

      ! ------------------------------------------------------------- !

      ! THE HYDROGEN NUCLEI THE TRANSPORTED CARRIERS OF A CELL MAY HOLD
      ! [cm^-3]: the hydrogen a step can move, less the nuclei sitting in
      ! the stages it holds FROZEN and rescales rather than re-solves.
      !
      ! It is NOT the hydrogen element total.  A write-back that does not
      ! re-solve H I, H II, H2+ and H3+ but rescales them in proportion to
      ! what they already held leaves those stages nothing when the carriers
      ! spend the whole element, the free electron density with them, and
      ! the ionization sweep that follows a cell it can find no root in.
      !
      ! WHICH STAGES ARE FROZEN DEPENDS ON THE CALCULATION.  With the
      ! ionization state transported the proton is a carrier, solved by the
      ! same system as H2, so its nuclei are on the carriers' side of the
      ! budget and the caller accounts for them there; subtracting them here
      ! as well would charge one nucleus twice and cap the carriers below
      ! the hydrogen the cell has.
      double precision function carrier_hydrogen_budget(nH_free, n_hii,   &
                               n_h2p, n_h3p, proton_carried) result(nH)
      real(dp), intent(in) :: nH_free, n_hii, n_h2p, n_h3p
      logical,  intent(in) :: proton_carried
      if (proton_carried) then
         nH = max(nH_free - 2.0d0*n_h2p - 3.0d0*n_h3p, 0.0d0)
      else
         nH = max(nH_free - n_hii - 2.0d0*n_h2p                          &
                          - 3.0d0*n_h3p, 0.0d0)
      endif
      end function carrier_hydrogen_budget

      ! ------------------------------------------------------------- !

      ! WHETHER THE SHARED CONSTRAINT IS WHAT A STATE IS JUDGED BY.
      ! EXHALE_INVENTORY_SHARED=0 judges by the separate box sides of the
      ! carriers instead, which is what the shared constraint is measured
      ! against.  It is a measurement switch of the classification and not a
      ! physics option: no production path reads it, and the constraint is a
      ! property of the state either way.
      logical function inventory_judges_the_shared_constraint() result(on)
      character(len=8) :: env
      call get_environment_variable('EXHALE_INVENTORY_SHARED', env)
      on = (trim(env) .ne. '0')
      end function inventory_judges_the_shared_constraint

      ! ------------------------------------------------------------- !

      ! THE SHARED CONSTRAINT OF ONE CELL AS A SIGNED QUANTITY, AND ITS
      ! DERIVATIVES, so that a step control can carry the constraint in its
      ! own model instead of discovering it by refusing an evaluated state.
      !
      !     c(ie) = (nuclei of ie the carriers of the cell hold)
      !             - (budget of ie)                              <= 0
      !
      ! in cm^-3.  c <= 0 is the feasible half-space and the relative breach
      ! budgets_and_breach returns is c/budget wherever the budget is
      ! positive.  carrier_budget_frozen is the budget the calculation holds
      ! fixed for its step -- the conditional budget of context (b): for
      ! hydrogen the nuclei left over after the stages the step does not
      ! re-solve (carrier_hydrogen_budget), for oxygen and carbon the
      ! element densities themselves.
      !
      ! THE THREE DERIVATIVES ARE THE THREE WAYS A CELL'S UNKNOWNS ENTER IT.
      !
      !   dc_dn(isp,ie)  the density of species isp [cm^-3].  A carried
      !     species adds nuclei_per_particle(isp,ie) of them to the demand
      !     and nothing else does, the budget being frozen data of the step.
      !     A step whose unknown is ln n chains this with dn/d ln n = n.
      !
      !   dc_dnd(ie)  the cell's own particle density [dimensionless], at
      !     fixed mass fractions and with the frozen budget read as a fixed
      !     number of nuclei per unit mass of the cell, which is what the
      !     budget of a step is: the composition it was formed from does not
      !     move with the density.  Every term of c is then homogeneous of
      !     degree one in nd, so dc/dnd = c/nd exactly, and the constraint
      !     surface c = 0 CONTAINS the ray through the origin -- the density
      !     direction is tangent to it, and an active constraint therefore
      !     does not push the cell's density.
      !
      !   dc_dbudget(ie)  the budget itself [dimensionless], -1 for every
      !     element a carrier holds and 0 for helium.  A step that carries
      !     an element as an unknown chains this with the derivative of
      !     that element's density with respect to its own unknown.
      !
      ! Helium is held by no carrier, so c(ien_He) is identically zero and
      ! so are its derivatives.
      !
      ! The constraint is affine in the species densities and in the frozen
      ! budget and homogeneous of degree one in nd, so a finite difference
      ! of it reproduces these derivatives exactly and a test row needs no
      ! tolerance beyond the rounding of the two sums.
      subroutine element_constraint_derivatives(nd, f_cell, carried,      &
                                        carrier_budget_frozen, c, dc_dn,  &
                                        dc_dnd, dc_dbudget)
      real(dp), intent(in)  :: nd
      real(dp), intent(in)  :: f_cell(n_species)
      logical,  intent(in)  :: carried(n_species)
      real(dp), intent(in)  :: carrier_budget_frozen(n_inv_element)
      real(dp), intent(out) :: c(n_inv_element)
      real(dp), intent(out) :: dc_dn(n_species,n_inv_element)
      real(dp), intent(out) :: dc_dnd(n_inv_element)
      real(dp), intent(out) :: dc_dbudget(n_inv_element)
      type(element_inventory_cell) :: inv
      integer  :: ie, ib, im, isp, nu

      call element_inventory_of_cell(nd, f_cell, carried, inv)
      c = 0.0d0;  dc_dn = 0.0d0;  dc_dnd = 0.0d0;  dc_dbudget = 0.0d0

      do ie = 1, n_inv_element
         if (ie .eq. ien_He) cycle
         c(ie) = inv%n_carried(ie) - carrier_budget_frozen(ie)
         if (nd .gt. 0.0d0) dc_dnd(ie) = c(ie)/nd
         dc_dbudget(ie) = -1.0d0
      enddo

      ! The species columns, from the one stoichiometric source of the map.
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         isp = bsp_fsp(ib)
         if (isp .eq. isp_HeHp) cycle
         if (.not. carried(isp)) cycle
         do ie = 1, n_inv_element
            if (ie .eq. ien_He) cycle
            nu = nuclei_per_particle(isp, ie)
            if (nu .gt. 0) dc_dn(isp,ie) = dc_dn(isp,ie) + dble(nu)
         enddo
      enddo
      do im = 1, n_mion
         if (mion_elem(im) .ne. iel_O .and. mion_elem(im) .ne. iel_C) cycle
         isp = mion_fsp(im)
         if (.not. carried(isp)) cycle
         ie  = ien_O
         if (mion_elem(im) .eq. iel_C) ie = ien_C
         dc_dn(isp,ie) = dc_dn(isp,ie) + 1.0d0
      enddo

      end subroutine element_constraint_derivatives

      ! ------------------------------------------------------------- !

      ! THE BREACH OF ONE CARRIER AGAINST ITS OWN BOX SIDE, as a fraction of
      ! that side.  This is the separate bound each carrier carries today:
      ! H2 and H+ against the available hydrogen, OH and H2O against the
      ! free oxygen, CO against the smaller of the free oxygen and the free
      ! carbon, since it needs one nucleus of each.
      double precision function carrier_box_breach(isp, inv) result(b)
      integer, intent(in) :: isp
      type(element_inventory_cell), intent(in) :: inv
      real(dp) :: side, occ
      b = 0.0d0
      if (isp .eq. isp_H2 .or. isp .eq. isp_HII) then
         side = element_box_side(isp, ien_H, inv%carrier_budget(ien_H))
      else if (isp .eq. isp_OH .or. isp .eq. isp_H2O) then
         side = element_box_side(isp, ien_O, inv%carrier_budget(ien_O))
      else if (isp .eq. isp_CO) then
         side = element_box_side(isp, ien_O,                              &
                   min(inv%carrier_budget(ien_O),                         &
                       inv%carrier_budget(ien_C)))
      else
         return
      endif
      call relative_breach(inv%n_of_species(isp), side, b, occ)
      end function carrier_box_breach

      ! ------------------------------------------------------------- !

      ! WHETHER A STATE IS ADMISSIBLE, by the size of its worst breach and
      ! not by a count of cells.
      logical function inventory_is_feasible(inv) result(ok)
      type(element_inventory_cell), intent(in) :: inv
      integer  :: ie
      real(dp) :: worst
      if (inventory_judges_the_shared_constraint()) then
         ok = (inventory_worst_violation(inv, ie) .le.                    &
               inventory_feasible_slack)
      else
         ! The separate sides, whose conjunction admits a state holding two
         ! hydrogen nuclei where the cell has one.
         worst = 0.0d0
         worst = max(worst, carrier_box_breach(isp_H2,  inv))
         worst = max(worst, carrier_box_breach(isp_HII, inv))
         worst = max(worst, carrier_box_breach(isp_OH,  inv))
         worst = max(worst, carrier_box_breach(isp_H2O, inv))
         worst = max(worst, carrier_box_breach(isp_CO,  inv))
         ok = (worst .le. inventory_feasible_slack)
      endif
      end function inventory_is_feasible

      ! ------------------------------------------------------------- !

      ! THE WORST BREACH OF A CELL and the element that carries it.
      double precision function inventory_worst_violation(inv, ie_worst)  &
                               result(worst)
      type(element_inventory_cell), intent(in) :: inv
      integer, intent(out) :: ie_worst
      integer :: ie
      worst    = 0.0d0
      ie_worst = 0
      do ie = 1, n_inv_element
         if (inv%violation(ie) .gt. worst) then
            worst    = inv%violation(ie)
            ie_worst = ie
         endif
      enddo
      end function inventory_worst_violation

      ! ------------------------------------------------------------- !

      ! THE HELIUM TO HYDROGEN NUCLEUS RATIO OF A CELL, from the same
      ! inventory the carrier budgets come from.  Zero where the cell has no
      ! hydrogen, which is not a state the atmosphere reaches and is
      ! reported rather than divided by.
      double precision function helium_to_hydrogen(inv) result(y)
      type(element_inventory_cell), intent(in) :: inv
      y = 0.0d0
      if (inv%n_total(ien_H) .gt. 0.0d0)                                  &
         y = inv%n_total(ien_He)/inv%n_total(ien_H)
      end function helium_to_hydrogen

      ! ------------------------------------------------------------- !

      ! THE LARGEST DEPARTURE OF AN ELEMENT RATIO FROM THE BASE CELL'S, and
      ! the cell that carries it.  ratio(1) is the base cell's, which is the
      ! reservoir the element operator pins its Dirichlet row to.
      !
      ! WHETHER A DEPARTURE IS A BREACH IS A QUESTION ABOUT THE RUN.  With
      ! no element transport active no operator of this code may move the
      ! ratio of two elements: the ionization sweep repartitions stages, the
      ! carrier step moves nuclei between species of ONE element, and the
      ! write-back restores the element totals it was handed, so a departure
      ! is an element-conservation defect of the state.  With helium
      ! diffusion on the helium mass fraction is a transported unknown and a
      ! gradient IS the solution, so the same number is a measurement.
      subroutine ratio_departure_from_the_base(ratio, nc, worst, jworst)
      real(dp), intent(in)  :: ratio(nc)
      integer,  intent(in)  :: nc
      real(dp), intent(out) :: worst
      integer,  intent(out) :: jworst
      real(dp) :: d
      integer  :: j
      worst  = 0.0d0
      jworst = 0
      if (nc .lt. 2) return
      if (.not. (ratio(1) .gt. 0.0d0)) return
      do j = 2, nc
         d = abs(ratio(j)/ratio(1) - 1.0d0)
         if (d .gt. worst) then
            worst  = d
            jworst = j
         endif
      enddo
      end subroutine ratio_departure_from_the_base

      ! ------------------------------------------------------------- !

      ! WHETHER AN ELEMENT RATIO OF THIS RUN MAY DEPART FROM THE RESERVOIR.
      ! True only while an operator transports an element: with helium
      ! diffusion off, the ratio of helium to hydrogen nuclei is fixed data
      ! of the run and a departure is a defect.
      logical function element_ratio_may_depart() result(may)
      may = he_diffusion
      end function element_ratio_may_depart

      ! ------------------------------------------------------------- !

      logical function element_inventory_report_on() result(on)
      ! EXHALE_INVENTORY_REPORT=1 prints the inventory verdict of a state.
      ! Default off: the report costs one pass over the grid.
      character(len=8) :: env
      call get_environment_variable('EXHALE_INVENTORY_REPORT', env)
      on = (trim(env) .eq. '1')
      end function element_inventory_report_on

      ! ------------------------------------------------------------- !

      ! THE VERDICT OF A STATE, one block of lines under
      ! EXHALE_INVENTORY_REPORT=1: the worst breach of each element with the
      ! cell that carries it, the tightest occupancy of each budget, and
      ! whether the state is admissible.  A breach is a magnitude; the
      ! number of cells that carry one is printed beside it and is not the
      ! verdict.
      subroutine element_inventory_report(label, rho, f_sp, carried)
      character(len=*), intent(in) :: label
      real(dp), dimension(1-Ng:N+Ng),           intent(in) :: rho
      real(dp), dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      logical,  intent(in) :: carried(n_species)
      type(element_inventory_cell) :: inv
      real(dp) :: vmax(n_inv_element), omax(n_inv_element)
      real(dp) :: vres(n_inv_element)
      real(dp) :: yhe(N), dhe
      integer  :: jv(n_inv_element), nbad(n_inv_element)
      integer  :: j, ie, jhe
      if (.not. element_inventory_report_on()) return
      vmax = 0.0d0;  omax = 0.0d0;  vres = 0.0d0
      jv   = 0;      nbad = 0
      do j = 1, N
         call element_inventory_of_cell(rho(j)*n0, f_sp(j,:), carried, inv)
         yhe(j) = helium_to_hydrogen(inv)
         do ie = 1, n_inv_element
            if (inv%violation(ie) .gt. vmax(ie)) then
               vmax(ie) = inv%violation(ie)
               jv(ie)   = j
            endif
            omax(ie) = max(omax(ie), inv%occupancy(ie))
            vres(ie) = max(vres(ie), inv%violation_less_reserved(ie))
            if (inv%violation(ie) .gt. inventory_feasible_slack)          &
               nbad(ie) = nbad(ie) + 1
         enddo
      enddo
      write(*,'(a)') ' (inventory) '//trim(label)
      do ie = 1, n_inv_element
         if (ie .eq. ien_He) cycle
         write(*,'(a,a,a,es10.3,a,i5,a,es10.3)')                          &
            '   ', trim(inv_element_name(ie)),                            &
            ': worst breach ', vmax(ie), ' at cell ', jv(ie),             &
            ', tightest occupancy ', omax(ie)
         write(*,'(a,es10.3,a,i5)')                                       &
            '       breach against the budget less the frozen stages ',   &
            vres(ie), ', cells breaching ', nbad(ie)
      enddo
      ! THE RATIO OF TWO ELEMENTS AGAINST THE RESERVOIR IT IS PINNED TO.
      ! It is a different constraint from the carrier budgets above: those
      ! bound the partition of ONE element among its species, and this one
      ! says that the amounts of two elements stand in the ratio the run
      ! was given.  A calculation that writes a species density without
      ! restoring its element total breaks this and no carrier budget
      ! notices.
      call ratio_departure_from_the_base(yhe, N, dhe, jhe)
      if (element_ratio_may_depart()) then
         write(*,'(a,es10.3,a,i5,a)')                                     &
            '   He/H departure from the base reservoir ', dhe,            &
            ' at cell ', jhe,                                             &
            ' (helium is transported, so a gradient is the solution)'
      else
         write(*,'(a,es10.3,a,i5,a)')                                     &
            '   He/H departure from the base reservoir ', dhe,            &
            ' at cell ', jhe,                                             &
            ' (no operator of this run may move it)'
      endif
      if (maxval(vmax) .le. inventory_feasible_slack .and.                &
          (element_ratio_may_depart() .or.                                &
           dhe .le. element_ratio_gate)) then
         write(*,'(a,es10.3)') '   verdict: feasible, worst breach ',      &
            maxval(vmax)
      else if (maxval(vmax) .gt. inventory_feasible_slack) then
         ie = maxloc(vmax, 1)
         write(*,'(a,es10.3,a,a,a,i0)')                                   &
            '   verdict: infeasible by ', vmax(ie), ' in ',               &
            trim(inv_element_name(ie)), ' at cell ', jv(ie)
      else
         write(*,'(a,es10.3,a,i0)')                                       &
            '   verdict: infeasible, the He/H of the state departs from'//&
            ' its reservoir by ', dhe, ' at cell ', jhe
      endif
      end subroutine element_inventory_report

      end module element_inventory

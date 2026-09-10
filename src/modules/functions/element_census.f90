      module element_census
      ! ONE definition of the elemental and charge content of a composition
      ! state (rho, f_sp), and the invariant checks built on it.
      !
      ! WHAT IT COMPUTES.  For every element the code carries -- H, He and the
      ! ten metal elements of species_table -- the NUCLEUS number density of
      ! that element in each cell, with every carrier counted at its exact
      ! stoichiometric multiplicity: H2 and H2+ two H nuclei, H3+ three,
      ! HeH+ one H and one He, OH one H and one O, H2O two H and one O, CO one
      ! C and one O, and each metal ion stage one nucleus of its element.
      ! Beside them it computes the positive charge density, which under
      ! neutrality is the free electron density.
      !
      ! The multiplicities are NOT written here: they are read from the
      ! bsp_nH / bsp_nHe / bsp_nO / bsp_nC / bsp_charge / mion_elem /
      ! mion_stage tables of species_table, so adding a carrier is adding a
      ! row there and this module follows.  Excited levels
      ! (bsp_is_excited_level, i.e. He 2^3S) are skipped: their population is
      ! already inside their parent species.
      !
      ! WHY IT EXISTS.  Every operator that moves nuclei -- the ionization
      ! equilibrium sweep, the molecular network, the carrier transport
      ! write-back, the restart equilibration -- conserves each element by
      ! construction, so any departure is a defect and not physics.  Before
      ! this module those departures could only be seen in the finished output
      ! (src/utils/element_budget.py), which says a run failed but not which
      ! refresh failed it.  element_census_take / element_census_verify put
      ! the same statement around a single call.
      !
      ! WHAT IS INVARIANT ACROSS WHAT.  Two different statements, because the
      ! operators differ:
      !
      !  * the RATIO n_El/n_H is invariant across every composition operator,
      !    including ioniz_eq -- which rewrites rho from its own calc_rho and
      !    therefore may move all densities by one common factor;
      !  * the ABSOLUTE nucleus density n_El is invariant across an operator
      !    that repartitions at fixed rho, which is what the carrier
      !    write-back and the transport step (whose write-back restores the
      !    entry element totals cell by cell) are.
      !
      ! verify() always gates the ratios; it gates the absolute densities as
      ! well when the caller states that its operator holds rho fixed.
      !
      ! TOLERANCE.  These are exact bookkeeping identities: the operators
      ! rescale stage populations onto element totals, so the only error is
      ! floating-point.  The default gate is 1.0e-9 -- three decades above the
      ! accumulated round-off measured on the marching states (<= 1e-12) and a
      ! decade below the cell solver's own xtol = sqrt(eps) = 1.5e-8, which is
      ! the largest departure a converged cell solve can leave.  It is NOT a
      ! percentage of abundance.  Override with EXHALE_ELEMENT_ASSERT_TOL.
      !
      ! COST AND CONTROL.  Off unless EXHALE_ELEMENT_ASSERT is 1 (report and
      ! continue) or 2 (report and stop).  When off, take() and verify() do
      ! nothing and allocate nothing.

      use global_parameters
      use species_table, only: n_bsp, bsp_fsp, bsp_nH, bsp_nHe, bsp_nO,   &
                               bsp_nC, bsp_charge, bsp_mass, bsp_name,    &
                               bsp_is_excited_level,                      &
                               n_mion, mion_fsp, mion_elem, mion_stage,   &
                               mion_name, n_melem, melem_A, melem_name,  &
                               iel_O, iel_C

      implicit none
      private

      ! H and He first, then the metal elements in species_table order.
      integer, parameter, public :: n_element = 2 + n_melem
      integer, parameter, public :: ie_H  = 1
      integer, parameter, public :: ie_He = 2
      ! Metal element im of species_table is element ie_metal0 + im here.
      integer, parameter, public :: ie_metal0 = 2

      type, public :: element_census_state
         ! Nucleus density of each element [cm^-3] and the positive charge
         ! density [cm^-3] of the state this was taken from.
         real*8, allocatable :: n_nuc(:,:)
         real*8, allocatable :: n_chg(:)
         ! Mass density the same state implies through calc_rho's weights
         ! [m_H cm^-3], and the density it was taken at [same units].
         real*8, allocatable :: rho_comp(:)
         real*8, allocatable :: rho_in(:)
         logical             :: taken = .false.
         character(len=64)   :: label = ''
      end type element_census_state

      ! THE GATE ON AN ELEMENT RATIO, one definition for the whole code.
      ! An element ratio n_El/n_H is an exact bookkeeping identity carried
      ! across every operator of a run, so what separates a conserving state
      ! from a defect is the accumulated round-off of thousands of
      ! repartitions: 1e-9 relative, three decades above the round-off
      ! measured on marching states and a decade below the cell solver's own
      ! xtol = sqrt(eps). element_census_tolerance is the same number with
      ! the measurement override applied; the parameter is what the callers
      ! that need a compile-time constant read (the species floor of the
      ! temporal error estimate, the element inventory's ratio gate).
      real*8, parameter, public :: element_ratio_gate = 1.0d-9

      public :: element_census_on, element_census_fatal
      public :: element_nuclei_and_charge, element_census_take
      public :: element_census_verify, element_census_reservoir
      public :: element_name, element_census_tolerance

      contains

      ! ------------------------------------------------------!

      logical function element_census_on()
      ! EXHALE_ELEMENT_ASSERT = 1 (report) or 2 (report and stop).
      character(len=8) :: env
      call get_environment_variable('EXHALE_ELEMENT_ASSERT', env)
      element_census_on = (trim(env) .eq. '1') .or. (trim(env) .eq. '2')
      end function element_census_on

      ! ------------------------------------------------------!

      logical function element_census_fatal()
      character(len=8) :: env
      call get_environment_variable('EXHALE_ELEMENT_ASSERT', env)
      element_census_fatal = (trim(env) .eq. '2')
      end function element_census_fatal

      ! ------------------------------------------------------!

      real*8 function element_census_tolerance()
      character(len=32) :: env
      real*8 :: v
      integer :: ios
      element_census_tolerance = element_ratio_gate
      call get_environment_variable('EXHALE_ELEMENT_ASSERT_TOL', env)
      if (len_trim(env) .gt. 0) then
         read(env,*,iostat=ios) v
         if (ios .eq. 0 .and. v .gt. 0.0d0) element_census_tolerance = v
      endif
      end function element_census_tolerance

      ! ------------------------------------------------------!

      character(len=5) function element_name(ie)
      integer, intent(in) :: ie
      if (ie .eq. ie_H) then
         element_name = 'H    '
      else if (ie .eq. ie_He) then
         element_name = 'He   '
      else
         element_name = melem_name(ie - ie_metal0)//'   '
      endif
      end function element_name

      ! ------------------------------------------------------!

      subroutine element_nuclei_and_charge(rho, f_sp, n_nuc, n_chg,       &
                                           rho_comp)
      ! Nucleus density of every element, the positive charge density, and
      ! the mass density the composition implies, all in cm^-3 / m_H cm^-3.
      !
      ! rho is the code's adimensional mass density, so rho*n0 is the number
      ! of m_H per cm^3 that f_sp is a fraction OF -- the same nd = rho*n0 the
      ! carrier write-back and ioniz_eq use.
      !
      ! rho_comp is calc_rho's own weighted sum rebuilt from the same table
      ! (bsp_mass for the H/He/molecular/carrier species, melem_A per metal
      ! nucleus under the eos_metals policy).  Comparing it with rho*n0 is the
      ! mass closure that catches a moved n_H denominator directly, which is
      ! the signature a pure ratio test reports as every element moving by the
      ! same amount.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng,n_element), intent(out) :: n_nuc
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: n_chg
      real*8, dimension(1-Ng:N+Ng), optional, intent(out) :: rho_comp

      real*8, dimension(1-Ng:N+Ng) :: nd, c, rc
      integer :: ib, im

      nd    = rho*n0
      n_nuc = 0.0d0
      n_chg = 0.0d0
      rc    = 0.0d0

      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         c = f_sp(:,bsp_fsp(ib))*nd
         if (bsp_nH(ib)  .gt. 0)                                          &
            n_nuc(:,ie_H)  = n_nuc(:,ie_H)  + dble(bsp_nH(ib)) *c
         if (bsp_nHe(ib) .gt. 0)                                          &
            n_nuc(:,ie_He) = n_nuc(:,ie_He) + dble(bsp_nHe(ib))*c
         if (bsp_nO(ib)  .gt. 0)                                          &
            n_nuc(:,ie_metal0+iel_O) = n_nuc(:,ie_metal0+iel_O)           &
                                     + dble(bsp_nO(ib))*c
         if (bsp_nC(ib)  .gt. 0)                                          &
            n_nuc(:,ie_metal0+iel_C) = n_nuc(:,ie_metal0+iel_C)           &
                                     + dble(bsp_nC(ib))*c
         n_chg = n_chg + dble(bsp_charge(ib))*c
         rc    = rc    + bsp_mass(ib)*c
      enddo

      ! Metal ion stages: one nucleus of their element each, charge = stage.
      ! The mass they carry enters rho only under the eos_metals policy, and
      ! the C and O nuclei bound in the carriers above already carried theirs
      ! through bsp_mass (which is why bsp_mass(CO) = 1 + melem_A(C) +
      ! melem_A(O) - 1: see the species_table comment).
      do im = 1, n_mion
         c = f_sp(:,mion_fsp(im))*nd
         n_nuc(:,ie_metal0+mion_elem(im)) =                               &
            n_nuc(:,ie_metal0+mion_elem(im)) + c
         n_chg = n_chg + dble(mion_stage(im))*c
         if (eos_include_metals .and. thereis_metals)                     &
            rc = rc + melem_A(mion_elem(im))*c
      enddo

      if (present(rho_comp)) rho_comp = rc

      end subroutine element_nuclei_and_charge

      ! ------------------------------------------------------!

      subroutine element_census_take(label, rho, f_sp, snap)
      ! Record the elemental content of (rho, f_sp) under a label.  A no-op
      ! unless the check is switched on.
      character(len=*),                       intent(in)    :: label
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: rho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)    :: f_sp
      type(element_census_state),             intent(inout) :: snap

      snap%taken = .false.
      if (.not. element_census_on()) return
      if (.not. allocated(snap%n_nuc)) then
         allocate(snap%n_nuc(1-Ng:N+Ng,n_element))
         allocate(snap%n_chg(1-Ng:N+Ng))
         allocate(snap%rho_comp(1-Ng:N+Ng))
         allocate(snap%rho_in(1-Ng:N+Ng))
      endif
      call element_nuclei_and_charge(rho, f_sp, snap%n_nuc, snap%n_chg,   &
                                     snap%rho_comp)
      snap%rho_in = rho*n0
      snap%label  = label
      snap%taken  = .true.

      end subroutine element_census_take

      ! ------------------------------------------------------!

      subroutine element_census_verify(snap, rho, f_sp, rho_is_fixed)
      ! Compare the state (rho, f_sp) against the recorded snapshot.
      !
      ! Gated always: the ratio n_El/n_H of every element other than hydrogen.
      ! Gated when rho_is_fixed: the absolute nucleus density of every
      ! element, hydrogen included.
      ! Reported always: the mass closure sum(species mass)/rho of the state,
      ! which is where a moved n_H denominator shows up.
      !
      ! On failure the element, the cell, its radius, the before/after totals
      ! and every species that carries the element in that cell are printed,
      ! so the offending refresh names itself.
      type(element_census_state),             intent(in) :: snap
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      logical, optional,                      intent(in) :: rho_is_fixed

      real*8, dimension(1-Ng:N+Ng,n_element) :: n_nuc
      real*8, dimension(1-Ng:N+Ng)           :: n_chg, rho_comp
      real*8  :: tol, dev, dev_worst, r_bef, r_aft, mclos, mclos_worst
      integer :: ie, j, jw, ie_worst, n_bad
      logical :: absolute

      if (.not. element_census_on()) return
      if (.not. snap%taken)          return

      absolute = .false.
      if (present(rho_is_fixed)) absolute = rho_is_fixed
      tol = element_census_tolerance()

      call element_nuclei_and_charge(rho, f_sp, n_nuc, n_chg, rho_comp)

      n_bad       = 0
      dev_worst   = 0.0d0
      ie_worst    = 0
      jw          = 0
      mclos_worst = 0.0d0

      do j = 1-Ng, N+Ng
         if (snap%n_nuc(j,ie_H) .le. 0.0d0) cycle
         mclos = abs(rho_comp(j) - rho(j)*n0)/max(rho(j)*n0, 1.0d-99)
         if (mclos .gt. mclos_worst) mclos_worst = mclos
         do ie = 1, n_element
            if (snap%n_nuc(j,ie) .le. 0.0d0) cycle
            dev = 0.0d0
            if (ie .ne. ie_H) then
               r_bef = snap%n_nuc(j,ie)/snap%n_nuc(j,ie_H)
               r_aft = n_nuc(j,ie)/max(n_nuc(j,ie_H), 1.0d-99)
               dev   = abs(r_aft - r_bef)/r_bef
            endif
            if (absolute) dev = max(dev, abs(n_nuc(j,ie)                  &
                                 - snap%n_nuc(j,ie))/snap%n_nuc(j,ie))
            if (dev .gt. dev_worst) then
               dev_worst = dev;  ie_worst = ie;  jw = j
            endif
            if (dev .gt. tol) n_bad = n_bad + 1
         enddo
      enddo

      if (n_bad .eq. 0) then
         if (element_census_fatal())                                      &
            write(*,'(A,A,A,ES9.2,A,ES9.2)')                              &
               ' (element_census) ', trim(snap%label),                    &
               ': closed, worst ', dev_worst, ', mass closure ', mclos_worst
         return
      endif

      write(*,'(A)') ' ================================================='//&
                     '================='
      write(*,'(A,A)') ' (element_census) ELEMENT BUDGET BROKEN across: ', &
           trim(snap%label)
      write(*,'(A,I0,A,ES10.3,A,ES10.3)') '   ', n_bad,                    &
           ' (element,cell) pair(s) outside tol =', tol,                   &
           ';  worst departure ', dev_worst
      write(*,'(A,ES10.3)') '   mass closure |rho(composition)/rho - 1|'// &
           ' worst on the grid: ', mclos_worst
      if (ie_worst .gt. 0) then
         write(*,'(A,A,A,I5,A,F10.5)') '   worst element ',                &
              trim(element_name(ie_worst)), ' at cell ', jw,               &
              '  r = ', r(jw)
         write(*,'(A,ES16.8,A,ES16.8)') '     n_El before ',               &
              snap%n_nuc(jw,ie_worst), '   after ', n_nuc(jw,ie_worst)
         write(*,'(A,ES16.8,A,ES16.8)') '     n_H  before ',               &
              snap%n_nuc(jw,ie_H), '   after ', n_nuc(jw,ie_H)
         write(*,'(A,ES16.8,A,ES16.8)') '     rho*n0    ',                 &
              snap%rho_in(jw), ' -> ', rho(jw)*n0
         write(*,'(A,ES16.8,A,ES16.8)') '     rho(comp) ',                 &
              snap%rho_comp(jw), ' -> ', rho_comp(jw)
         call report_carriers(jw, ie_worst, snap, f_sp, rho)
      endif
      ! Every element that missed, so a common shift (all of them, equally)
      ! is told apart from a real repartition defect (one or two of them).
      write(*,'(A)') '   worst departure of each element over the grid:'
      do ie = 1, n_element
         dev_worst = 0.0d0
         jw = 0
         do j = 1-Ng, N+Ng
            if (snap%n_nuc(j,ie) .le. 0.0d0) cycle
            if (snap%n_nuc(j,ie_H) .le. 0.0d0) cycle
            dev = 0.0d0
            if (ie .ne. ie_H) then
               r_bef = snap%n_nuc(j,ie)/snap%n_nuc(j,ie_H)
               r_aft = n_nuc(j,ie)/max(n_nuc(j,ie_H), 1.0d-99)
               dev   = abs(r_aft - r_bef)/r_bef
            endif
            if (absolute) dev = max(dev, abs(n_nuc(j,ie)                  &
                                 - snap%n_nuc(j,ie))/snap%n_nuc(j,ie))
            if (dev .gt. dev_worst) then
               dev_worst = dev;  jw = j
            endif
         enddo
         if (jw .ne. 0) write(*,'(A,A,A,ES10.3,A,I5,A,F9.5)')             &
              '     ', trim(element_name(ie)), '  ', dev_worst,           &
              '  at cell ', jw, ' r =', r(jw)
      enddo
      write(*,'(A)') ' ================================================='//&
                     '================='

      if (element_census_fatal()) then
         write(*,'(A)') ' (element_census) EXHALE_ELEMENT_ASSERT=2: stop.'
         stop 1
      endif

      end subroutine element_census_verify

      ! ------------------------------------------------------!

      subroutine element_census_reservoir(label, rho, f_sp)
      ! The ABSOLUTE statement, against the reservoirs the run resolved --
      ! the same test src/utils/element_budget.py applies to the finished
      ! output, made here on the live state so that the refresh that broke it
      ! can still be named:
      !
      !     n_El/n_H  =  (El/H)_resolved      (He/H, and melem_ab per metal)
      !     n_H       =  rho / mass_per_H
      !
      ! Helium is exempt when the element-diffusion operator is on: He/H is
      ! then a solved profile, not a column invariant (element_budget.py makes
      ! the same exemption and checks the base cell instead).
      character(len=*),                       intent(in) :: label
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp

      real*8, dimension(1-Ng:N+Ng,n_element) :: n_nuc
      real*8, dimension(1-Ng:N+Ng)           :: n_chg
      real*8  :: tol, res, dev, dev_worst
      integer :: ie, j, jw, n_bad

      if (.not. element_census_on()) return
      tol = element_census_tolerance()
      call element_nuclei_and_charge(rho, f_sp, n_nuc, n_chg)

      n_bad = 0
      do ie = 1, n_element
         if (ie .eq. ie_H) cycle
         if (ie .eq. ie_He) then
            if (.not. thereis_He) cycle
            if (he_diffusion)     cycle
            res = HeH
         else
            if (.not. thereis_metals) cycle
            res = melem_ab(ie - ie_metal0)
         endif
         if (res .le. 0.0d0) cycle
         dev_worst = 0.0d0
         jw = 0
         do j = 1, N
            if (n_nuc(j,ie_H) .le. 0.0d0) cycle
            dev = abs(n_nuc(j,ie)/n_nuc(j,ie_H) - res)/res
            if (dev .gt. dev_worst) then
               dev_worst = dev;  jw = j
            endif
         enddo
         if (dev_worst .gt. tol) then
            n_bad = n_bad + 1
            write(*,'(A,A,A,A,ES10.3,A,ES10.3,A,I5,A,F9.5)')             &
                 ' (element_census) ', trim(label), ': ',                &
                 trim(element_name(ie)), ' misses its reservoir ', res,  &
                 ' by ', dev_worst, ' at cell ', jw, ' r =', r(jw)
         endif
      enddo

      ! Hydrogen against the mass density, the closure element_budget.py
      ! states as rho/mass_per_H = n_H.
      dev_worst = 0.0d0
      jw = 0
      do j = 1, N
         if (n_nuc(j,ie_H) .le. 0.0d0) cycle
         dev = abs(rho(j)*n0/mass_per_H - n_nuc(j,ie_H))/n_nuc(j,ie_H)
         if (dev .gt. dev_worst) then
            dev_worst = dev;  jw = j
         endif
      enddo
      if (dev_worst .gt. tol .and. .not. he_diffusion) then
         n_bad = n_bad + 1
         write(*,'(A,A,A,ES10.3,A,I5,A,F9.5)')                           &
              ' (element_census) ', trim(label),                         &
              ': H misses rho/mass_per_H by ', dev_worst,                &
              ' at cell ', jw, ' r =', r(jw)
      endif

      if (n_bad .gt. 0 .and. element_census_fatal()) then
         write(*,'(A)') ' (element_census) EXHALE_ELEMENT_ASSERT=2: stop.'
         stop 1
      endif

      end subroutine element_census_reservoir

      ! ------------------------------------------------------!

      subroutine report_carriers(j, ie, snap, f_sp, rho)
      ! Every species that carries element ie in cell j, before and after.
      integer,                                intent(in) :: j, ie
      type(element_census_state),             intent(in) :: snap
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in) :: rho
      integer :: ib, im, mult
      real*8  :: nd

      nd = rho(j)*n0
      write(*,'(A)') '     contributing species [cm^-3, after]:'
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         mult = 0
         if (ie .eq. ie_H)  mult = bsp_nH(ib)
         if (ie .eq. ie_He) mult = bsp_nHe(ib)
         if (ie .eq. ie_metal0+iel_O) mult = bsp_nO(ib)
         if (ie .eq. ie_metal0+iel_C) mult = bsp_nC(ib)
         if (mult .le. 0) cycle
         write(*,'(A,A,A,I2,A,ES16.8)') '       ', bsp_name(ib),          &
              '  x', mult, '  ', dble(mult)*f_sp(j,bsp_fsp(ib))*nd
      enddo
      do im = 1, n_mion
         if (ie_metal0 + mion_elem(im) .ne. ie) cycle
         write(*,'(A,A,A,ES16.8)') '       ', mion_name(im),              &
              '       ', f_sp(j,mion_fsp(im))*nd
      enddo
      ! Silence an unused-argument warning while keeping the snapshot in the
      ! interface: the before-state totals are printed by the caller.
      if (.false.) write(*,*) snap%taken

      end subroutine report_carriers

      end module element_census

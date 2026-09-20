      program carrier_helium_inventory
      ! THE HELIUM A CELL HOLDS IS THE HELIUM IT HAD, whatever the carrier
      ! solve does to the partition of that helium among its stages.
      !
      ! The two ionized stages are carried as fractions of the HELIUM
      ! NUCLEUS COUNT of the cell, and that count includes the nucleus
      ! inside HeH+.  The nucleus in a molecule is not theirs to take: the
      ! write-back leaves the molecule where it is, because rescaling it
      ! would move a hydrogen nucleus with it.  Neither is the neutral
      ! helium held in the He 2^3S metastable, which this step freezes and
      ! whose population the chemistry rows subtract from the neutral stage
      ! to form the ground singlet they are written in.  So the admissible
      ! partition is
      !
      !    x(He II) + x(He III) <= 1 - (n_HeH+ + n_He(2^3S))/n_He,nuc,
      !
      ! and a projection that stops at one admits states that hold more
      ! helium than the cell has, or whose ground singlet is negative.
      !
      ! WHAT IS MEASURED HERE IS THE STATE THE PRODUCTION ROUTINES RETURN,
      ! species by species: the helium nuclei, the hydrogen nuclei, the
      ! mass and the charge formed from f_sp after limit_to_element_budget
      ! and carrier_write_back have run on it.  An algebraically closed sum
      ! of a remaining stage and two fluxes would state none of these.
      !
      ! The four cases:
      !
      !  (1) THE INVENTORY.  A cell whose HeH+ holds a tenth of its helium
      !      and whose two ionized stages are handed 0.60 and 0.35 of it.
      !      The sum 0.95 is inside the interval [0,1] and outside the
      !      admissible set, and the state written from it holds 1.05 of
      !      the helium the cell has.
      !
      !  (2) THE LEVEL.  A cell whose neutral helium is nearly exhausted
      !      and almost all of it metastable.  A partition the interval
      !      admits leaves less neutral helium than the frozen level holds,
      !      so the ground singlet the rows are evaluated at is negative;
      !      the old closure returned zero there and the trial stood.
      !
      !  (3) THE CENSUS of the column: hydrogen nuclei, helium nuclei and
      !      mass unchanged by the write-back in every cell it owns, the
      !      lower ghosts -- which hold the inflow composition and which
      !      this operator does not own -- untouched, and the electron
      !      density of the returned species the one the species imply.
      !
      !  (4) THE SOURCE LEDGER.  sprod - sloss is the source the residual
      !      is assembled from, for every carrier, in molecular and in
      !      atomic gas.  The proton row of a molecular cell is split by
      !      mol_heh_rows and then MOVED by the He <-> H charge-exchange
      !      pair, so a record taken from the split alone is not the row.
      !
      ! The gas is the metal-free hydrogen and helium column of
      ! carrier_reference_scales, with the molecular chemistry, the He 2^3S
      ! level and the transported ionization state on.  The frozen cell
      ! state is defined rather than realistic; nothing asserted here
      ! depends on the size of a rate except that the charge-exchange pair
      ! is not zero, which case (4) states on its own.

      use global_parameters
      use species_table
      use utils, only: calc_ne
      use charge_exchange, only: he_h_charge_exchange
      use ionization_equilibrium, only: bg_cell, ioniz_eq_allocate_arrays
      use diffusive_photochemistry, only: carrier_set_init, carrier_state, &
                                          carrier_source,                  &
                                          carrier_write_back,              &
                                          limit_to_element_budget,         &
                                          carrier_nucleus_reference,       &
                                          carrier_helium_singlet_breach,   &
                                          n_carrier, n_carrier_max,        &
                                          carrier_solved, carrier_name,    &
                                          ic_H2, ic_OH, ic_H2O, ic_CO,     &
                                          ic_Hp, ic_HeII, ic_HeIII
      use assertion_report
      implicit none

      ! The cell holding the inventory counterexample and the cell holding
      ! the exhausted neutral helium.  They are not neighbours: every
      ! routine under test is cell-local, so a defect that leaked would
      ! leak to a neighbour first.
      integer, parameter :: j_mol   = 2
      integer, parameter :: j_tight = 4
      ! The counterexample, as fractions of the helium nuclei of the cell
      ! (the review's arithmetic: 0.10 of the helium in HeH+, and the two
      ! ionized stages handed 0.60 and 0.35 of the element).
      real*8,  parameter :: x_heii_over  = 0.60d0
      real*8,  parameter :: x_heiii_over = 0.35d0
      ! The stage fractions handed to the exhausted cell: inside [0,1] and
      ! outside the admissible set of a cell whose neutral helium is almost
      ! all metastable.
      real*8,  parameter :: x_heii_tight  = 0.98d0
      real*8,  parameter :: x_heiii_tight = 0.01d0
      real*8,  parameter :: T_gas = 6.0d3

      real*8, allocatable :: rho(:), f_sp(:,:), f_entry(:,:), fc(:,:)
      real*8, allocatable :: ntot(:), TK(:), mbar(:), nrho(:), wfac(:)
      real*8, allocatable :: nH_free(:), nO_free(:), nC_free(:)
      real*8, allocatable :: nucH0(:), nucHe0(:), mass0(:), chg0(:)
      real*8 :: nc(n_carrier_max), src(n_carrier_max)
      real*8 :: sprod(n_carrier_max), sloss(n_carrier_max)
      real*8 :: xsum, xbound, singlet, dH, dHe, dM, ghost_move
      integer :: j, ic, jw
      ! The cell-state arrays of the ionization module are allocated once;
      ! the second configuration reuses them.
      logical :: ion_arrays_ready = .false.

      call setup_globals(.true.)
      call build_helium_column()

      ! The pair whose omission from the split the ledger case is about is
      ! on by default; a run with it off would make case (4) vacuous.
      call check_absolute('he_h_charge_exchange_is_on',                  &
                          logical_as_double(he_h_charge_exchange),       &
                          1.0d0, 0.0d0)

      ! ---- (4) the source ledger of a molecular cell ------------------ !
      ! Taken before anything moves the state, at the entry partition of a
      ! cell whose helium and hydrogen are both well inside their budgets.
      call ledger_identity_of_cell(1, 'molecular')

      ! ---- the two cells are set up as the cases describe -------------- !
      fc(j_mol,ic_HeII)    = x_heii_over
      fc(j_mol,ic_HeIII)   = x_heiii_over
      fc(j_tight,ic_HeII)  = x_heii_tight
      fc(j_tight,ic_HeIII) = x_heiii_tight
      ! The bound is formed HERE, from the composition this suite built,
      ! and not taken from the module under test.
      xsum   = fc(j_mol,ic_HeII) + fc(j_mol,ic_HeIII)
      xbound = 1.0d0 - (f_entry(j_mol,isp_HeHp)                          &
                      + f_entry(j_mol,isp_HeTR))/nucHe0(j_mol)
      call check_positive('the_counterexample_is_inside_the_interval',   &
                          1.0d0 - xsum)
      call check_positive('the_counterexample_is_outside_the_set',       &
                          xsum - xbound)

      ! ---- the production path: the limiter, then the write-back ------ !
      call limit_to_element_budget(fc, nrho, TK, nH_free, nO_free,       &
                                   nC_free)

      ! (1) the projection returns the counterexample to the admissible set
      xsum = fc(j_mol,ic_HeII) + fc(j_mol,ic_HeIII)
      call check_at_most('projected_stage_sum_within_the_helium_available',&
                         xsum - xbound, 1.0d-12)

      ! (2) the ground singlet the exhausted cell's partition leaves
      singlet = helium_singlet_of_returned_partition(j_tight)
      call check_at_least('ground_singlet_of_the_exhausted_cell',        &
                          singlet, -1.0d-12)
      singlet = helium_singlet_of_returned_partition(j_mol)
      call check_at_least('ground_singlet_of_the_counterexample_cell',   &
                          singlet, -1.0d-12)

      call carrier_write_back(rho, f_sp, fc, nrho, nH_free, nO_free,     &
                              nC_free)

      ! (1) the helium the write-back handed back, species by species
      call worst_census_departure(dH, dHe, dM, jw)
      call check_at_most('helium_nuclei_of_the_written_state', dHe,      &
                         1.0d-12)
      call check_at_most('hydrogen_nuclei_of_the_written_state', dH,     &
                         1.0d-12)
      call check_at_most('mass_of_the_written_state', dM, 1.0d-12)

      ! (2) the level stays inside the stage that contains it
      call check_at_least('metastable_inside_the_neutral_helium',        &
                          f_sp(j_tight,isp_HeI) - f_sp(j_tight,isp_HeTR),&
                          0.0d0)

      ! the operator's own measure of the breach, on the state it wrote
      call check_at_most('write_back_helium_breach_measure',             &
                         carrier_helium_singlet_breach(), 1.0d-12)

      ! (3) the ghosts, by their ownership: the lower ones hold the inflow
      ! composition and are not this operator's to write
      ghost_move = 0.0d0
      do j = 1-Ng, 0
         ghost_move = max(ghost_move,                                     &
                          maxval(abs(f_sp(j,:) - f_entry(j,:))))
      enddo
      call check_absolute('lower_ghosts_are_not_written', ghost_move,    &
                          0.0d0, 0.0d0)

      ! (3) the electron budget of the returned species
      call check_at_most('electron_density_of_the_returned_species',     &
                         electron_density_departure(), 1.0d-14)

      ! ---- (4) the same ledger in an atomic gas ----------------------- !
      deallocate(rho, f_sp, f_entry, fc, ntot, TK, mbar, nrho, wfac,     &
                 nH_free, nO_free, nC_free, nucH0, nucHe0, mass0, chg0)
      call setup_globals(.false.)
      call build_helium_column()
      call ledger_identity_of_cell(1, 'atomic')

      if (assertion_failures .gt. 0) stop 1

      contains

      !--------------!

      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double

      !--------------!

      subroutine ledger_identity_of_cell(jcell, gas)
      ! sprod - sloss = src for every solved carrier of one cell, measured
      ! against the size of the two entries so that a row whose halves
      ! cancel is judged on its own scale.
      integer,          intent(in) :: jcell
      character(len=*), intent(in) :: gas
      real*8 :: dmax, d, scale
      integer :: icw
      call cell_densities(jcell, nc)
      call carrier_source(jcell, nc, nH_free(jcell), nO_free(jcell),     &
                          src, sprod = sprod, sloss = sloss)
      dmax = 0.0d0
      icw  = 1
      do ic = 1, n_carrier
         if (.not. carrier_solved(ic)) cycle
         scale = max(abs(sprod(ic)), abs(sloss(ic)), abs(src(ic)),       &
                     1.0d-300)
         d = abs(sprod(ic) - sloss(ic) - src(ic))/scale
         if (d .gt. dmax) then
            dmax = d
            icw  = ic
         endif
      enddo
      call check_at_most('source_ledger_'//trim(gas)//'_worst_'//        &
                         trim(carrier_name(icw)), dmax, 1.0d-12)
      end subroutine ledger_identity_of_cell

      !--------------!

      double precision function helium_singlet_of_returned_partition(j)  &
                               result(x)
      ! The ground singlet the returned stage fractions leave, as a
      ! fraction of the helium nuclei of the cell: one less the two ionized
      ! stages, the helium inside HeH+ and the frozen metastable, all four
      ! read from the state the routines were handed.  Negative says the
      ! returned partition is one no cell can be in.
      integer, intent(in) :: j
      x = 1.0d0 - fc(j,ic_HeII) - fc(j,ic_HeIII)                         &
        - (f_entry(j,isp_HeHp) + f_entry(j,isp_HeTR))/nucHe0(j)
      end function helium_singlet_of_returned_partition

      !--------------!

      subroutine worst_census_departure(dH, dHe, dM, jworst)
      ! The largest relative change of the hydrogen nuclei, the helium
      ! nuclei and the mass of a cell this operator owns, between the state
      ! it was handed and the state it wrote.
      real*8,  intent(out) :: dH, dHe, dM
      integer, intent(out) :: jworst
      real*8 :: nh, nhe, m, c, fs(n_species)
      integer :: j
      dH = 0.0d0;  dHe = 0.0d0;  dM = 0.0d0;  jworst = 0
      do j = 1, N+Ng
         fs = f_sp(j,:)
         call census_of_cell(fs, nh, nhe, m, c)
         if (abs(nhe/nucHe0(j) - 1.0d0) .gt. dHe) jworst = j
         dH  = max(dH,  abs(nh /nucH0(j)  - 1.0d0))
         dHe = max(dHe, abs(nhe/nucHe0(j) - 1.0d0))
         dM  = max(dM,  abs(m  /mass0(j)  - 1.0d0))
      enddo
      end subroutine worst_census_departure

      !--------------!

      double precision function electron_density_departure() result(d)
      ! The free electrons the returned species carry, taken cell by cell
      ! from the charges of the species table, against the column the
      ! production routine calc_ne builds from the same species.  Charge is
      ! not conserved by the write-back -- the operator transports the
      ! ionization partition, which is what a charge is -- so what is
      ! asserted is that the electron budget the next sweep reads is the
      ! one the written species imply.
      real*8, dimension(1-Ng:N+Ng) :: ne
      real*8, dimension(1-Ng:N+Ng,4) :: nmol
      real*8 :: nh, nhe, m, c, fs(n_species)
      integer :: j
      nmol(:,1) = f_sp(:,isp_H2)
      nmol(:,2) = f_sp(:,isp_H2p)
      nmol(:,3) = f_sp(:,isp_H3p)
      nmol(:,4) = f_sp(:,isp_HeHp)
      call calc_ne(f_sp(:,isp_HII), f_sp(:,isp_HeII), f_sp(:,isp_HeIII), &
                   ne, nmol = nmol)
      d = 0.0d0
      do j = 1, N+Ng
         fs = f_sp(j,:)
         call census_of_cell(fs, nh, nhe, m, c)
         d = max(d, abs(c - ne(j))/max(abs(ne(j)), 1.0d-300))
         if (.not. (c .ge. 0.0d0)) d = 1.0d0
      enddo
      end function electron_density_departure

      !--------------!

      subroutine census_of_cell(fs, nucH, nucHe, mass, charge)
      ! The hydrogen nuclei, the helium nuclei, the mass and the free
      ! electrons of one cell's composition, per unit mass, from the
      ! stoichiometry of the species table.  The He 2^3S column is a level
      ! inside He I and its nucleus, mass and charge are already counted
      ! there, so every sum skips it.
      real*8, intent(in)  :: fs(1:n_species)
      real*8, intent(out) :: nucH, nucHe, mass, charge
      integer :: ib
      nucH = 0.0d0;  nucHe = 0.0d0;  mass = 0.0d0;  charge = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         nucH   = nucH   + dble(bsp_nH(ib))    *fs(bsp_fsp(ib))
         nucHe  = nucHe  + dble(bsp_nHe(ib))   *fs(bsp_fsp(ib))
         mass   = mass   +      bsp_mass(ib)   *fs(bsp_fsp(ib))
         charge = charge + dble(bsp_charge(ib))*fs(bsp_fsp(ib))
      enddo
      end subroutine census_of_cell

      !--------------!

      subroutine cell_densities(j, nc)
      ! The species density of every carrier of cell j, each against its
      ! own reference density: the element nucleus density for an
      ! ionization stage and the density of the mass row for the rest.
      integer, intent(in)  :: j
      real*8,  intent(out) :: nc(n_carrier_max)
      integer :: ic
      do ic = 1, n_carrier_max
         nc(ic) = fc(j,ic)*carrier_nucleus_reference(ic, j, nrho(j))
      enddo
      end subroutine cell_densities

      !--------------!

      subroutine setup_globals(molecular)
      ! A five-cell hydrogen and helium column with the He 2^3S level and
      ! the transported ionization state on, so that both helium stages are
      ! carried, and with the molecular network selected by the argument.
      logical, intent(in) :: molecular
      integer :: j
      N   = 5
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.1d0
      p_base_bar         = 1.0d-6
      thereis_He         = .true.
      thereis_HeITR      = .true.
      thereis_metals     = .false.
      thereis_mol        = molecular
      thereis_oxychem    = .false.
      eos_include_metals = .false.
      he_diffusion       = .false.
      carrier_transport    = .true.
      ionization_transport = .true.
      if (.not. allocated(r)) then
         allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
         do j = 1-Ng, N+Ng
            r(j)     = 1.0d0 + 0.01d0*dble(j)
            r_edg(j) = 1.0d0 + 0.01d0*(dble(j) + 0.5d0)
            dr_j(j)  = 0.01d0
         enddo
      endif
      if (.not. allocated(melem_ab)) then
         allocate(melem_ab(n_melem))
         melem_ab = 0.0d0
      endif
      call carrier_set_init()
      if (.not. ion_arrays_ready) then
         call ioniz_eq_allocate_arrays()
         ion_arrays_ready = .true.
      endif
      end subroutine setup_globals

      !--------------!

      subroutine build_helium_column()
      ! A partly ionized molecular base gas.  A tenth of the helium of
      ! every cell is bound in HeH+ and a fiftieth of it sits in the
      ! metastable level inside He I; cell j_tight instead carries nearly
      ! all its helium ionized, with what neutral helium is left almost
      ! entirely metastable.  The outer ghost carries a composition of its
      ! own, so that a census taken there is not a copy of cell N's.
      real*8  :: fs(n_species)
      integer :: j
      allocate(rho(1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(f_entry(1-Ng:N+Ng,n_species))
      allocate(fc(1-Ng:N+Ng,n_carrier_max))
      allocate(ntot(1-Ng:N+Ng), TK(1-Ng:N+Ng), mbar(1-Ng:N+Ng))
      allocate(nrho(1-Ng:N+Ng), wfac(1-Ng:N+Ng))
      allocate(nH_free(1-Ng:N+Ng), nO_free(1-Ng:N+Ng), nC_free(1-Ng:N+Ng))
      allocate(nucH0(1-Ng:N+Ng), nucHe0(1-Ng:N+Ng))
      allocate(mass0(1-Ng:N+Ng), chg0(1-Ng:N+Ng))
      rho  = 1.0d0
      f_sp = 0.0d0
      f_sp(:,isp_HI)    = 0.20d0
      f_sp(:,isp_HII)   = 0.02d0
      f_sp(:,isp_HeI)   = 0.065d0
      f_sp(:,isp_HeTR)  = 0.002d0
      f_sp(:,isp_HeII)  = 0.020d0
      f_sp(:,isp_HeIII) = 0.005d0
      f_sp(:,isp_HeHp)  = 0.010d0
      if (thereis_mol) then
         f_sp(:,isp_H2)  = 0.30d0
         f_sp(:,isp_H2p) = 1.0d-6
         f_sp(:,isp_H3p) = 1.0d-6
      else
         f_sp(:,isp_HI)  = 0.80d0
         f_sp(:,isp_HeHp) = 0.0d0
         f_sp(:,isp_HeI) = 0.075d0
      endif
      ! The cell whose neutral helium is nearly exhausted and nearly all
      ! metastable.  Its helium nucleus total is that of every other cell.
      f_sp(j_tight,isp_HeHp)  = 0.0d0
      f_sp(j_tight,isp_HeI)   = 0.005d0
      f_sp(j_tight,isp_HeTR)  = 0.0045d0
      f_sp(j_tight,isp_HeII)  = 0.090d0
      f_sp(j_tight,isp_HeIII) = 0.005d0
      ! The outer ghost: the same elements in a different partition.
      do j = N+1, N+Ng
         f_sp(j,isp_HI)    = f_sp(j,isp_HI)    + 0.02d0
         f_sp(j,isp_HII)   = f_sp(j,isp_HII)   - 0.02d0
         f_sp(j,isp_HeI)   = f_sp(j,isp_HeI)   + 0.004d0
         f_sp(j,isp_HeIII) = f_sp(j,isp_HeIII) - 0.004d0
      enddo
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
      enddo
      call carrier_state(rho, f_sp, fc, ntot, nrho, wfac, TK, mbar,      &
                         nH_free, nO_free, nC_free)
      TK = T_gas
      f_entry = f_sp
      do j = 1-Ng, N+Ng
         fs = f_sp(j,:)
         call census_of_cell(fs, nucH0(j), nucHe0(j), mass0(j), chg0(j))
      enddo
      end subroutine build_helium_column

      !--------------!

      subroutine set_frozen_cell_rates(j)
      ! Every field of the frozen cell state, so that none is read
      ! undefined.  The charge-exchange pair is the one rate this suite
      ! needs nonzero: it is what moves the proton row after the molecular
      ! rows have split it.
      integer, intent(in) :: j
      bg_cell(j)%P_HI        = 1.0d-6      ! [1/s]
      bg_cell(j)%P_HeI       = 1.0d-7
      bg_cell(j)%P_HeII      = 1.0d-9
      bg_cell(j)%P_HeITR     = 1.0d-6
      bg_cell(j)%P_H2        = 1.0d-8
      bg_cell(j)%P_H2_di     = 0.0d0
      bg_cell(j)%P_H2_dd     = 0.0d0
      bg_cell(j)%P_H2_nd     = 0.0d0
      bg_cell(j)%k_LW        = 0.0d0
      bg_cell(j)%rchiiB      = 2.6d-13     ! [cm^3/s] at ~1e4 K
      bg_cell(j)%rcheiiB     = 4.3d-13
      bg_cell(j)%rcheiiiB    = 2.2d-12
      bg_cell(j)%rcheiTR     = 1.0d-13
      bg_cell(j)%a_ion_HI    = 0.0d0
      bg_cell(j)%a_ion_HeI   = 0.0d0
      bg_cell(j)%a_ion_HeII  = 0.0d0
      bg_cell(j)%a_ion_HeITR = 0.0d0
      bg_cell(j)%q13         = 0.0d0
      bg_cell(j)%q31a        = 0.0d0
      bg_cell(j)%q31b        = 0.0d0
      bg_cell(j)%Q31         = 0.0d0
      bg_cell(j)%A31         = 1.272d-4    ! [1/s], He 2^3S radiative decay
      bg_cell(j)%kcx_He0_Hp  = 1.0d-14     ! [cm^3/s]
      bg_cell(j)%kcx_Hep_H0  = 1.0d-9
      bg_cell(j)%nh          = 0.0d0
      bg_cell(j)%nhe         = 0.0d0
      bg_cell(j)%n_ofam      = 0.0d0
      bg_cell(j)%n_co        = 0.0d0
      bg_cell(j)%x_h2_fixed  = .false.
      bg_cell(j)%x_ox_fixed  = .false.
      bg_cell(j)%x_hp_fixed  = .false.
      bg_cell(j)%x_h2_fix    = 0.0d0
      bg_cell(j)%x_oh_fix    = 0.0d0
      bg_cell(j)%x_h2o_fix   = 0.0d0
      bg_cell(j)%x_hp_fix    = 0.0d0
      bg_cell(j)%T_K         = T_gas
      bg_cell(j)%ntot        = 0.5d0*rho(j)*n0
      end subroutine set_frozen_cell_rates

      end program carrier_helium_inventory

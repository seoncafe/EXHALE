      program transported_proton_constraint_layout
      ! The layout of the constrained molecular network: which species are
      ! unknowns of the continuation, and whether the system it builds is
      ! square.
      !
      ! Production routine exercised: constrained_network_layout_of_cell
      ! (src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90),
      ! which runs the production set_molecular_network_layout and
      ! set_rung_partition on the cell standing in ieq_cell and reports what
      ! they built.
      !
      ! REFERENCE.  ion_cell_state carries three flags saying that a
      ! partition of the cell is owned from outside it: x_h2_fixed (the
      ! lower-boundary reservoir), x_ox_fixed (carrier transport of OH and
      ! H2O) and x_hp_fixed (the transported ionization state, the key
      ! `Ionization transport`).  For each of them the same statement has to
      ! hold: the species is handed its density and is therefore NOT an
      ! unknown, and its own balance row is not a row of the system.  The
      ! fraction system the continuation's candidate is finally judged
      ! against replaces the H+ balance row by the constraint x(1) = x_hp_fix
      ! (impose_transported_ionization_fractions), so a continuation that
      ! solves the local H+ balance instead returns a composition measured
      ! against a row it was never asked to satisfy.
      !
      ! Squareness is the second statement: unknowns = reaction rows +
      ! conservation rows.  Removing an unknown without removing its row, or
      ! the other way round, leaves a system that cannot be solved at all.
      !
      ! The tolerances are exact (integer counts, and a flag read as 0 or 1).
      use global_parameters,  only: thereis_HeITR, thereis_metals,        &
                                    thereis_oxychem
      use ion_cell_state,     only: ieq_cell
      use constrained_chemical_equilibrium, only:                         &
           constrained_network_layout_of_cell, is_HII, is_H2, n_species_max
      use assertion_report
      implicit none

      ! The H/He block without the metastable, the metals or the oxygen
      ! carriers: seven balance rows (H+, He+, He++, H2, H2+, H3+, HeH+) and
      ! two conservation rows (H, He).
      integer, parameter :: n_row_HHe = 7
      real*8  :: sden(n_species_max)
      logical :: unknown_local(n_species_max), unknown_pinned(n_species_max)
      logical :: unknown_h2_pinned(n_species_max)
      integer :: nu_local, nrow_local, ncons_local
      integer :: nu_pinned, nrow_pinned, ncons_pinned
      integer :: nu_h2, nrow_h2, ncons_h2

      thereis_HeITR  = .false.
      thereis_metals = .false.
      thereis_oxychem = .false.

      ieq_cell%nh     = 1.0d10
      ieq_cell%nhe    = 8.0d8
      ieq_cell%n_ofam = 0.0d0
      ieq_cell%T_K    = 1.5d3
      ieq_cell%x_ox_fixed = .false.

      ! A composition far above the trace-holdout threshold (eps_hold = 1e-10
      ! of the element), so that the partition holds nothing out and what is
      ! measured is the cell layout alone.
      sden(:) = 1.0d8

      ! --- the H+ partition solved locally ---
      ieq_cell%x_h2_fixed = .false.
      ieq_cell%x_hp_fixed = .false.
      call constrained_network_layout_of_cell(n_row_HHe, n_row_HHe+1,     &
           n_row_HHe+1, sden, nu_local, nrow_local, ncons_local,          &
           unknown_local)

      ! --- the H+ partition imposed by the transported ionization state ---
      ieq_cell%x_hp_fixed = .true.
      call constrained_network_layout_of_cell(n_row_HHe, n_row_HHe+1,     &
           n_row_HHe+1, sden, nu_pinned, nrow_pinned, ncons_pinned,       &
           unknown_pinned)

      ! --- the H2 partition imposed by the lower-boundary reservoir: the
      !     control, the case the code already treats this way ---
      ieq_cell%x_hp_fixed = .false.
      ieq_cell%x_h2_fixed = .true.
      call constrained_network_layout_of_cell(n_row_HHe, n_row_HHe+1,     &
           n_row_HHe+1, sden, nu_h2, nrow_h2, ncons_h2, unknown_h2_pinned)

      call check_absolute('local_proton_is_an_unknown',                   &
           flag(unknown_local(is_HII)), 1.0d0, 0.0d0)
      call check_absolute('transported_proton_is_not_an_unknown',         &
           flag(unknown_pinned(is_HII)), 0.0d0, 0.0d0)
      call check_absolute('reservoir_h2_is_not_an_unknown',               &
           flag(unknown_h2_pinned(is_H2)), 0.0d0, 0.0d0)

      call check_absolute('transported_proton_unknown_count',             &
           dble(nu_pinned), dble(nu_local - 1), 0.0d0)
      call check_absolute('transported_proton_reaction_row_count',        &
           dble(nrow_pinned), dble(nrow_local - 1), 0.0d0)
      call check_absolute('transported_proton_conservation_row_count',    &
           dble(ncons_pinned), dble(ncons_local), 0.0d0)

      call check_absolute('local_proton_layout_is_square',                &
           dble(nrow_local + ncons_local - nu_local), 0.0d0, 0.0d0)
      call check_absolute('transported_proton_layout_is_square',          &
           dble(nrow_pinned + ncons_pinned - nu_pinned), 0.0d0, 0.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'transported_proton_constraint_layout: ',    &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)')                                                      &
           'transported_proton_constraint_layout: all assertions passed'

      contains

      double precision function flag(l)
      ! A logical read as a number, so that it can be compared like any
      ! other measured quantity.
      logical, intent(in) :: l
      if (l) then
         flag = 1.0d0
      else
         flag = 0.0d0
      endif
      end function flag

      end program transported_proton_constraint_layout

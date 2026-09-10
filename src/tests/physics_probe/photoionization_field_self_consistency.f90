      program photoionization_field_self_consistency
      ! QUANTITY UNDER TEST
      !   the field ONE CELL is solved in, against the composition that
      !   cell returns (ionization_equilibrium.f90,
      !   xuv_self_field_passes; util_ion_eq.f90,
      !   photoionization_field_at_cell_H / _HHe).
      !
      !   A cell's attenuated field is exp(-tau_out) times the cell mean of
      !   its OWN attenuation, and that own depth,
      !   dtau = sum_abs sigma_nu n_abs dr, is built from the very densities
      !   the cell is being solved for.  One field build and one solve
      !   therefore return a composition whose own opacity is not the
      !   opacity the field was built with.  Four things are asserted.
      !
      !     (1) THE FIELD OF A CELL IS THE B4-6 CELL MEAN OF ITS OWN
      !         COLUMN, and nothing else: with no column outside it, the
      !         rate the production routine returns is the unattenuated rate
      !         times (1 - exp(-dtau))/dtau, the single mean of
      !         utils::cell_mean_attenuation.
      !
      !     (2) NO SECOND ATTENUATION IS COMPOSED ON IT.  With a column
      !         outside as well, the rate is the unattenuated one times
      !         exp(-tau_out) (1 - exp(-dtau))/dtau -- the ONE expression
      !         the routine forms from the two depths together, not a
      !         face value multiplied by a cell mean.
      !
      !     (3) THE FIELD IS A FUNCTION OF THE STATE AND OF NOTHING ELSE.
      !         Two evaluations at one composition return one field, bit for
      !         bit, so an iterating caller is iterating a map and not a
      !         history.
      !
      !     (4) ONE PASS IS NOT SELF-CONSISTENT, AND THE ITERATION CLOSES
      !         IT.  Solving the cell's chemistry in the field of the
      !         composition it was HANDED leaves a defect: re-forming the
      !         field from what it returned and solving again moves the
      !         answer.  The defect after one pass is measured and asserted
      !         to be real; the defect after the passes this iteration takes
      !         is asserted to be at round-off, and the contraction is
      !         measured.  This is the statement xuv_self_field_passes
      !         rests on, and the statement that its default of one leaves
      !         open.
      !
      !   THE FIXED POINT IS NOT MOVED BY ANY OF IT: (5) the cell that
      !   iterates its own field to convergence and the cell that is handed
      !   the converged composition and solved once land on the same state,
      !   because at the fixed point the entry and the returned composition
      !   are one state.
      !
      ! CONFIGURATION.  The slab of xuv_cell_mean_attenuation and
      ! photoionization_field_substitution: H I alone, thereis_He false, no
      ! metals, no secondary ionization, opa_pf = 1, a_tau = 0, a
      ! MONOCHROMATIC 20 eV spectrum so that one cell's optical depth is one
      ! number and the closed forms above can be written down.  The cell is
      ! made ONE optical depth thick at the neutral density, which is where
      ! the cell mean departs from the face value by tens of per cent (B4-6)
      ! and where the self-consistency defect is largest.  The chemistry is
      ! photoionization equilibrium at a fixed recombination coefficient,
      !     P_HI (1-x) n_H = alpha n_H^2 x^2 ,
      ! solved in closed form, as in photoionization_field_substitution: it
      ! is not EXHALE's network, it is the simplest composition that
      ! responds to the field the way the real one does.

      use global_parameters
      use species_table, only: n_mion, n_mphot
      use energy_vectors_construct, only: set_energy_vectors
      use utils,        only: cell_mean_attenuation
      use utils_ion_eq, only: advance_starward_columns,                    &
                              photoionization_field_at_cell_H
      use assertion_report

      implicit none

      real*8, parameter :: alpha_B  = 2.59d-13   ! case B at 1e4 K [cm^3/s]
      ! The cell is one optical depth thick at the neutral density: the
      ! regime the cell mean was introduced for, and the regime in which a
      ! cell's own attenuation responds to its own composition at order one.
      real*8, parameter :: dtau_target = 1.0d0
      real*8, parameter :: n_H_tot = 1.0d9       ! H nuclei [cm^-3]
      ! Ionization parameter of the unattenuated cell, chosen so that the
      ! equilibrium fraction sits near a half, where dx/dP is largest.
      real*8, parameter :: ion_par = 1.0d0
      ! A column outside the cell, for assertion (2).
      real*8, parameter :: tau_out_target = 0.7d0
      real*8, parameter :: tol_fix  = 1.0d-13
      integer, parameter :: max_pass = 200

      real*8, dimension(:), allocatable :: nhi, nzero
      real*8, dimension(:,:), allocatable :: nm_zero
      real*8  :: dr_geo, sigma_cm2, alpha_use
      real*8  :: P_free, P_self, P_out, P_again
      real*8  :: h1, hh, qq, h1b, hhb, qqb
      real*8  :: dtau_cell, tau_out, N1_out
      real*8  :: x_hand, x_one, x_conv, x_check
      real*8  :: defect_one, defect_conv, contraction
      real*8  :: x_prev, move, move_prev
      integer :: jc, npass

      ! --- the run-wide input, as input_read resolves it ------------------
      is_PL_sed    = .false.
      do_read_sed  = .false.
      is_monochr   = .true.
      e_low        = 20.0d0
      LEUV         = 30.42d0
      LX           = -3.0d2
      a_orb        = 0.02544d0*AU
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0
      thereis_He   = .false.
      thereis_HeITR       = .false.
      thereis_Xray        = .false.
      thereis_lowIP_metal = .false.
      use_sec_ion    = .false.
      sec_ion_active = .false.

      call set_energy_vectors

      ! The photoabsorption cross section of the single bin [cm^2], and the
      ! cell width that makes the neutral cell one optical depth thick.
      sigma_cm2 = s_hi(1)*1.0d-18
      dr_geo    = dtau_target/(sigma_cm2*n_H_tot)

      N  = 3
      R0 = 1.0d10
      call allocate_grid_arrays
      dr_j   = dr_geo/R0
      opa_pf = 1.0d0
      jc     = 1

      allocate(nhi(1-Ng:N+Ng), nzero(1-Ng:N+Ng))
      allocate(nm_zero(1-Ng:N+Ng,n_mion))
      nzero   = 0.0d0
      nm_zero = 0.0d0

      ! ---------------------------------------------------------------- !
      ! (1) the field of a cell is the cell mean of its own column
      ! ---------------------------------------------------------------- !
      ! The unattenuated rate: an empty cell with nothing above it.
      nhi = 0.0d0
      call photoionization_field_at_cell_H(jc, 0.0d0, 0.0d0, 0.0d0,        &
               .false., P_free, h1, hh, qq)
      call check_positive('unattenuated_rate_is_positive', P_free)

      ! The same cell filled with neutral hydrogen, still nothing above it.
      nhi       = n_H_tot
      dtau_cell = sigma_cm2*n_H_tot*dr_geo
      call photoionization_field_at_cell_H(jc, 0.0d0, nhi(jc), 0.0d0,      &
               .false., P_self, h1, hh, qq)
      write(*,'(a,es13.6,a,es13.6)')                                       &
        '  DIAGNOSTIC own depth = ', dtau_cell,                            &
        ' , cell mean = ', cell_mean_attenuation(0.0d0, dtau_cell)
      call check_relative('one_cell_field_is_the_cell_mean_of_its_own_'//  &
                          'column',                                        &
                          P_self,                                          &
                          P_free*cell_mean_attenuation(0.0d0, dtau_cell),  &
                          1.0d-14)

      ! ---------------------------------------------------------------- !
      ! (2) the two depths enter ONE mean, not two attenuations
      ! ---------------------------------------------------------------- !
      tau_out = tau_out_target
      N1_out  = tau_out/sigma_cm2
      call photoionization_field_at_cell_H(jc, N1_out, nhi(jc), 0.0d0,     &
               .false., P_out, h1, hh, qq)
      call check_relative('field_is_one_cell_mean_over_both_depths',       &
                          P_out,                                           &
                          P_free*cell_mean_attenuation(tau_out,dtau_cell), &
                          1.0d-14)
      ! What a face value times a cell mean would give, i.e. the attenuation
      ! composed twice.  It must NOT be what the routine returns; the check
      ! is that the two differ by far more than the tolerance above.
      call check_at_least('a_second_attenuation_would_be_visible',         &
                          abs(P_out - P_free*exp(-tau_out)                 &
                              *cell_mean_attenuation(0.0d0,dtau_cell)      &
                              *cell_mean_attenuation(0.0d0,dtau_cell))     &
                          /P_out, 1.0d-2)

      ! ---------------------------------------------------------------- !
      ! (3) determinism: one composition, one field
      ! ---------------------------------------------------------------- !
      call photoionization_field_at_cell_H(jc, N1_out, nhi(jc), 0.0d0,     &
               .false., P_again, h1b, hhb, qqb)
      call check_absolute('two_evaluations_at_one_state_return_one_rate',  &
                          P_again, P_out, 0.0d0)
      call check_absolute('two_evaluations_at_one_state_return_one_'//     &
                          'heating_rate', h1b, h1, 0.0d0)
      call check_absolute('two_evaluations_at_one_state_return_one_'//     &
                          'efficiency', qqb, qq, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! (4) one pass is not self-consistent; the iteration closes it
      ! ---------------------------------------------------------------- !
      ! The recombination coefficient that puts the unattenuated cell at the
      ! stated ionization parameter.
      alpha_use = P_free/(n_H_tot*ion_par)

      ! The composition the cell is HANDED: fully neutral, the coldest
      ! possible entry state, so that the defect of one pass is the whole
      ! response of the field to the composition.
      x_hand = 0.0d0
      x_one  = solve_at_field_of(x_hand)
      ! The defect: what the cell's own field says its composition should be
      ! once that field is re-formed from the composition just returned.
      x_check    = solve_at_field_of(x_one)
      defect_one = abs(x_check - x_one)
      write(*,'(a,f10.6,a,f10.6,a,es10.3)')                                &
        '  DIAGNOSTIC one pass x = ', x_one, ' , second pass x = ',         &
        x_check, ' , defect = ', defect_one
      ! MEASURED on this cell: one pass leaves a defect of order a per cent
      ! of the fraction itself.  The bound is set well below the measured
      ! value and well above round-off, so it states "one pass is not
      ! self-consistent" and nothing finer.
      call check_at_least('one_pass_leaves_a_self_consistency_defect',     &
                          defect_one, 1.0d-4)

      ! The iteration the cell runs when it is given more than one pass.
      x_conv    = x_hand
      x_prev    = x_hand
      move_prev = 0.0d0
      contraction = 0.0d0
      npass     = 0
      do
         npass  = npass + 1
         x_prev = x_conv
         x_conv = solve_at_field_of(x_conv)
         move   = abs(x_conv - x_prev)
         ! The contraction of the mode, taken once the iteration is in its
         ! linear regime (the second step onward).
         if (npass .ge. 2 .and. move_prev .gt. 0.0d0)                      &
            contraction = move/move_prev
         move_prev = move
         if (move .le. tol_fix .or. npass .ge. max_pass) exit
      enddo
      defect_conv = abs(solve_at_field_of(x_conv) - x_conv)
      write(*,'(a,i0,a,es10.3,a,f8.4)')                                    &
        '  DIAGNOSTIC self-consistent in ', npass,                         &
        ' passes, remaining defect = ', defect_conv,                       &
        ' , contraction = ', contraction
      call check_at_most('self_consistent_cell_defect_is_at_round_off',    &
                         defect_conv, 1.0d-12)
      call check_at_most('self_consistent_cell_passes_are_bounded',        &
                         dble(npass), 60.0d0)
      call check_at_most('self_field_mode_contracts', contraction, 0.9d0)

      ! ---------------------------------------------------------------- !
      ! (5) the fixed point is the fixed point of the one-pass map too
      ! ---------------------------------------------------------------- !
      ! A cell handed the converged composition and solved ONCE returns it
      ! again: the two schemes differ in the path, not in the answer.
      call check_absolute('one_pass_from_the_fixed_point_returns_it',      &
                          solve_at_field_of(x_conv), x_conv, 1.0d-12)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'photoionization_field_self_consistency: ',   &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'photoionization_field_self_consistency: '//          &
                     'all assertions passed'

      contains

      !--------------!

      double precision function solve_at_field_of(x_field) result(x_new)
      ! One pass of the cell: build the cell's field from the composition
      ! x_field through the production routine, then solve the cell's
      ! chemistry in it.  The column outside the cell is held at zero, so
      ! what this map iterates is the cell's OWN optical depth and nothing
      ! else -- which is the mode xuv_self_field_passes closes.
      real*8, intent(in) :: x_field
      real*8 :: P, hone, hcell, eff, a, b
      nhi(jc) = n_H_tot*(1.0d0 - x_field)
      call photoionization_field_at_cell_H(jc, 0.0d0, nhi(jc), 0.0d0,      &
               .false., P, hone, hcell, eff)
      ! P (1-x) n_H = alpha n_H^2 x^2 , x in [0,1]
      a = alpha_use*n_H_tot
      b = P
      x_new = (-b + sqrt(b*b + 4.0d0*a*b))/(2.0d0*a)
      end function solve_at_field_of

      end program photoionization_field_self_consistency

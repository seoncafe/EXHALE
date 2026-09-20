      program element_operator_tests
      ! WHOSE BOUNDARY THE ELEMENT TRANSPORT OPERATOR IS MEASURED UNDER
      ! (PLAN_20260909_rev1 items N26b and N26c).
      !
      ! The element operator has two spellings of one balance: the fixed-wind
      ! relaxation, which takes backward-Euler steps of it, and
      ! element_transport_residual, which measures it at an infinite step
      ! length.  At the top the condition both are posed under is zero
      ! gradient: the diffusive face flux vanishes at face N by the face
      ! coefficients, and with the upper ghosts carrying the outermost
      ! physical cell's composition the advective reconstruction of the two
      ! outermost rows sees no step across the boundary, so gas crossing the
      ! top carries the column's own composition.
      !
      ! THE TWO SPELLINGS DO NOT WRITE THAT CONDITION IN THE SAME PLACE, and
      ! this suite is the standing measurement of the difference.  The
      ! relaxation writes it on its own working arrays before every pass.
      ! The residual writes nothing: it reads the ghosts out of the species
      ! vector, so its outer boundary is the CALLER's.  Both production
      ! callers do write it -- the marching path through the relaxation
      ! itself, the stationary system through
      ! write_species_rows_into_composition (steady_newton.f90) -- so in
      ! production the two agree to the last bits and no row of a state is
      ! measured under a boundary its own caller did not pose.  A caller
      ! that writes neither gets a different operator, and the rows below
      ! say by how much rather than assuming it away.
      !
      ! ROWS.
      !   * zero_gradient_upper_ghost_leaves_the_flat_column_row_at_zero:
      !     with the ghosts the production callers write, a column with no
      !     interior gradient is an exact stationary state of the operator
      !     and its two outermost rows are zero.  This is the boundary the
      !     balance is posed under, measured.
      !   * stale_upper_ghost_changes_the_element_rows: the same physical
      !     column with upper ghosts five decades away gives DIFFERENT
      !     outermost rows.  Rows N-1 and N are both measured: the WENO3
      !     reconstruction of face N-1 reaches the first upper ghost.
      !   * the_stale_ghost_is_really_stale: the two states compared do
      !     differ in the ghosts, so the row above is not vacuous.
      !   * outermost_element_row_of_a_flat_column_reads_a_stale_ghost: an
      !     exactly stationary column is not read as stationary at the top
      !     when the ghosts are not its continuation.
      !   * the_relaxation_moved_the_column: the fixed point below is not
      !     the state the relaxation started from, so the row after it is
      !     measured on a column the operator shaped.
      !   * relaxation_fixed_point_is_not_read_as_stationary_under_a_stale
      !     _ghost: the fixed point of the relaxation, with stale ghosts put
      !     back, is read as satisfied in the interior and not at the top.
      !
      ! THE SECOND THING MEASURED HERE IS CONSERVATION.  The operator hands
      ! back a species vector, and the density the hydrodynamics evolves is
      ! not among the quantities it is allowed to change: what its species
      ! weigh under the code's own mass policy (calc_rho) has to be what they
      ! weighed when it was handed them.
      !   * the_column_is_mass_closed_as_built / _with_molecules: the columns
      !     the two rows below run on carry their own density exactly, so a
      !     closure measured after the operator is the operator's.
      !   * the_relaxation_returns_the_mass_it_was_given, with the trace
      !     metals transported, with them frozen, and with a molecular
      !     hydrogen column.
      !   * the_relaxation_moved_the_composition: the fixed point is not the
      !     state it started from, so the closure rows are not vacuous.
      !   * the_advective_write_back_returns_the_mass_it_was_given: the other
      !     projection, the one that writes advected element mass fractions
      !     back into the species vector, on a state whose helium and metals
      !     have both been moved.
      !
      ! THE THIRD IS THE FAILURE CONTRACT.  The operator solves a discrete
      ! transport equation at a fixed density, and neither an unsolved
      ! nonlinear iteration nor a composition that has stopped carrying that
      ! density is visible in a finite, nonnegative species vector.  What it
      ! hands back is therefore either a solution of its own step or the
      ! composition it was given, and the rows of
      ! an_inadmissible_composition_is_not_handed_back say which by name.
      !
      ! Gravity is off, the ambipolar field is off and the hydrogen is
      ! atomic, so the settling coefficient G vanishes identically and the
      ! cell coefficients of the two outermost rows carry no ghost of their
      ! own.  What is left in those rows is the boundary and nothing else.
      !
      ! Every row prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line.  Lines beginning with two spaces or with DIAGNOSTIC are
      ! context, not verdicts.  Exit status is nonzero if any row fails.

      use global_parameters
      use species_table
      use element_inventory
      use diffusive_photochemistry, only: carrier_set_init
      use composition, only: mass_per_H_nucleus_without_He
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use grid_construction, only: define_grid,                          &
                                   spherical_face_area_and_cell_volume
      use binary_element_diffusion, only: element_transport_residual,     &
                                          element_nucleus_face_flux,      &
                                          mixture_mass_sum,               &
                                          relax_element_composition,      &
                                          element_diffusion_step,         &
                                          element_mass_fractions,         &
                                          project_element_mass_fractions, &
                                          element_step_accepted,          &
                                          element_step_solve_failed,      &
                                          element_step_mass_closure_failed,&
                                          element_relaxation_converged,   &
                                          element_mass_closure_departure, &
                                          element_transport_residual_norm,&
                                          element_row_scale_floor,        &
                                          element_step_last_status,       &
                                          element_base_flux,              &
                                          element_base_flux_report
      ! The certification's own reduction of a row, so that the number the
      ! progress control reads and the number the certification reports for
      ! the same rows are asserted to be one measure and not two.
      use certification, only: certification_row_measure, cert_scale_floor
      use test_columns, only: column_carrying_its_own_density,            &
                              column_mass_closure
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      implicit none
      integer :: nf
      nf = 0
      ! First, because one of its rows is the state of the base flux record
      ! before any operator of this program has written one.
      call the_base_flux_record_belongs_to_its_evaluation(nf)
      call the_upper_ghost_of_the_residual_is_the_callers(nf)
      call the_relaxation_fixed_point_under_a_callers_ghost(nf)
      call the_element_transport_returns_the_mass_it_was_given(nf)
      call the_advective_write_back_returns_the_mass_it_was_given(nf)
      call an_inadmissible_composition_is_not_handed_back(nf)
      call the_enforcement_does_not_read_the_caller(nf)
      call the_progress_measures_of_one_element_pass(nf)
      call the_element_row_is_the_divergence_of_one_flux(nf)
      call the_spherical_divergence_is_second_order(nf)
      write(*,'(A)') ''
      if (nf .gt. 0) then
         write(*,'(A,I0,A)') 'element_operator: ', nf, ' row(s) failed'
         stop 1
      endif
      write(*,'(A)') 'element_operator: every row passed'

      contains

      ! ================================================================= !

      subroutine bound_row(name, measured, tol, nf)
      ! One row: a measured number against zero, to a stated bound.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, tol
      integer,          intent(inout) :: nf
      if (measured .le. tol) then
         write(*,'(A,A,A,ES12.5,A,ES12.5)') 'PASS ', name,                &
              ' measured=', measured, ' reference=0 tol=', tol
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5)') 'FAIL ', name,                &
              ' measured=', measured, ' reference=0 tol=', tol
         nf = nf + 1
      endif
      end subroutine bound_row

      ! ================================================================= !

      subroutine exceeds_row(name, measured, floor_v, nf)
      ! One row: a measured number that has to be LARGER than a stated
      ! floor, which is how a non-vacuity guard reads.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, floor_v
      integer,          intent(inout) :: nf
      if (measured .gt. floor_v) then
         write(*,'(A,A,A,ES12.5,A,ES12.5)') 'PASS ', name,                &
              ' measured=', measured, ' reference=>', floor_v, ' tol=0'
      else
         write(*,'(A,A,A,ES12.5,A,ES12.5)') 'FAIL ', name,                &
              ' measured=', measured, ' reference=>', floor_v, ' tol=0'
         nf = nf + 1
      endif
      end subroutine exceeds_row

      ! ================================================================= !

      subroutine outcome_row(name, measured, expected, nf)
      ! One row: a named outcome against the one the contract requires.
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: measured, expected
      integer,          intent(inout) :: nf
      if (measured .eq. expected) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=',        &
              measured, ' reference=', expected, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=',        &
              measured, ' reference=', expected, ' tol=0'
         nf = nf + 1
      endif
      end subroutine outcome_row

      ! ================================================================= !

      subroutine synthetic_element_column(nc)
      ! The configuration and the grid both rows run on: helium and trace
      ! metals transported, WENO3 reconstruction, no gravity, no ambipolar
      ! field, atomic hydrogen.
      integer, intent(in) :: nc
      integer :: j
      real*8  :: dr_u

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      b0 = 0.0d0
      spherical_domain     = .true.
      thereis_He           = .true.
      thereis_HeITR        = .false.
      thereis_mol          = .false.
      carrier_transport    = .false.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_metals       = .true.
      eos_include_metals   = .true.
      he_diffusion         = .true.
      he_metal_diffusion   = .true.
      he_ambipolar         = .false.
      he_alphaT            = 0.0d0
      he_kzz               = 1.0d11
      HeH                  = 0.0833333333333333d0
      use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
      recon_lambda_on = .false.
      call carrier_set_init()

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      dr_u = 1.0d0/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u
      j_min = 1
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 1.0d-4
      end subroutine synthetic_element_column

      ! ================================================================= !

      subroutine production_grid_element_column(nc)
      ! The same configuration as synthetic_element_column, on the grid the
      ! PRODUCTION constructor builds: the Mixed grid of the LHS 1140 b
      ! catalog, nc cells with 50 uniform base cells of 2e-4 R_p under a
      ! geometric stretch out to 30 R_p.  The cell width then varies by four
      ! decades down the column, which is where a divergence weighted by
      ! r_j^2 dr_j and one weighted by the exact shell volume stop agreeing
      ! from one cell to the next.
      !
      ! A temperature gradient and a thermal-diffusion coefficient are put
      ! in so that the settling drift B is not identically zero and the
      ! Peclet switch has both branches to choose between; gravity and the
      ! ambipolar field stay off, so G is the thermal term alone.
      integer, intent(in) :: nc

      N  = nc
      T0 = 1.0d3
      R0 = 1.0d10
      n0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      b0 = 0.0d0
      spherical_domain     = .true.
      thereis_He           = .true.
      thereis_HeITR        = .false.
      thereis_mol          = .false.
      carrier_transport    = .false.
      thereis_oxychem      = .false.
      ionization_transport = .false.
      thereis_metals       = .true.
      eos_include_metals   = .true.
      he_diffusion         = .true.
      he_metal_diffusion   = .true.
      he_ambipolar         = .false.
      he_alphaT            = 3.0d-1
      he_kzz               = 1.0d11
      HeH                  = 0.0833333333333333d0
      use_plm = .false.;  use_weno3 = .true.;  rec_method = 'WENO3'
      recon_lambda_on = .false.
      call carrier_set_init()

      grid_type   = 'Mixed'
      N_low_cells = 50
      dr_base     = 2.0d-4
      r_max       = 30.0d0
      r_esc       = 2.0d0
      r_flux      = 1.2d0
      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      call define_grid
      call eddy_diffusion_on_grid
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 1.0d-4
      end subroutine production_grid_element_column

      ! ================================================================= !

      subroutine the_element_row_is_the_divergence_of_one_flux(nf)
      ! WHETHER THE ELEMENT ROW IS THE DIVERGENCE OF A SINGLE FACE FLUX
      ! (docs/PLAN_20260917.md item L30).
      !
      ! The conservative spherical divergence is
      ! [A_+ F_+ - A_- F_-]/V_j with A = r_edg^2 and V = (r_+^3 - r_-^3)/3,
      ! and every contribution to one conserved quantity has to divide by
      ! that same V_j.  When one half of the row divides by r_j^2 dr_j
      ! instead, its face terms carry the extra factor V_j/(r_j^2 dr_j),
      ! which differs between neighbours on a stretched grid: the internal
      ! faces stop cancelling and V_j times the row is no longer a
      ! difference of two face quantities at all.
      !
      ! ROWS.
      !  * element_row_is_the_divergence_of_the_exposed_flux: cell by cell,
      !    V_j R_j equals A_+ Phi_+ - A_- Phi_- built from the arrays
      !    element_nucleus_face_flux returns, so there is one spelling of
      !    the flux and a caller that needs it reads the operator's own.
      !  * element_column_sum_is_the_boundary_flux_difference: summed down
      !    the column the internal faces cancel and what is left is the
      !    difference of the two boundary face fluxes, measured against the
      !    size of the terms the rows balance (the operator's own row
      !    scale), which is the relative measure every other row of this
      !    operator is read on.
      !  * a_column_at_rest_with_uniform_composition_has_a_zero_element_row
      !    and a_constant_nucleus_flux_column_has_a_zero_element_row: the
      !    two states whose exact row is zero, so that the rows above are
      !    not measuring an operator that returns zero for everything.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 500
      real*8, dimension(:),   allocatable :: rho_c, T_c, Frho, msum
      real*8, dimension(:,:), allocatable :: f_c
      real*8, dimension(:),   allocatable :: res_he_a, sc_he_a
      real*8, dimension(:,:), allocatable :: res_tr, sc_tr
      real*8, dimension(:),   allocatable :: Fadv, Jdif, n_el, m_one
      real*8, dimension(:),   allocatable :: fa, cv
      real*8  :: famp, Kj, sL, sR, cadv, R0sq, R0cb
      real*8  :: rr, worst, wscale, csum, cscale, bnd, dvj
      logical :: ok_he, ok_tr
      integer :: j

      call production_grid_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng), Frho(1-Ng:N+Ng),         &
               msum(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(res_he_a(1:N), sc_he_a(1:N))
      allocate(res_tr(1:N,n_melem), sc_tr(1:N,n_melem))
      allocate(Fadv(0:N), Jdif(0:N), n_el(0:N), m_one(0:N))
      allocate(fa(0:N), cv(1:N))

      famp = 1.0d-4
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         T_c(j)   = 1.0d0 + 2.0d0*(r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 1.0d-1)

      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he_a,    &
               sc_he_a, ok_he, res_tr, sc_tr, ok_tr)
      call element_nucleus_face_flux(rho_c, T_c, f_c, Frho, Fadv, Jdif,   &
               n_el, m_one)
      call mixture_mass_sum(f_c, msum)
      call spherical_face_area_and_cell_volume(fa, cv)
      R0sq = R0*R0
      R0cb = R0sq*R0

      ! Cell by cell: the row rebuilt from the exposed face flux, in the
      ! order composition_residual assembles it.
      worst  = 0.0d0
      wscale = 0.0d0
      csum   = 0.0d0
      cscale = 0.0d0
      do j = 2, N
         Kj   = 1.0d0/(cv(j)*R0cb)
         sL   = fa(j-1)*R0sq
         sR   = fa(j)*R0sq
         cadv = n0*mu*msum(j)*v0/R0
         dvj  = (fa(j)*Fadv(j) - fa(j-1)*Fadv(j-1))/cv(j)
         rr   = Kj*(sR*Jdif(j) - sL*Jdif(j-1))
         rr   = rr + dvj*cadv
         worst  = max(worst, abs(rr - res_he_a(j)))
         wscale = max(wscale, sc_he_a(j))
         csum   = csum   + cv(j)*res_he_a(j)
         cscale = cscale + cv(j)*sc_he_a(j)
      enddo
      write(*,'(A,ES12.5,A,ES12.5)')                                      &
           '  DIAGNOSTIC largest |row - divergence of the exposed flux| ', &
           worst, ', largest row scale ', wscale
      call bound_row('element_row_is_the_divergence_of_the_exposed_'//    &
           'flux', worst/max(wscale, 1.0d-300), 1.0d-14, nf)

      ! Down the column: the internal faces cancel and the two boundary
      ! face fluxes are left.  The unit factor of the advective half is the
      ! cell's own mixture mass, so the cancellation is exact up to the
      ! variation of that mass between neighbours, which is the state's
      ! mass closure and not a property of the grid.
      bnd = (fa(N)*Jdif(N) - fa(1)*Jdif(1))/R0                            &
          + n0*mu*msum(N)*v0/R0*(fa(N)*Fadv(N) - fa(1)*Fadv(1))
      write(*,'(A,ES16.9,A,ES16.9)')                                      &
           '  DIAGNOSTIC column sum ', csum, ', boundary difference ', bnd
      write(*,'(A,ES12.5)')                                               &
           '  DIAGNOSTIC summed row scale ', cscale
      call bound_row('element_column_sum_is_the_boundary_flux_'//         &
           'difference', abs(csum - bnd)/max(cscale, 1.0d-300),           &
           1.0d-13, nf)

      ! A column at rest with no composition gradient: no face carries a
      ! flux of either kind, so every row is zero.
      Frho = 0.0d0
      do j = 1-Ng, N+Ng
         T_c(j) = 1.0d0
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .false., 0.0d0)
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he_a,    &
               sc_he_a, ok_he, res_tr, sc_tr, ok_tr)
      worst  = 0.0d0
      wscale = 0.0d0
      do j = 2, N
         worst  = max(worst, abs(res_he_a(j)))
         wscale = max(wscale, sc_he_a(j))
      enddo
      call bound_row('a_column_at_rest_with_uniform_composition_has_a_'// &
           'zero_element_row', worst/max(wscale, 1.0d-300), 1.0d-14, nf)

      ! A wind whose nucleus flux A F_rho is the same through every face,
      ! at a uniform composition: each cell gains through one face exactly
      ! what it loses through the other, so the row is zero again -- and
      ! this time the advective term of every row is large.
      do j = 1-Ng, N+Ng
         Frho(j) = famp/(r_edg(j)*r_edg(j))
      enddo
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he_a,    &
               sc_he_a, ok_he, res_tr, sc_tr, ok_tr)
      worst  = 0.0d0
      wscale = 0.0d0
      do j = 2, N
         worst  = max(worst, abs(res_he_a(j)))
         wscale = max(wscale, sc_he_a(j))
      enddo
      write(*,'(A,ES12.5,A,ES12.5)')                                      &
           '  DIAGNOSTIC constant-flux column: largest row ', worst,      &
           ', largest row scale ', wscale
      call bound_row('a_constant_nucleus_flux_column_has_a_zero_'//       &
           'element_row', worst/max(wscale, 1.0d-300), 1.0d-12, nf)

      deallocate(rho_c, T_c, Frho, msum, f_c, res_he_a, sc_he_a, res_tr,  &
                 sc_tr, Fadv, Jdif, n_el, m_one, fa, cv)
      end subroutine the_element_row_is_the_divergence_of_one_flux

      ! ================================================================= !

      subroutine the_spherical_divergence_is_second_order(nf)
      ! THE ORDER OF THE DISCRETE SPHERICAL DIVERGENCE the transport
      ! operators divide by.  For the analytic face flux
      !
      !    F(r) = sin(3 r)/r^2 ,   (1/r^2) d(r^2 F)/dr = 3 cos(3 r)/r^2 ,
      !
      ! the finite-volume divergence [A_+ F_+ - A_- F_-]/V_j is a
      ! second-order approximation of the right-hand side at the cell
      ! centre on a uniform grid.  Measured at 250, 500 and 1000 cells.
      !
      ! This row does not separate the two volume weights: they differ at
      ! O(dr^2) themselves, so an accuracy test cannot see the difference.
      ! It is the guard that the exact volume did not cost the order.
      integer, intent(inout) :: nf
      integer, parameter :: nres = 3
      integer, dimension(nres) :: ncell = (/ 250, 500, 1000 /)
      real*8, dimension(nres)  :: emax
      real*8, dimension(:), allocatable :: fa, cv
      real*8  :: dvr, ex, order_lo, order_hi
      integer :: k, j

      do k = 1, nres
         call synthetic_element_column(ncell(k))
         allocate(fa(0:N), cv(1:N))
         call spherical_face_area_and_cell_volume(fa, cv)
         emax(k) = 0.0d0
         do j = 1, N
            dvr = (fa(j)*sin(3.0d0*r_edg(j))/(r_edg(j)*r_edg(j))          &
                 - fa(j-1)*sin(3.0d0*r_edg(j-1))                          &
                   /(r_edg(j-1)*r_edg(j-1)))/cv(j)
            ex  = 3.0d0*cos(3.0d0*r(j))/(r(j)*r(j))
            emax(k) = max(emax(k), abs(dvr - ex))
         enddo
         deallocate(fa, cv)
      enddo
      order_lo = log(emax(1)/emax(2))/log(2.0d0)
      order_hi = log(emax(2)/emax(3))/log(2.0d0)
      write(*,'(A,3ES12.5)') '  DIAGNOSTIC max error at 250/500/1000: ',  &
           emax(1), emax(2), emax(3)
      write(*,'(A,2F8.4)') '  DIAGNOSTIC observed orders: ',              &
           order_lo, order_hi
      call bound_row('spherical_divergence_order_250_to_500',             &
           max(0.0d0, 2.0d0 - order_lo), 2.0d-1, nf)
      call bound_row('spherical_divergence_order_500_to_1000',            &
           max(0.0d0, 2.0d0 - order_hi), 2.0d-1, nf)
      end subroutine the_spherical_divergence_is_second_order


      ! ================================================================= !

      subroutine the_progress_measures_of_one_element_pass(nf)
      ! THE THREE QUANTITIES ONE ELEMENT PASS RETURNS ARE THREE DIFFERENT
      ! QUANTITIES (docs/PLAN_20260916_rev3.md section 3 step 1).
      !
      !   * the residual of the elemental transport balance at the
      !     composition the pass hands back, which is what the outer
      !     iteration reads as the composition's distance from the fixed
      !     point of its own operator.  The row below states that it is the
      !     certification's own measure of the same rows: the reduction
      !     certification_row_measure forms, on the rows
      !     element_transport_residual returns for the same state.  Neither
      !     spelling may drift from the other, because the outer iteration
      !     judges progress by one of them and the acceptance by the other.
      !   * the distance to the endpoint of the relaxation's inner map,
      !     formed before the damping.
      !   * the displacement the pass actually applied, which the damping
      !     shortens: with omega = 1 the two coincide exactly, and with
      !     omega = 1/2 on a column whose only transported element is
      !     helium the displacement is half the map distance, because the
      !     damping is the only thing standing between the two.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho
      real*8, dimension(:,:), allocatable :: f_c
      real*8, dimension(:),   allocatable :: res_he_a, sc_he_a
      real*8, dimension(:,:), allocatable :: res_tr, sc_tr
      logical, dimension(n_melem) :: carried
      real*8  :: famp, map_d, disp, rnorm, res_abs, sc_abs
      real*8  :: rmax_cert, rmax_one, ratio
      logical :: ok_he, ok_tr, ok_fin, measured
      integer :: j, im, nstep, jworst, ielem, jw

      call synthetic_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(res_he_a(1:N), sc_he_a(1:N))
      allocate(res_tr(1:N,n_melem), sc_tr(1:N,n_melem))

      famp = 1.0d-4
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)

      ! ONE UNDAMPED PASS, and the state it hands back.
      call relax_element_composition(rho_c, v_c, T_c, f_c, Frho, 1.0d0,   &
                                     map_d, nstep, displacement = disp)
      write(*,'(A,I0,A,ES12.5,A,ES12.5)')                                 &
           '  DIAGNOSTIC undamped pass: steps = ', nstep,                 &
           ', map distance = ', map_d, ', displacement = ', disp
      call bound_row('an_undamped_pass_displacement_equals_its_map_'//    &
           'distance', abs(disp - map_d)/max(map_d, 1.0d-300), 1.0d-14, nf)

      ! THE RESIDUAL OF THAT STATE, BOTH WAYS.  The reference is the
      ! certification's reduction of the same rows, element by element.
      call element_transport_residual_norm(rho_c, T_c, f_c, Frho, rnorm,  &
               jworst, ielem, res_abs, sc_abs, measured)
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he_a,    &
               sc_he_a, ok_he, res_tr, sc_tr, ok_tr, tr_carried = carried)
      call certification_row_measure(N, j_min, res_he_a, sc_he_a,         &
                                     rmax_cert, jw, ok_fin)
      do im = 1, n_melem
         if (.not. carried(im)) cycle
         call certification_row_measure(N, j_min, res_tr(:,im),           &
                                        sc_tr(:,im), rmax_one, jw, ok_fin)
         rmax_cert = max(rmax_cert, rmax_one)
      enddo
      write(*,'(A,ES12.5,A,ES12.5,A,I0,A,I0)')                            &
           '  DIAGNOSTIC element residual norm = ', rnorm,                &
           ', certification row measure = ', rmax_cert, ' at cell ',      &
           jworst, ', element ', ielem
      call outcome_row('the_element_residual_of_this_state_is_measured',  &
           merge(1, 0, measured), 1, nf)
      call bound_row('element_residual_norm_is_the_certification_row_'//  &
           'measure', abs(rnorm - rmax_cert)/max(rmax_cert, 1.0d-300),    &
           1.0d-14, nf)
      ! The two floors are one floor: the certification divides by
      ! cert_scale_floor and this module by element_row_scale_floor, and a
      ! row whose terms vanish must read the same number in both.
      call bound_row('the_two_row_scale_floors_are_one_number',           &
           abs(element_row_scale_floor - cert_scale_floor)                &
           /cert_scale_floor, 0.0d0, nf)
      ! The row above is not vacuous: the state carries a residual the
      ! measure can tell from zero.
      call exceeds_row('the_relaxed_state_still_carries_an_element_row',  &
           rnorm, 0.0d0, nf)

      ! A DAMPED PASS, on a column whose only transported element is
      ! helium, so the measure runs over the one element the damping acts
      ! on.  The relaxation is entered again from the state above, which is
      ! not its own fixed point at the boundary the caller left.
      he_metal_diffusion = .false.
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)
      call relax_element_composition(rho_c, v_c, T_c, f_c, Frho, 5.0d-1,  &
                                     map_d, nstep, displacement = disp)
      ratio = disp/max(map_d, 1.0d-300)
      write(*,'(A,ES12.5,A,ES12.5,A,F10.6)')                              &
           '  DIAGNOSTIC damped pass: map distance = ', map_d,            &
           ', displacement = ', disp, ', ratio = ', ratio
      call bound_row('a_damped_pass_displacement_is_omega_times_its_'//   &
           'map_distance', abs(ratio - 5.0d-1)/5.0d-1, 1.0d-10, nf)
      he_metal_diffusion = .true.

      deallocate(rho_c, v_c, T_c, Frho, f_c, res_he_a, sc_he_a, res_tr,   &
                 sc_tr)
      end subroutine the_progress_measures_of_one_element_pass

      ! ================================================================= !

      subroutine put_the_stale_ghosts_back(f_c, mpH)
      ! The upper ghosts a caller that never wrote them leaves behind: five
      ! decades of metal and three times the helium of the cell they close,
      ! the distance the atomic element reload measured (item N25).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_c
      real*8,                                 intent(in)    :: mpH
      integer :: j, ie
      do j = N+1, N+Ng
         f_c(j,:)       = 0.0d0
         f_c(j,isp_HI)  = 1.0d0/mpH
         f_c(j,isp_HeI) = 3.0d0*HeH/mpH
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = 1.0d5*melem_ab(ie)/mpH
         enddo
      enddo
      end subroutine put_the_stale_ghosts_back

      ! ================================================================= !

      subroutine copy_cell_n_into_the_ghosts(f_c)
      ! The zero-gradient upper ghosts the operator poses, written into the
      ! species vector itself so that the caller hands them over.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_c
      integer :: j
      do j = N+1, N+Ng
         f_c(j,:) = f_c(N,:)
      enddo
      end subroutine copy_cell_n_into_the_ghosts

      ! ================================================================= !

      subroutine outer_rows(rho_c, T_c, f_c, Frho, rows)
      ! The two outermost element rows of a state, each normalized by its
      ! own scale: helium and every trace element at cells N-1 and N.  The
      ! WENO3 face of cell N-1 reaches the first upper ghost, so both rows
      ! are boundary rows.
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: rho_c, T_c
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_c
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Frho
      real*8, dimension(2*(1+n_melem)),       intent(out) :: rows
      real*8, dimension(1:N)         :: res_he, sc_he
      real*8, dimension(1:N,n_melem) :: res_tr, sc_tr
      logical :: ok_he, ok_tr
      integer :: im, k
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he,      &
                                      sc_he, ok_he, res_tr, sc_tr, ok_tr)
      rows = 0.0d0
      rows(1) = res_he(N-1)/max(sc_he(N-1), 1.0d-300)
      rows(2) = res_he(N)  /max(sc_he(N),   1.0d-300)
      k = 2
      do im = 1, n_melem
         rows(k+1) = res_tr(N-1,im)/max(sc_tr(N-1,im), 1.0d-300)
         rows(k+2) = res_tr(N,im)  /max(sc_tr(N,im),   1.0d-300)
         k = k + 2
      enddo
      end subroutine outer_rows

      ! ================================================================= !

      subroutine the_upper_ghost_of_the_residual_is_the_callers(nf)
      ! THE RESIDUAL'S OUTER BOUNDARY IS THE CALLER'S, MEASURED: with the
      ! ghosts the production callers write, the flat column's outermost
      ! rows are zero; with ghosts that are not the column's continuation,
      ! they are not.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho
      real*8, dimension(:,:), allocatable :: f_c
      real*8, dimension(:),   allocatable :: rows_flat, rows_stale
      real*8  :: m_1, mpH, famp, dmax, dghost, flat_row, zero_grad_row
      integer :: j, ie, k, nr

      call synthetic_element_column(nc)
      nr = 2*(1 + n_melem)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(rows_flat(nr), rows_stale(nr))

      ! r^2 F_rho constant, the wind a steady mass row carries: with a flat
      ! composition the material divergence is then exactly zero cell by
      ! cell and the balance the row measures is empty.
      famp = 1.0d-4
      m_1  = mass_per_H_nucleus_without_He()
      mpH  = m_1 + m_He_over_m_H*HeH
      f_c  = 0.0d0
      do j = 1-Ng, N+Ng
         f_c(j,isp_HI)  = 1.0d0/mpH
         f_c(j,isp_HeI) = HeH/mpH
         do ie = 1, n_melem
            f_c(j,mion_fsp(melem_i0(ie))) = melem_ab(ie)*f_c(j,isp_HI)
         enddo
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
      enddo

      call copy_cell_n_into_the_ghosts(f_c)
      call outer_rows(rho_c, T_c, f_c, Frho, rows_flat)

      ! The boundary the balance is posed under, and the one both production
      ! callers hand over: an exactly stationary column reads zero at the
      ! top as it does inside.
      zero_grad_row = 0.0d0
      do k = 1, nr
         zero_grad_row = max(zero_grad_row, abs(rows_flat(k)))
      enddo
      call bound_row(                                                     &
           'zero_gradient_upper_ghost_leaves_the_flat_column_row_at_zero',&
           zero_grad_row, 1.0d-12, nf)

      ! The distance between the two ghosts this row compares, in the
      ! helium mass fraction the operator transports: the guard that the
      ! comparison below is between two different inputs.
      dghost = 0.0d0
      call put_the_stale_ghosts_back(f_c, mpH)
      do j = N+1, N+Ng
         dghost = max(dghost, abs(f_c(j,isp_HeI) - f_c(N,isp_HeI))        &
                              /max(abs(f_c(N,isp_HeI)), 1.0d-300))
      enddo
      call outer_rows(rho_c, T_c, f_c, Frho, rows_stale)

      dmax = 0.0d0
      do k = 1, nr
         dmax = max(dmax, abs(rows_stale(k) - rows_flat(k)))
      enddo
      ! The size of the caller dependence, in the rows' own measure.
      call exceeds_row('stale_upper_ghost_changes_the_element_rows',      &
                       dmax, 1.0d-3, nf)
      call exceeds_row('the_stale_ghost_is_really_stale', dghost, 1.0d0,  &
                       nf)

      flat_row = 0.0d0
      do k = 1, nr
         flat_row = max(flat_row, abs(rows_stale(k)))
      enddo
      call exceeds_row(                                                   &
           'outermost_element_row_of_a_flat_column_reads_a_stale_ghost',  &
           flat_row, 1.0d-3, nf)

      deallocate(rho_c, v_c, T_c, Frho, f_c, rows_flat, rows_stale)
      end subroutine the_upper_ghost_of_the_residual_is_the_callers

      ! ================================================================= !

      subroutine the_relaxation_fixed_point_under_a_callers_ghost(nf)
      ! THE FIXED POINT OF THE RELAXATION IS THE ZERO OF THE ROW INSIDE, AND
      ! AT THE OUTERMOST CELL ONLY WHILE THE CALLER'S GHOSTS ARE THE
      ! COLUMN'S CONTINUATION.
      !
      ! The column starts with an interior composition gradient, so the
      ! operator has work to do and the fixed point is not the initial
      ! state.  After the relaxation the ghosts are overwritten with values
      ! that are not that column's, and the outermost row then reads a face
      ! flux the column does not carry while the interior stays at the
      ! relaxation's own stopping measure.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho
      real*8, dimension(:,:), allocatable :: f_c, Y0, Y1
      real*8, dimension(:),   allocatable :: res_he_a, sc_he_a
      real*8, dimension(:,:), allocatable :: res_tr, sc_tr
      real*8  :: m_1, mpH, famp, drift, moved
      real*8  :: row_top, row_interior
      logical :: ok_he, ok_tr
      integer :: j, im, nstep

      call synthetic_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(Y0(1-Ng:N+Ng,1+n_melem), Y1(1-Ng:N+Ng,1+n_melem))
      allocate(res_he_a(1:N), sc_he_a(1:N))
      allocate(res_tr(1:N,n_melem), sc_tr(1:N,n_melem))

      famp = 1.0d-4
      m_1  = mass_per_H_nucleus_without_He()
      mpH  = m_1 + m_He_over_m_H*HeH
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
      enddo
      ! An interior gradient in the transported elements, in a column whose
      ! species carry its own density: helium falls by a factor two from the
      ! base to the outermost cell, so the operator has work to do, and the
      ! composition it is handed is one the transport is allowed to advance
      ! (a column that already misses rho is refused before the boundary
      ! this row is about can be reached: the mass fractions this row was
      ! first written with miss it by 1.321e-1 at the outermost cell,
      ! MEASURED, because they carry helium and the metals at a radial
      ! factor the hydrogen they are taken against does not).
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)

      call element_mass_fractions(f_c, Y0)
      call relax_element_composition(rho_c, v_c, T_c, f_c, Frho, 1.0d0,   &
                                     drift, nstep)
      call element_mass_fractions(f_c, Y1)
      moved = 0.0d0
      do j = 2, N
         moved = max(moved, abs(Y1(j,1) - Y0(j,1))/max(Y0(1,1), 1.0d-300))
      enddo
      write(*,'(A,I0,A,ES12.5,A,ES12.5)')                                 &
           '  DIAGNOSTIC relaxation steps = ', nstep,                     &
           ', drift = ', drift, ', helium moved by ', moved
      ! The fixed point is not the state it started from, so the row below
      ! is measured on a column the operator actually shaped.
      call exceeds_row('the_relaxation_moved_the_column', moved, 1.0d-3,  &
                       nf)

      ! The caller's ghosts, put back after the relaxation wrote its own.
      call put_the_stale_ghosts_back(f_c, mpH)
      call element_transport_residual(rho_c, T_c, f_c, Frho, res_he_a,    &
                                      sc_he_a, ok_he, res_tr, sc_tr,      &
                                      ok_tr)
      row_top      = abs(res_he_a(N))/max(sc_he_a(N), 1.0d-300)
      row_interior = 0.0d0
      do j = 2, N-1
         row_interior = max(row_interior,                                 &
                            abs(res_he_a(j))/max(sc_he_a(j), 1.0d-300))
      enddo
      do im = 1, n_melem
         row_top = max(row_top, abs(res_tr(N,im))                         &
                                /max(sc_tr(N,im), 1.0d-300))
         do j = 2, N-1
            row_interior = max(row_interior, abs(res_tr(j,im))            &
                                             /max(sc_tr(j,im), 1.0d-300))
         enddo
      enddo
      write(*,'(A,ES12.5)') '  DIAGNOSTIC worst interior element row = ', &
           row_interior
      ! The floor is six decades above the relaxation's own stopping measure
      ! (1e-12 of the reservoir composition, which the interior row reads),
      ! so what this row measures is the boundary and not the fixed point's
      ! own accuracy.
      call exceeds_row(                                                   &
           'relaxation_fixed_point_is_not_read_as_stationary_under_a'//   &
           '_stale_ghost', row_top, 1.0d-6, nf)

      deallocate(rho_c, v_c, T_c, Frho, f_c, Y0, Y1, res_he_a, sc_he_a,   &
                 res_tr, sc_tr)
      end subroutine the_relaxation_fixed_point_under_a_callers_ghost

      ! ================================================================= !

      subroutine the_base_flux_record_belongs_to_its_evaluation(nf)
      ! THE BASE BUDGET ENTRY IS A STATEMENT ABOUT ONE EVALUATION.  The
      ! operator leaves the helium element mass flux the base carries in
      ! element_base_flux, and a report written later is about the state
      ! that wrote it and no other, so the record carries the evaluation it
      ! came from and says so when none has run.
      !
      ! What the rows assert, on a column with a nonzero bulk velocity so
      ! that the advective half is not trivially zero:
      !   * before any evaluation the record is absent, and the report says
      !     so instead of printing a flux of zero;
      !   * one step later the record is that step's, and a second step
      !     moves the count, so a caller can tell the two apart;
      !   * the diffusive half at the base face is exactly zero, which is a
      !     property of the discretization (the gradient and drift
      !     coefficients of faces 0 and N are zero) and not of the state,
      !     while at the first solved face it is the flux the operator's own
      !     face loop returned there;
      !   * the advective half is the face mass flux of the state carrying
      !     the face helium mass fraction, the same number the one public
      !     spelling of the element flux returns for that face.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, dt_c
      real*8, dimension(:),   allocatable :: Frho_f, rho_phys
      real*8, dimension(:,:), allocatable :: f_c
      real*8, dimension(:),   allocatable :: Jf_out, Fadv, Jdif, n_el, m_one
      real*8, dimension(:),   allocatable :: msum
      character(len=300) :: line
      real*8  :: dev, adv_ref
      integer :: j, jup, ev0

      call synthetic_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               dt_c(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(Frho_f(1-Ng:N+Ng), rho_phys(1-Ng:N+Ng))
      allocate(msum(1-Ng:N+Ng))
      allocate(Jf_out(0:N), Fadv(0:N), Jdif(0:N), n_el(0:N), m_one(0:N))

      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 1.0d-1
         T_c(j)   = 1.0d0
         dt_c(j)  = 1.0d-3
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)

      ! Before any evaluation of this program.
      ev0 = element_base_flux%evaluation
      call outcome_row('the_base_flux_record_starts_absent', ev0, 0, nf)
      call element_base_flux_report(line)
      call outcome_row('an_absent_base_flux_record_is_reported_as_absent',&
           min(index(line, 'no element operator evaluation'), 1), 1, nf)

      call element_diffusion_step(rho_c, v_c, T_c, f_c, dt_c,             &
                                  Jface_out=Jf_out)
      call outcome_row('the_base_flux_record_is_this_evaluations',        &
                       element_base_flux%evaluation, ev0 + 1, nf)
      call element_base_flux_report(line, ev0 + 1)
      call outcome_row('the_record_of_the_evaluation_asked_for_is_not'//  &
           '_flagged', index(line, 'NOT the evaluation asked for'), 0, nf)

      ! The base face carries no diffusive flux in this discretization, and
      ! the first solved face carries the operator's own.
      call bound_row('the_base_face_diffusive_half_is_zero',              &
                     abs(element_base_flux%diffusive_base_face),          &
                     0.0d0, nf)
      dev = abs(element_base_flux%diffusive_solved_face                   &
                - Jf_out(element_base_flux%first_solved_face))
      call bound_row('the_solved_face_diffusive_half_is_the_operators',   &
                     dev, 0.0d0, nf)

      ! The advective half against the one public spelling of the element
      ! flux, evaluated on the composition the step handed back.
      call mixture_mass_sum(f_c, msum)
      rho_phys = rho_c*n0*mu*msum
      do j = 1-Ng, N+Ng
         jup = min(j+1, N+Ng)
         Frho_f(j) = 0.5d0*(rho_phys(j) + rho_phys(jup))                  &
                    *0.5d0*(v_c(j) + v_c(jup))*v0
      enddo
      call element_nucleus_face_flux(rho_c, T_c, f_c, Frho_f,             &
                                     Fadv, Jdif, n_el, m_one)
      adv_ref = Fadv(0)
      call exceeds_row('the_base_face_advective_half_is_not_zero',        &
                       abs(element_base_flux%advective_base_face),        &
                       0.0d0, nf)
      dev = abs(element_base_flux%advective_base_face - adv_ref)          &
            /max(abs(adv_ref), 1.0d-300)
      call bound_row('the_base_face_advective_half_is_the_exposed_flux',  &
                     dev, 1.0d-12, nf)

      ! A second evaluation, and a report asked for the first one.
      call element_diffusion_step(rho_c, v_c, T_c, f_c, dt_c,             &
                                  Jface_out=Jf_out)
      call outcome_row('a_second_evaluation_moves_the_record',            &
                       element_base_flux%evaluation, ev0 + 2, nf)
      call element_base_flux_report(line, ev0 + 1)
      call outcome_row('a_record_of_another_evaluation_is_flagged',       &
           min(index(line, 'NOT the evaluation asked for'), 1), 1, nf)

      deallocate(rho_c, v_c, T_c, dt_c, f_c, Frho_f, rho_phys, msum,      &
                 Jf_out, Fadv, Jdif, n_el, m_one)
      end subroutine the_base_flux_record_belongs_to_its_evaluation

      ! ================================================================= !

      subroutine the_element_transport_returns_the_mass_it_was_given(nf)
      ! WHAT THE OPERATOR HANDS BACK WEIGHS WHAT IT WAS HANDED.  The
      ! composition it returns is its own, the density is the
      ! hydrodynamics', and the two have to describe one gas: a projection
      ! that gives hydrogen a mass the cell's metals did not release leaves
      ! a state whose species and whose rho are two different atmospheres,
      ! and no row of the stationary inventory sees it (the mass row is on
      ! rho, the element rows on ratios).
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, Frho
      real*8, dimension(:,:), allocatable :: f_c, Y0, Y1
      real*8  :: famp, cl0, cl1, moved, drift
      integer :: j, nstep, ic

      do ic = 1, 3
         call synthetic_element_column(nc)
         ! ic = 1 the trace metals transported, ic = 2 the metals frozen at
         ! their hydrogen ratio, ic = 3 a molecular hydrogen column.
         if (ic .eq. 2) he_metal_diffusion = .false.
         if (ic .eq. 3) then
            thereis_mol = .true.
            call carrier_set_init()
         endif
         allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),       &
                  Frho(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
         allocate(Y0(1-Ng:N+Ng,1+n_melem), Y1(1-Ng:N+Ng,1+n_melem))
         famp = 1.0d-4
         do j = 1-Ng, N+Ng
            rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
            v_c(j)   = 0.0d0
            T_c(j)   = 1.0d0
            Frho(j)  = famp/(r_edg(j)*r_edg(j))
         enddo
         if (ic .eq. 3) then
            call column_carrying_its_own_density(f_c, 8.0d-1,             &
                                                .true., 0.0d0)
         else
            call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)
         endif

         cl0 = column_mass_closure(rho_c, f_c)
         call element_mass_fractions(f_c, Y0)
         call relax_element_composition(rho_c, v_c, T_c, f_c, Frho,       &
                                        1.0d0, drift, nstep)
         call element_mass_fractions(f_c, Y1)
         cl1 = column_mass_closure(rho_c, f_c)
         moved = 0.0d0
         do j = 2, N
            moved = max(moved, abs(Y1(j,1) - Y0(j,1))                     &
                               /max(Y0(1,1), 1.0d-300))
         enddo

         if (ic .eq. 1) then
            call bound_row('the_column_is_mass_closed_as_built', cl0,     &
                           1.0d-14, nf)
            call exceeds_row('the_relaxation_moved_the_composition',      &
                             moved, 1.0d-3, nf)
            call bound_row(                                               &
                 'the_relaxation_returns_the_mass_it_was_given',          &
                 cl1, 1.0d-14, nf)
         else if (ic .eq. 2) then
            call bound_row(                                               &
                 'the_relaxation_returns_the_mass_it_was_given_with'//    &
                 '_the_metals_frozen', cl1, 1.0d-14, nf)
         else
            call bound_row('the_column_is_mass_closed_as_built_with'//    &
                           '_molecules', cl0, 1.0d-14, nf)
            call bound_row(                                               &
                 'the_relaxation_returns_the_mass_it_was_given_with'//    &
                 '_molecules', cl1, 1.0d-14, nf)
         endif

         deallocate(rho_c, v_c, T_c, Frho, f_c, Y0, Y1)
      enddo
      end subroutine the_element_transport_returns_the_mass_it_was_given

      ! ================================================================= !

      subroutine the_advective_write_back_returns_the_mass_it_was_given(nf)
      ! THE OTHER PROJECTION.  project_element_mass_fractions writes the
      ! element mass fractions the Runge-Kutta stages advected back into the
      ! species vector; helium and each trace element arrive on their own
      ! face fluxes, so the two do not move together and the write-back has
      ! to close the mixture around both.  The fractions here stand in for
      ! that advection: helium up by one percent, the metals down by two,
      ! with a radial shape, which is the size of the departure a marching
      ! step leaves.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c
      real*8, dimension(:,:), allocatable :: f_c, Y
      real*8  :: cl0, cl1, sh
      integer :: j, im

      call synthetic_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), f_c(1-Ng:N+Ng,n_species))
      allocate(Y(1-Ng:N+Ng,1+n_melem))
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)
      cl0 = column_mass_closure(rho_c, f_c)

      call element_mass_fractions(f_c, Y)
      do j = 1-Ng, N+Ng
         sh = (r(j) - 1.0d0)/max(r(N) - 1.0d0, 1.0d-30)
         Y(j,1) = Y(j,1)*(1.0d0 + 1.0d-2*sh)
         do im = 1, n_melem
            Y(j,1+im) = Y(j,1+im)*(1.0d0 - 2.0d-2*sh)
         enddo
      enddo
      call project_element_mass_fractions(f_c, Y)
      cl1 = column_mass_closure(rho_c, f_c)

      call bound_row('the_write_back_column_is_mass_closed_as_built',     &
                     cl0, 1.0d-14, nf)
      call bound_row(                                                     &
           'the_advective_write_back_returns_the_mass_it_was_given',      &
           cl1, 1.0d-14, nf)

      deallocate(rho_c, f_c, Y)
      end subroutine the_advective_write_back_returns_the_mass_it_was_given

      ! ================================================================= !

      subroutine an_inadmissible_composition_is_not_handed_back(nf)
      ! WHAT THE TRANSPORT HANDS BACK IS EITHER A SOLUTION OF ITS OWN STEP OR
      ! THE COMPOSITION IT WAS GIVEN.
      !
      ! The step solves a discrete transport equation at a FIXED density.
      ! Two things can be wrong with what it produces and neither is visible
      ! in a finite, nonnegative species vector: the nonlinear iteration may
      ! not have solved the row, and the species may no longer weigh the
      ! density that was held fixed.  A composition that fails either is not
      ! a state of this atmosphere, and the chemistry that reads it next has
      ! no way to tell.
      !
      ! ROWS.
      !   * a_solved_step_is_accepted / _moves_the_composition /
      !     _carries_the_density_it_was_given: the ordinary step, so that the
      !     refusals below are refusals and not a routine that never
      !     advances anything.
      !   * an_unsolved_step_is_refused / _leaves_the_entry_composition: a
      !     step whose residual is not a number at any iterate.  The line
      !     search cannot descend, so no trial along the direction is an
      !     iterate of the equation; the outcome says so and the species
      !     vector is the entry one to the last bit.
      !   * a_step_on_a_column_that_misses_rho_is_still_solved /
      !     a_step_returns_the_mass_it_was_handed: the departure from rho is
      !     an invariant of the step, so a column whose species weigh 1.2 rho
      !     is transported (the helium fraction is a ratio and does not see
      !     the scaling) and comes back weighing 1.2 rho.  Establishing the
      !     closure belongs to whoever builds the composition and preserving
      !     it to the transport; it is the increase that is refused, and the
      !     reproduced failure of D1 (0.177 of rho from an entry state closed
      !     to 1.0e-15, docs/solver_partition_experiment_20260911.md section
      !     5.3) is an increase of the whole of it.
      !   * the_relaxation_reports_the_fixed_point_it_reached /
      !     _carries_the_density_it_was_given: a relaxation that meets its
      !     movement measure reports the fixed point of the operator, which
      !     is a different outcome from running out of steps, and its
      !     returned composition closes against rho.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, T_nan, Frho
      real*8, dimension(:),   allocatable :: dt_c
      real*8, dimension(:,:), allocatable :: f_c, f_keep, Y0, Y1
      real*8  :: famp, cl0, cl1, moved, dfmax, drift
      integer :: j, st, nstep
      ! One diffusion time of the synthetic column: dr^2/K_zz = 6.6e5 s
      ! against the code time R0/v0 = 3.5e4 s (MEASURED from the column's
      ! own numbers), so a step of this length moves the composition by a
      ! finite fraction of its gradient.
      real*8, parameter :: dt_diffusive = 2.0d1

      call synthetic_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               T_nan(1-Ng:N+Ng), Frho(1-Ng:N+Ng), dt_c(1-Ng:N+Ng))
      allocate(f_c(1-Ng:N+Ng,n_species), f_keep(1-Ng:N+Ng,n_species))
      allocate(Y0(1-Ng:N+Ng,1+n_melem), Y1(1-Ng:N+Ng,1+n_melem))

      famp = 1.0d-4
      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         Frho(j)  = famp/(r_edg(j)*r_edg(j))
         dt_c(j)  = dt_diffusive
      enddo
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)
      f_keep = f_c
      cl0    = column_mass_closure(rho_c, f_c)

      ! --- the ordinary step
      call element_mass_fractions(f_c, Y0)
      call element_diffusion_step(rho_c, v_c, T_c, f_c, dt_c,             &
                                  Frho_in = Frho, status = st)
      call element_mass_fractions(f_c, Y1)
      cl1   = column_mass_closure(rho_c, f_c)
      moved = 0.0d0
      do j = 2, N
         moved = max(moved, abs(Y1(j,1) - Y0(j,1))/max(Y0(1,1), 1.0d-300))
      enddo
      call outcome_row('a_solved_step_is_accepted', st,                   &
                       element_step_accepted, nf)
      call exceeds_row('a_solved_step_moves_the_composition', moved,      &
                       1.0d-6, nf)
      call bound_row('a_solved_step_carries_the_density_it_was_given',    &
                     cl1, 1.0d-14, nf)
      write(*,'(A,ES12.5,A,ES12.5)')                                      &
           '  DIAGNOSTIC entry closure = ', cl0,                          &
           ', closure the step reported = ',                              &
           element_mass_closure_departure

      ! --- a step whose residual is not a number at any iterate
      f_c   = f_keep
      T_nan = T_c
      T_nan(N/2) = ieee_value(1.0d0, ieee_quiet_nan)
      call element_diffusion_step(rho_c, v_c, T_nan, f_c, dt_c,           &
                                  Frho_in = Frho, status = st)
      dfmax = maxval(abs(f_c - f_keep))
      call outcome_row('an_unsolved_step_is_refused', st,                 &
                       element_step_solve_failed, nf)
      call bound_row('an_unsolved_step_leaves_the_entry_composition',     &
                     dfmax, 0.0d0, nf)

      ! --- the mass the species carry is the step's invariant, measured on a
      ! column that does not carry rho at all: the departure comes back the
      ! size it went in, which is what the closure refusal is measured
      ! against and what the reproduced failure of D1 broke.
      f_c    = 1.2d0*f_keep
      f_keep = f_c
      cl0    = column_mass_closure(rho_c, f_c)
      call element_diffusion_step(rho_c, v_c, T_c, f_c, dt_c,             &
                                  Frho_in = Frho, status = st)
      cl1 = column_mass_closure(rho_c, f_c)
      call outcome_row('a_step_on_a_column_that_misses_rho_is_still'//    &
                       '_solved', st, element_step_accepted, nf)
      call bound_row('a_step_returns_the_mass_it_was_handed',             &
                     abs(cl1 - cl0), 1.0d-14, nf)
      write(*,'(A,ES12.5,A,ES12.5)') '  DIAGNOSTIC departure in = ',      &
           cl0, ', out = ', cl1

      ! --- the relaxation to the fixed point
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)
      call element_mass_fractions(f_c, Y0)
      call relax_element_composition(rho_c, v_c, T_c, f_c, Frho, 1.0d0,   &
                                     drift, nstep, status = st)
      call element_mass_fractions(f_c, Y1)
      cl1   = column_mass_closure(rho_c, f_c)
      moved = 0.0d0
      do j = 2, N
         moved = max(moved, abs(Y1(j,1) - Y0(j,1))/max(Y0(1,1), 1.0d-300))
      enddo
      write(*,'(A,I0,A,ES12.5,A,ES12.5)')                                 &
           '  DIAGNOSTIC relaxation steps = ', nstep, ', drift = ',       &
           drift, ', helium moved by ', moved
      call outcome_row('the_relaxation_reports_the_fixed_point_it'//      &
                       '_reached', st, element_relaxation_converged, nf)
      call exceeds_row('the_relaxation_of_this_column_moved_it', moved,   &
                       1.0d-3, nf)
      call bound_row('the_relaxation_carries_the_density_it_was_given',   &
                     cl1, 1.0d-14, nf)

      deallocate(rho_c, v_c, T_c, T_nan, Frho, dt_c, f_c, f_keep, Y0, Y1)
      end subroutine an_inadmissible_composition_is_not_handed_back

      ! ================================================================= !

      subroutine the_enforcement_does_not_read_the_caller(nf)
      ! THE SAME FOUR TESTS RUN WHETHER OR NOT status IS PRESENT.  R4 of
      ! docs/Update_EXHALE_stage2_review_20260912.md: the marching call at
      ! EXHALE_main.f90:2297 passes no status argument, and on the entry
      ! text (before this item) that call carried judged = .false., so an
      ! unconverged inner solve was adopted there uninspected.  The two rows
      ! below repeat the ordinary step and the residual-not-a-number step of
      ! an_inadmissible_composition_is_not_handed_back with NO status
      ! argument at the call site, and read element_step_last_status
      ! afterward instead: the entry text moves the composition on the
      ! second row (a NaN residual is adopted because nothing there ever
      ! judges it), so this row is RED before the fix and GREEN after.
      !
      ! A candidate that the converged solve itself pushes outside [0,1] by
      ! more than element_fraction_bound_tol is not reproduced here: the
      ! drift flux vanishes at both ends of the composition axis in the
      ! discrete operator (module header), so a solve that reports solved
      ! leaves X inside [0,1] to round-off on every column this suite can
      ! build, and forcing it further (an oversized dt_code, a reservoir
      ! X_base outside [0,1]) only ever reaches the solve_failed branch
      ! instead, because the same nonlinearity that could push X outside
      ! the bound is what stops the line search from descending first. The
      ! P1 report's one measured out-of-bounds excursion is on the strongly
      ! nonstationary atomic background of a production reload
      ! (atomic_elem_newton, not reachable on a synthetic column of this
      ! size); element_step_out_of_bounds is exercised there, not here.
      integer, intent(inout) :: nf
      integer, parameter :: nc = 40
      real*8, dimension(:),   allocatable :: rho_c, v_c, T_c, T_nan, Frho
      real*8, dimension(:),   allocatable :: dt_c
      real*8, dimension(:,:), allocatable :: f_c, f_keep, Y0, Y1
      real*8  :: moved, dfmax
      integer :: j
      real*8, parameter :: dt_diffusive = 2.0d1

      call synthetic_element_column(nc)
      allocate(rho_c(1-Ng:N+Ng), v_c(1-Ng:N+Ng), T_c(1-Ng:N+Ng),          &
               T_nan(1-Ng:N+Ng), Frho(1-Ng:N+Ng), dt_c(1-Ng:N+Ng))
      allocate(f_c(1-Ng:N+Ng,n_species), f_keep(1-Ng:N+Ng,n_species))
      allocate(Y0(1-Ng:N+Ng,1+n_melem), Y1(1-Ng:N+Ng,1+n_melem))

      do j = 1-Ng, N+Ng
         rho_c(j) = exp(-3.0d0*(r(j) - 1.0d0))
         v_c(j)   = 0.0d0
         T_c(j)   = 1.0d0
         Frho(j)  = 1.0d-4/(r_edg(j)*r_edg(j))
         dt_c(j)  = dt_diffusive
      enddo

      ! --- a healthy step, no status: accepted, and the composition moves
      call column_carrying_its_own_density(f_c, 0.0d0, .true., 0.0d0)
      f_keep = f_c
      call element_mass_fractions(f_c, Y0)
      call element_diffusion_step(rho_c, v_c, T_c, f_c, dt_c, Frho_in = Frho)
      call element_mass_fractions(f_c, Y1)
      moved = 0.0d0
      do j = 2, N
         moved = max(moved, abs(Y1(j,1) - Y0(j,1))/max(Y0(1,1), 1.0d-300))
      enddo
      call outcome_row('an_unread_healthy_step_is_accepted',                &
                       element_step_last_status, element_step_accepted, nf)
      call exceeds_row('an_unread_healthy_step_moves_the_composition',     &
                       moved, 1.0d-6, nf)

      ! --- a residual that is not a number at any iterate, no status: the
      ! composition the caller reads back has to be the one it handed in,
      ! bit for bit, and element_step_last_status has to name the refusal.
      f_c   = f_keep
      T_nan = T_c
      T_nan(N/2) = ieee_value(1.0d0, ieee_quiet_nan)
      call element_diffusion_step(rho_c, v_c, T_nan, f_c, dt_c,            &
                                  Frho_in = Frho)
      dfmax = maxval(abs(f_c - f_keep))
      call outcome_row('an_unread_unsolved_step_is_still_refused',         &
                       element_step_last_status, element_step_solve_failed,&
                       nf)
      call bound_row(                                                      &
           'an_unread_unsolved_step_still_leaves_the_entry_composition',   &
           dfmax, 0.0d0, nf)

      deallocate(rho_c, v_c, T_c, T_nan, Frho, dt_c, f_c, f_keep, Y0, Y1)
      end subroutine the_enforcement_does_not_read_the_caller

      end program element_operator_tests

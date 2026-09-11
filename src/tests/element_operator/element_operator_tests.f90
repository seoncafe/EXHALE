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
      use binary_element_diffusion, only: element_transport_residual,     &
                                          relax_element_composition,      &
                                          element_diffusion_step,         &
                                          element_mass_fractions,         &
                                          project_element_mass_fractions, &
                                          element_step_accepted,          &
                                          element_step_solve_failed,      &
                                          element_step_mass_closure_failed,&
                                          element_relaxation_converged,   &
                                          element_mass_closure_departure
      use test_columns, only: column_carrying_its_own_density,            &
                              column_mass_closure
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      implicit none
      integer :: nf
      nf = 0
      call the_upper_ghost_of_the_residual_is_the_callers(nf)
      call the_relaxation_fixed_point_under_a_callers_ghost(nf)
      call the_element_transport_returns_the_mass_it_was_given(nf)
      call the_advective_write_back_returns_the_mass_it_was_given(nf)
      call an_inadmissible_composition_is_not_handed_back(nf)
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

      end program element_operator_tests

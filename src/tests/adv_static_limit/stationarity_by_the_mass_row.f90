      program stationarity_by_the_mass_row
      ! WHICH TEST DECIDES THAT THE STATE THE POST-PROCESS WAS HANDED IS
      ! STATIONARY IN A CELL.
      !
      ! Every equation the advection post-process solves is a STEADY
      ! equation integrated along the recorded flow, so a cell of the output
      ! is a correction only where the mass flux through it is stationary:
      ! without a mass source, div(rho v) = 0 is part of what stationary
      ! means, and an energy equation solved on a prescribed non-stationary
      ! mass profile returns a temperature the gas does not have.
      !
      ! The mass balance of the run is a FINITE-VOLUME one: the numerical
      ! face fluxes of the Riemann solve, differenced over the cell's own
      ! control volume and divided by the largest term the row itself
      ! contains (steady_residual.f90, mass_flux_row_scale). That ratio is
      ! the fractional change of the face mass flux across the cell, and it
      ! is what the stationary certification measures the mass row by
      ! (certification.f90, hydro_row_entry, cert_tol_mass).
      !
      ! THE CORRECTION IS FIRST ORDER IN THAT RATIO: the terms the steady
      ! equations drop are the ones the mass divergence puts into them, each
      ! of them the ratio times a term they keep. So a row at or below
      ! adv_conditional_tol is a CONDITIONAL correction accurate to that
      ! fraction of itself, and a row above it is not corrected at all. That
      ! fraction and the tolerance of a CERTIFIED stationary state
      ! (cert_tol_mass) are two different statements about one measure, and
      ! the rows below keep them apart.
      !
      ! A center-to-center difference of rho v r^2 is the same quantity in
      ! the continuum and a DIFFERENT discrete test. The rows below measure
      ! the two on one state: an analytic outflow whose center-sampled
      ! rho v r^2 is constant to round-off, so that the enthalpy-flux term
      ! ratio built on it is essentially zero and accepts every cell,
      ! carries a face-flux mass row five decades above cert_tol_mass. Its
      ! discrete departure from stationary is real, and what it costs the
      ! correction is that departure itself: two parts in ten million.
      !
      ! Everything measured here is PRODUCTION: define_grid,
      ! set_gravity_grid, set_base_reservoir, Apply_BC, assemble_residual
      ! (Reconstruct + Num_flux + RK_rhs + Source) and residual_row_scale.
      ! Nothing re-implements a flux, a source or a scale.
      !
      ! One line per assertion:
      !   PASS|FAIL <name> measured= reference= tol=
      use global_parameters
      use grid_construction,          only: define_grid
      use gravity_grid_construction,  only: set_gravity_grid
      use base_boundary,              only: set_base_reservoir
      use BC_Apply,                   only: Apply_BC
      use steady_residual_mod,        only: assemble_residual,            &
                                            residual_row_scale
      use certification,              only: cert_tol_mass, cert_scale_floor
      use post_processing,            only: enthalpy_flux_term_ratio,     &
                                            enthalpy_ratio_report_level
      use output_write,               only: adv_conditional_tol
      implicit none

      ! A small spherical domain with a uniform grid, so that the only
      ! departure from a discrete steady state is the operator's own
      ! truncation error on the analytic profile.
      integer, parameter :: Ncell = 200
      real*8,  parameter :: nhat  = 1.0d0   ! particles per unit mass
      real*8,  parameter :: v_ref = 1.0d-2  ! flux constant, rho v r^2
      real*8,  parameter :: T_base = 1.0d0  ! base temperature, code units
      real*8,  parameter :: T_grad = 0.5d0  ! dT/dr, code units
      ! Growth rate of the mass flux of the second fixture: F = rho v r^2
      ! proportional to exp(flux_growth*(r-1)), whose fractional change
      ! across a cell of width dr is flux_growth*dr in every cell, so the
      ! whole column sits at one measure by construction.
      real*8,  parameter :: flux_growth = 3.0d0

      real*8, dimension(:,:), allocatable :: u, Res
      real*8, dimension(:),   allocatable :: n_part, heat, cool
      real*8, dimension(:),   allocatable :: rho_s, v_s, p_s, Tc_s
      real*8  :: row, row_max, r_row_max, rest_row_max, row_min
      real*8  :: dlnF, dlnF_max, ratio, ratio_max
      real*8  :: Fup, Flo, dVc, divf, e_up, e_lo, w_up, h_up, drc
      integer :: j, ja, jb, nfail

      nfail = 0

      ! ---- the configuration of the state, in the keys input_read sets ----
      spherical_domain = .true.
      grid_type        = 'Uniform'
      r_max            = 3.0d0
      r_esc            = 2.0d0
      r_flux           = 1.20d0
      CFL              = 0.6d0
      rec_method       = 'WENO3'
      use_plm          = .false.
      use_weno3        = .true.
      flux             = 'ROE'
      thereis_He       = .false.
      N                = Ncell
      T0               = 1.0d3
      R0               = 1.0d10
      Mp               = 1.0d30
      v0               = sqrt(kb_erg*T0/mu)
      n0               = 1.0d0
      q0               = 1.0d0

      call allocate_grid_arrays
      call define_grid

      allocate(u(3,1-Ng:N+Ng), Res(3,1-Ng:N+Ng), n_part(1-Ng:N+Ng),         &
               heat(1-Ng:N+Ng), cool(1-Ng:N+Ng), rho_s(1-Ng:N+Ng),        &
               v_s(1-Ng:N+Ng), p_s(1-Ng:N+Ng), Tc_s(1-Ng:N+Ng))

      heat = 0.0d0
      cool = 0.0d0

      ! The interior cells, whose stencil reads no ghost row: the two tests
      ! are compared where neither of them is answering for a boundary
      ! closure.
      ja = 5
      jb = N - 4

      ! ------------------------------------------------------------------ !
      ! Row 1: A STATE AT REST IS STATIONARY. Uniform, no gravity, no flow:
      ! every face flux is zero, so the mass row is an exact zero and the
      ! decision accepts the cell. This is the row that says the operator
      ! cannot refuse a state that carries no mass flux at all.
      ! ------------------------------------------------------------------ !
      b0 = 0.0d0
      call set_gravity_grid
      do j = 1-Ng,N+Ng
         rho_s(j) = 1.0d0
         v_s(j)   = 0.0d0
         Tc_s(j)   = T_base
         p_s(j)   = nhat*rho_s(j)*Tc_s(j)
      enddo
      call build_state
      call assemble_residual(u, n_part, heat, cool, Res)
      rest_row_max = 0.0d0
      do j = ja,jb
         rest_row_max = max(rest_row_max, mass_row_of(j))
      enddo
      call chk_le('state_at_rest_is_stationary', rest_row_max,            &
                  adv_conditional_tol, nfail)

      ! ------------------------------------------------------------------ !
      ! Rows 2 to 4: THE ANALYTIC OUTFLOW.
      !
      !   rho = 1,  v = v_ref/r^2,  T = T_base + T_grad (r - 1)
      !
      ! Its center-sampled mass flux rho v r^2 is the constant v_ref to
      ! round-off, so the center-to-center divergence the enthalpy-flux term
      ! ratio is built on is essentially zero and that ratio accepts every
      ! cell (row 2). The temperature gradient keeps the advected energy
      ! term large, which is the other half of why the ratio is small: the
      ! ratio is a size comparison and not a continuity test.
      !
      ! The face fluxes of the same state do not cancel: the profile is not
      ! a discrete steady state of this operator, and the mass row measures
      ! by how much. That measure is 2e-7, four decades under the fraction a
      ! corrected row is accurate to, so every cell of this column IS
      ! corrected and each correction is good to 2e-7 of itself (row 3); the
      ! same measure is five decades ABOVE cert_tol_mass, so the column is
      ! not a certified stationary state (row 4).
      ! ------------------------------------------------------------------ !
      do j = 1-Ng,N+Ng
         rho_s(j) = 1.0d0
         v_s(j)   = v_ref/(r(j)*r(j))
         Tc_s(j)   = T_base + T_grad*(r(j) - 1.0d0)
         p_s(j)   = nhat*rho_s(j)*Tc_s(j)
      enddo
      call build_state
      call assemble_residual(u, n_part, heat, cool, Res)

      ! The center-to-center test, formed exactly as post_process_adv forms
      ! it for the sensitivity diagnostic.
      ratio_max = 0.0d0
      dlnF_max  = 0.0d0
      do j = ja,jb
         Fup  = rho_s(j)  *v_s(j)  *r(j)**2
         Flo  = rho_s(j-1)*v_s(j-1)*r(j-1)**2
         dVc  = (r(j)**3 - r(j-1)**3)/3.0d0
         divf = (Fup - Flo)/dVc
         drc  = r(j) - r(j-1)
         e_up = Tc_s(j)  /((gamma_ad - 1.0d0)*(1.0d0/nhat))
         e_lo = Tc_s(j-1)/((gamma_ad - 1.0d0)*(1.0d0/nhat))
         w_up = Tc_s(j)*nhat
         h_up = e_up + w_up
         ratio = enthalpy_flux_term_ratio(h_up, divf, drc, rho_s(j),       &
                    v_s(j), e_up, e_lo, w_up, rho_s(j-1))
         ratio_max = max(ratio_max, ratio)
         dlnF_max  = max(dlnF_max, abs(Fup - Flo)/max(abs(Fup),1.0d-99))
      enddo
      call chk_le('enthalpy_ratio_accepts_the_analytic_outflow',           &
                  ratio_max, enthalpy_ratio_report_level, nfail)
      write(*,'(a,es13.6)') '  largest center-to-center |dF/F| of the '//  &
         'same cells: ', dlnF_max

      ! The face-flux mass row of the same cells.
      row_max   = 0.0d0
      r_row_max = 0.0d0
      do j = ja,jb
         row = mass_row_of(j)
         if (row .gt. row_max) then
            row_max   = row
            r_row_max = r(j)
         endif
      enddo
      write(*,'(a,es13.6,a,f8.4,a)') '  largest face-flux mass row '//     &
         '|R_1|/s_1 of the same cells: ', row_max, ' at r = ',             &
         r_row_max, ' Rp'
      call chk_le('analytic_outflow_is_a_conditional_correction', row_max,  &
                  adv_conditional_tol, nfail)

      ! Row 4: and it is NOT a certified stationary state. The same measure
      ! against the other number, so that the file's two statements cannot
      ! be read as one.
      call chk_gt('analytic_outflow_is_not_certified_stationary', row_max, &
                  cert_tol_mass, nfail)

      ! ------------------------------------------------------------------ !
      ! Rows 5 and 6: THE TWO NUMBERS, both stated. adv_conditional_tol is
      ! the fraction of itself a corrected row is accurate to, and it is one
      ! percent; cert_tol_mass is the tolerance of a certified stationary
      ! state, and it is smaller, so the rows of a certified state are
      ! corrected while a corrected row is not thereby certified.
      ! ------------------------------------------------------------------ !
      call chk_exact('conditional_tolerance_is_one_percent',               &
                     adv_conditional_tol - 1.0d-2, nfail)
      call chk_gt('conditional_tolerance_is_above_the_certification',      &
                  adv_conditional_tol, cert_tol_mass, nfail)

      ! ------------------------------------------------------------------ !
      ! Row 7: A COLUMN AT THREE PERCENT IS RETAINED. The same analytic
      ! outflow with its mass flux growing as exp(flux_growth*(r-1)), whose
      ! fractional change across a cell is flux_growth*dr = 3e-2 in every
      ! cell: three times the fraction a corrected row is accurate to, so no
      ! cell of it is corrected. The row asserts the SMALLEST measure of the
      ! interior cells, which is the statement that every one of them is
      ! above the fraction.
      ! ------------------------------------------------------------------ !
      do j = 1-Ng,N+Ng
         rho_s(j) = 1.0d0
         v_s(j)   = v_ref*exp(flux_growth*(r(j) - 1.0d0))/(r(j)*r(j))
         Tc_s(j)  = T_base + T_grad*(r(j) - 1.0d0)
         p_s(j)   = nhat*rho_s(j)*Tc_s(j)
      enddo
      call build_state
      call assemble_residual(u, n_part, heat, cool, Res)
      row_min = huge(1.0d0)
      row_max = 0.0d0
      do j = ja,jb
         row     = mass_row_of(j)
         row_min = min(row_min, row)
         row_max = max(row_max, row)
      enddo
      write(*,'(a,es13.6,a,es13.6)') '  mass row of the growing-flux '//   &
         'column, smallest ', row_min, ' largest ', row_max
      call chk_gt('growing_flux_column_is_retained', row_min,              &
                  adv_conditional_tol, nfail)

      if (nfail .gt. 0) then
         write(*,'(a,i0,a)') 'stationarity_by_the_mass_row: ', nfail,      &
                             ' FAILED'
         call exit(1)
      endif

      contains

      ! ------------------------------------------------------------------ !

      subroutine build_state
      ! The primitive state into conserved variables, the particle count the
      ! residual is given, and the lower boundary the ghost fill needs.
      integer :: jc
      do jc = 1-Ng,N+Ng
         u(1,jc) = rho_s(jc)
         u(2,jc) = rho_s(jc)*v_s(jc)
         u(3,jc) = 0.5d0*rho_s(jc)*v_s(jc)**2                             &
                   + p_s(jc)/(gamma_ad - 1.0d0)
         n_part(jc) = nhat*rho_s(jc)
      enddo
      n_part_cell1 = n_part(1)
      call set_base_reservoir(p_s(1), Tc_s(1), nhat, 1.0d0)
      call Apply_BC(u)
      end subroutine build_state

      ! ------------------------------------------------------------------ !

      real*8 function mass_row_of(jc)
      ! The measure the post-process reads: the cell's mass row over the
      ! largest term that row holds, both from the production routines.
      integer, intent(in) :: jc
      mass_row_of = abs(Res(1,jc))                                          &
                    /max(residual_row_scale(1,jc,u), cert_scale_floor)
      end function mass_row_of

      ! ------------------------------------------------------------------ !

      subroutine chk_le(name, measured, level, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, level
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (measured .le. level)) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),   &
         ' measured=', measured, ' reference=<=', level, ' tol=0'
      end subroutine chk_le

      subroutine chk_gt(name, measured, level, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, level
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (measured .gt. level)) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),   &
         ' measured=', measured, ' reference=>', level, ' tol=0'
      end subroutine chk_gt

      subroutine chk_exact(name, measured, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (measured .ne. 0.0d0) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a)') trim(verdict), trim(name),            &
         ' measured=', measured, ' reference=0 tol=0'
      end subroutine chk_exact

      end program stationarity_by_the_mass_row

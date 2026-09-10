      program conduction_floor_tests
      ! Tests of the Crank-Nicolson transport stage of viscous_conduction
      ! (docs/PLAN_20260906_rev2.md step B2b): what it returns, what it
      ! refuses to return, and what it conserves.
      !
      ! The stage is driven directly, on a small uniform grid with the
      ! conduction key on and the viscosity key off, so the temperature
      ! equation it solves is
      !     C_j (T*_j - T_j)/dt = (1/2)[ Q_j(T) + Q_j(T*) ] ,
      !     Q_j = (1/dV_j)[ F_{j+1/2}(T_{j+1} - T_j)
      !                     - F_{j-1/2}(T_j - T_{j-1}) ] ,
      ! with F the face area times the interpolated conductivity over the
      ! centre spacing, the base ghost held fixed (Dirichlet) and the outer
      ! face carrying zero diffusive flux.  The driver assembles that same
      ! equation from the flux form and solves it with its own elimination,
      ! taking the conductivity from the module's public thermal_conductivity
      ! so that the reference tests the discretization and the solve rather
      ! than the coefficient fit.
      !
      ! Each assertion prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line. Lines beginning with two spaces or with DIAGNOSTIC are context.
      ! Exit status is nonzero if any assertion fails.
      use global_parameters
      use viscous_conduction
      use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      implicit none

      integer :: n_fail

      n_fail = 0
      call build_grid()
      ! The failure statuses are exercised on purpose here, so the module
      ! must report them instead of stopping the program.
      conduction_stop_on_failure = .false.

      call test_above_floor(n_fail)
      call test_floor(n_fail)
      call test_nonfinite_input(n_fail)
      call test_conservation(n_fail)

      write(*,*)
      if (n_fail .gt. 0) then
         write(*,'(a,i0,a)') 'conduction_floor: ', n_fail,                &
                             ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'conduction_floor: every assertion passed'

      contains

      ! ------------------------------------------------------!

      subroutine check(name, measured, reference, tol, nf)
      ! Relative comparison, absolute when the reference is zero.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      integer, intent(inout) :: nf
      real*8 :: d
      if (abs(reference) .gt. 0.0d0) then
         d = abs(measured - reference)/abs(reference)
      else
         d = abs(measured - reference)
      endif
      if (d .le. tol) then
         write(*,'(a,a,a,es14.7,a,es14.7,a,es9.2)') 'PASS ', name,        &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(a,a,a,es14.7,a,es14.7,a,es9.2)') 'FAIL ', name,        &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nf = nf + 1
      endif
      end subroutine check

      subroutine check_int(name, measured, reference, nf)
      character(len=*), intent(in) :: name
      integer, intent(in) :: measured, reference
      integer, intent(inout) :: nf
      if (measured .eq. reference) then
         write(*,'(a,a,a,i0,a,i0,a)') 'PASS ', name, ' measured=',        &
            measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(a,a,a,i0,a,i0,a)') 'FAIL ', name, ' measured=',        &
            measured, ' reference=', reference, ' tol=0'
         nf = nf + 1
      endif
      end subroutine check_int

      subroutine check_bitwise(name, measured, reference, nf)
      ! Equality of every entry of two state arrays, bit for bit.
      character(len=*), intent(in) :: name
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: measured, reference
      integer, intent(inout) :: nf
      integer :: i, j, ndiff
      ndiff = 0
      do j = 1-Ng, N+Ng
         do i = 1, 3
            if (measured(i,j) .ne. reference(i,j)) ndiff = ndiff + 1
         enddo
      enddo
      call check_int(name, ndiff, 0, nf)
      end subroutine check_bitwise

      subroutine check_above(name, measured, threshold, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, threshold
      integer, intent(inout) :: nf
      if (measured .gt. threshold) then
         write(*,'(a,a,a,es14.7,a,es14.7,a)') 'PASS ', name,              &
            ' measured=', measured, ' reference=>', threshold, ' tol=0'
      else
         write(*,'(a,a,a,es14.7,a,es14.7,a)') 'FAIL ', name,              &
            ' measured=', measured, ' reference=>', threshold, ' tol=0'
         nf = nf + 1
      endif
      end subroutine check_above

      ! ------------------------------------------------------!

      subroutine build_grid()
      ! A uniform radial grid of 20 cells between 1 and 2 planetary radii,
      ! and the normalizations the transport coefficients are expressed in.
      integer :: j
      real*8  :: dr
      N  = 20
      T0 = 1.0d3
      n0 = 1.0d13
      R0 = 1.0d10
      v0 = sqrt(kb_erg*T0/mu)
      cond_on  = .true.
      visc_on  = .false.
      visc_mu0 = 0.0d0
      ledger_family = ledger_family_init
      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng))
      dr = 1.0d0/dble(N)
      do j = 1-Ng, N+Ng
         r_edg(j) = 1.0d0 + dble(j)*dr
         r(j)     = 1.0d0 + (dble(j) - 0.5d0)*dr
      enddo
      end subroutine build_grid

      ! ------------------------------------------------------!

      subroutine face_coefficients(Tcell, Fm, Fp, dVol)
      ! F_{j-1/2} = A_{j-1/2} kappa_{j-1/2}/(r_j - r_{j-1}) and its partner
      ! at the right face, from the flux form of the conduction operator.
      ! Zero diffusive flux through the outer face, so F_{N+1/2} = 0.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(N),         intent(out) :: Fm, Fp, dVol
      real*8, dimension(1-Ng:N+Ng) :: kap
      real*8  :: rp, rm, Ap, Am, dpj, dmj, wp, wm, kp, km
      integer :: j
      call thermal_conductivity(Tcell, kap)
      do j = 1, N
         rp = r_edg(j);   rm = r_edg(j-1)
         Ap = rp*rp;      Am = rm*rm
         dVol(j) = (Ap*rp - Am*rm)/3.0d0
         if (j .lt. N) then
            dpj = 1.0d0/(r(j+1) - r(j))
            wp  = (rp - r(j))*dpj
            kp  = kap(j) + wp*(kap(j+1) - kap(j))
            Fp(j) = Ap*kp*dpj
         else
            Fp(j) = 0.0d0
         endif
         dmj = 1.0d0/(r(j) - r(j-1))
         wm  = (rm - r(j-1))*dmj
         km  = kap(j-1) + wm*(kap(j) - kap(j-1))
         Fm(j) = Am*km*dmj
      enddo
      end subroutine face_coefficients

      subroutine reference_operator(Tcell, Q)
      ! Q_j of the conduction operator, at the conductivity of Tcell.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell
      real*8, dimension(N),         intent(out) :: Q
      real*8, dimension(N) :: Fm, Fp, dVol
      integer :: j
      call face_coefficients(Tcell, Fm, Fp, dVol)
      do j = 1, N
         Q(j) = (Fp(j)*(Tcell(j+1) - Tcell(j))                           &
                 - Fm(j)*(Tcell(j) - Tcell(j-1)))/dVol(j)
      enddo
      end subroutine reference_operator

      subroutine reference_step(Tcell, n_part, dt, Tref)
      ! The Crank-Nicolson temperature the stage should return, from an
      ! independent assembly and an independent elimination.
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: Tcell, n_part, dt
      real*8, dimension(N),         intent(out) :: Tref
      real*8, dimension(N) :: Fm, Fp, dVol, Q, alo, adi, aup, rhs
      real*8, dimension(N) :: cp, dpv
      real*8  :: cap, den
      integer :: j
      call face_coefficients(Tcell, Fm, Fp, dVol)
      call reference_operator(Tcell, Q)
      do j = 1, N
         cap    = n_part(j)/(gamma_ad - 1.0d0)/dt(j)
         alo(j) = -0.5d0*Fm(j)/dVol(j)
         aup(j) = -0.5d0*Fp(j)/dVol(j)
         adi(j) = cap + 0.5d0*(Fp(j) + Fm(j))/dVol(j)
         rhs(j) = cap*Tcell(j) + 0.5d0*Q(j)
      enddo
      ! The base ghost is held fixed, so its column is a known term.
      rhs(1) = rhs(1) - alo(1)*Tcell(0)
      alo(1) = 0.0d0
      den    = adi(1)
      cp(1)  = aup(1)/den
      dpv(1) = rhs(1)/den
      do j = 2, N
         den    = adi(j) - alo(j)*cp(j-1)
         cp(j)  = aup(j)/den
         dpv(j) = (rhs(j) - alo(j)*dpv(j-1))/den
      enddo
      Tref(N) = dpv(N)
      do j = N-1, 1, -1
         Tref(j) = dpv(j) - cp(j)*Tref(j+1)
      enddo
      end subroutine reference_step

      ! ------------------------------------------------------!

      subroutine make_state(Tcell, n_part, u, W)
      ! A conserved state at rest with the given temperature and particle
      ! count: p = n_part T, no kinetic energy.
      real*8, dimension(1-Ng:N+Ng),   intent(in)  :: Tcell, n_part
      real*8, dimension(3,1-Ng:N+Ng), intent(out) :: u, W
      integer :: j
      do j = 1-Ng, N+Ng
         W(1,j) = 1.0d0
         W(2,j) = 0.0d0
         W(3,j) = n_part(j)*Tcell(j)
         u(1,j) = W(1,j)
         u(2,j) = 0.0d0
         u(3,j) = W(3,j)/(gamma_ad - 1.0d0)
      enddo
      end subroutine make_state

      ! ------------------------------------------------------!

      subroutine test_above_floor(nf)
      ! (a) A column whose Crank-Nicolson solution stays above the floor.
      ! The stage returns OK and the temperature of the independent solve.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng)   :: Tcell, n_part, dt
      real*8, dimension(3,1-Ng:N+Ng) :: u, W
      real*8, dimension(N) :: Tref
      real*8  :: worst, d
      integer :: j, stat, hits_before
      write(*,'(a)') '---- above_floor ----'
      do j = 1-Ng, N+Ng
         Tcell(j)  = 1.0d0 + 0.5d0*exp(-(r(j) - 1.0d0)/0.2d0)
         n_part(j) = 1.0d0
         dt(j)     = 1.0d5
      enddo
      ! The Dirichlet anchor of the base face, a ghost at the base value.
      Tcell(0) = Tcell(1)
      call make_state(Tcell, n_part, u, W)
      call reference_step(Tcell, n_part, dt, Tref)
      hits_before = n_conduction_floor_hits
      call viscous_conduction_step(u, W, Tcell, n_part, dt, 7, stat)
      call check_int('above_floor_status', stat, CONDUCTION_OK, nf)
      call check_int('above_floor_no_floor_hit',                          &
                     n_conduction_floor_hits - hits_before, 0, nf)
      worst = 0.0d0
      do j = 1, N
         d = abs(W(3,j)/n_part(j) - Tref(j))/abs(Tref(j))
         if (d .gt. worst) worst = d
      enddo
      call check('above_floor_temperature_vs_reference', worst, 0.0d0,    &
                 1.0d-12, nf)
      ! Every returned temperature is above the floor, which is what makes
      ! this the admissible branch.
      worst = 1.0d30
      do j = 1, N
         if (W(3,j)/n_part(j) .lt. worst) worst = W(3,j)/n_part(j)
      enddo
      call check_above('above_floor_min_temperature_over_floor',          &
                       worst, conduction_temperature_floor(), nf)
      ! The energy row follows the pressure row it was rebuilt from.
      worst = 0.0d0
      do j = 1, N
         d = abs(u(3,j) - W(3,j)/(gamma_ad - 1.0d0))/abs(u(3,j))
         if (d .gt. worst) worst = d
      enddo
      call check('above_floor_energy_row', worst, 0.0d0, 1.0d-15, nf)
      end subroutine test_above_floor

      ! ------------------------------------------------------!

      subroutine test_floor(nf)
      ! (b) A column drained against a cold base anchor, so that the
      ! implicit operator asks the base cell for a temperature below the
      ! floor. The stage must fail, name the cell, report the energy a
      ! floor would have injected, and return no state.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng)   :: Tcell, n_part, dt
      real*8, dimension(3,1-Ng:N+Ng) :: u, W, u_in, W_in
      real*8, dimension(N) :: Tref
      real*8  :: T_flr, deficit_ref, clamped
      integer :: j, stat, hits_before, n_below
      write(*,'(a)') '---- floor ----'
      do j = 1-Ng, N+Ng
         Tcell(j)  = 1.0d0
         n_part(j) = 1.0d0
         dt(j)     = 1.0d10
      enddo
      ! A base anchor at zero temperature: the base cell is drained toward
      ! it and, over this step, past the floor.
      Tcell(0) = 0.0d0
      call make_state(Tcell, n_part, u, W)
      u_in = u;  W_in = W
      call reference_step(Tcell, n_part, dt, Tref)
      T_flr = conduction_temperature_floor()
      hits_before = n_conduction_floor_hits
      call viscous_conduction_step(u, W, Tcell, n_part, dt, 11, stat)
      call check_int('floor_status', stat, CONDUCTION_FLOOR, nf)
      call check_int('floor_reason', conduction_last_reason,              &
                     CONDUCTION_REASON_FLOOR_REACHED, nf)
      call check_int('floor_cell', conduction_last_cell, 1, nf)
      ! Every cell the operator asks to go below the floor is one refused
      ! attempt, and the counters keep counting attempts.
      n_below = 0
      do j = 1, N
         if (Tref(j) .le. T_flr) n_below = n_below + 1
      enddo
      call check_int('floor_counted_as_an_attempt',                       &
                     n_conduction_floor_hits - hits_before, n_below, nf)
      call check_int('floor_failing_cell_count', conduction_last_nfail,   &
                     n_below, nf)
      ! THE CELL-BY-CELL RECORD, read the way energy_semi_implicit's
      ! energy_floor_cell_hits is read.  It is what the stage produced on
      ! this constructed failure: one count in every cell the operator
      ! asked to go below the floor, and none anywhere else.
      call check_int('floor_cell_record_is_allocated',                    &
                     merge(1, 0, allocated(conduction_floor_cell_hits)),  &
                     1, nf)
      if (allocated(conduction_floor_cell_hits)) then
         call check_int('floor_cell_record_sums_to_the_attempts',         &
                        sum(conduction_floor_cell_hits(1:N)) -            &
                        hits_before, n_below, nf)
         call check_int('floor_cell_record_marks_the_cells_below_the'//   &
                        '_floor', cells_marked_below(Tref, T_flr),        &
                        n_below, nf)
         call check_int('floor_cell_record_is_the_distinct_cell_count',   &
                        n_conduction_floor_cells(), n_below, nf)
      endif
      call check('floor_reported_solution', conduction_last_T_solution,   &
                 Tref(1), 1.0d-12, nf)
      ! c_v (T_floor - T_solved) per volume, in erg/cm3: the code energy
      ! density unit is n0 k_B T0.
      deficit_ref = (n_part(1)/(gamma_ad - 1.0d0))*(T_flr - Tref(1))      &
                    *n0*kb_erg*T0
      call check('floor_energy_deficit', conduction_last_energy_deficit,  &
                 deficit_ref, 1.0d-12, nf)
      call check_above('floor_energy_deficit_is_positive',                &
                       conduction_last_energy_deficit, 0.0d0, nf)
      ! No state is returned: the caller still holds what it passed in.
      call check_bitwise('floor_state_unchanged_W', W, W_in, nf)
      call check_bitwise('floor_state_unchanged_u', u, u_in, nf)
      ! The standing red measurement: what the clamp this step removed would
      ! have returned instead, and how far that is from the solved value.
      clamped = max(Tref(1), T_flr)
      call check_above('floor_clamp_would_have_moved_the_state_red',      &
                       abs(clamped - Tref(1))/T_flr, 1.0d-6, nf)
      write(*,'(a,es13.6,a,es13.6)') '  DIAGNOSTIC solved T [K] ',        &
         Tref(1)*T0, '   clamped T [K] ', clamped*T0
      end subroutine test_floor

      ! ------------------------------------------------------!

      integer function cells_marked_below(Tref, T_flr) result(n)
      ! How many of the cells the reference step puts at or below the floor
      ! carry a mark in the cell-by-cell record.  Equal to the number of
      ! such cells exactly when the record marks all of them and only them.
      real*8, intent(in) :: Tref(:), T_flr
      integer :: j
      n = 0
      do j = 1, size(Tref)
         if (Tref(j) .le. T_flr) then
            if (conduction_floor_cell_hits(j) .gt. 0) n = n + 1
         else
            if (conduction_floor_cell_hits(j) .gt. 0) n = n - 1
         endif
      enddo
      end function cells_marked_below

      ! ------------------------------------------------------!

      subroutine test_nonfinite_input(nf)
      ! (c) A non-finite temperature on entry is reported, not propagated.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng)   :: Tcell, n_part, dt
      real*8, dimension(3,1-Ng:N+Ng) :: u, W, u_in, W_in
      integer :: j, stat
      write(*,'(a)') '---- nonfinite_input ----'
      do j = 1-Ng, N+Ng
         Tcell(j)  = 1.0d0
         n_part(j) = 1.0d0
         dt(j)     = 1.0d5
      enddo
      call make_state(Tcell, n_part, u, W)
      Tcell(3) = ieee_value(1.0d0, ieee_quiet_nan)
      u_in = u;  W_in = W
      call viscous_conduction_step(u, W, Tcell, n_part, dt, 13, stat)
      call check_int('nonfinite_status', stat, CONDUCTION_NONFINITE, nf)
      call check_int('nonfinite_reason', conduction_last_reason,          &
                     CONDUCTION_REASON_NONFINITE_INPUT, nf)
      call check_int('nonfinite_cell', conduction_last_cell, 3, nf)
      call check_bitwise('nonfinite_state_unchanged_W', W, W_in, nf)
      call check_bitwise('nonfinite_state_unchanged_u', u, u_in, nf)
      end subroutine test_nonfinite_input

      ! ------------------------------------------------------!

      subroutine test_conservation(nf)
      ! (d) Conservation.
      !   (d1) The operator is in flux form, so on a column with no flux
      !        through either end the total thermal energy it moves is zero
      !        to round-off, whatever the interior profile is.
      !   (d2) The stage as boundary-conditioned in a run is NOT closed at
      !        the base: the ghost is a Dirichlet anchor. What it conserves
      !        is the total thermal energy up to the heat that crosses that
      !        one face, and that identity is asserted here.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng)   :: Tcell, n_part, dt
      real*8, dimension(3,1-Ng:N+Ng) :: u, W
      real*8, dimension(N) :: Q, Fm, Fp, dVol
      real*8  :: total, scale, dE, influx, Tnew1
      integer :: j, stat
      write(*,'(a)') '---- conservation ----'
      do j = 1-Ng, N+Ng
         Tcell(j)  = 1.0d0 + 0.3d0*sin(6.0d0*(r(j) - 1.0d0))
         n_part(j) = 1.0d0 + 0.5d0*(r(j) - 1.0d0)
         dt(j)     = 1.0d5
      enddo
      ! Zero flux through the base face as well: the ghost carries the base
      ! cell's own temperature, so the column is closed at both ends.
      Tcell(0) = Tcell(1)

      call reference_operator(Tcell, Q)
      call face_coefficients(Tcell, Fm, Fp, dVol)
      total = 0.0d0
      scale = 0.0d0
      do j = 1, N
         total = total + dVol(j)*Q(j)
         scale = scale + abs(dVol(j)*Q(j))
      enddo
      call check('closed_column_operator_conserves', total/scale, 0.0d0,  &
                 1.0d-13, nf)

      call make_state(Tcell, n_part, u, W)
      call viscous_conduction_step(u, W, Tcell, n_part, dt, 17, stat)
      call check_int('conservation_status', stat, CONDUCTION_OK, nf)
      dE     = 0.0d0
      scale  = 0.0d0
      do j = 1, N
         dE = dE + dVol(j)*(n_part(j)/(gamma_ad - 1.0d0))                 &
                   *(W(3,j)/n_part(j) - Tcell(j))
         scale = scale + dVol(j)*(n_part(j)/(gamma_ad - 1.0d0))           &
                   *abs(W(3,j)/n_part(j) - Tcell(j))
      enddo
      ! The heat that entered through the base face over the step, at the
      ! Crank-Nicolson average of the two times, with the ghost held fixed.
      Tnew1  = W(3,1)/n_part(1)
      influx = dt(1)*0.5d0*Fm(1)*((Tcell(0) - Tcell(1))                   &
                                  + (Tcell(0) - Tnew1))
      call check('stage_conserves_up_to_the_base_face_flux',              &
                 (dE - influx)/scale, 0.0d0, 1.0d-12, nf)
      write(*,'(a,es13.6,a,es13.6)') '  DIAGNOSTIC thermal energy change ',&
         dE, '   base face influx ', influx
      end subroutine test_conservation

      end program conduction_floor_tests

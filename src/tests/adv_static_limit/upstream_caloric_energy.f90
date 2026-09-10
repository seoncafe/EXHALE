      program upstream_caloric_energy
      ! THE UPSTREAM ENERGY OF THE ADVECTED ENERGY EQUATION (T_equation.f90).
      !
      ! The steady internal-energy equation of a supplied profile is
      !
      !     rho v de/dr - p v dln(rho)/dr + h div(rho v) = heating - cooling
      !
      ! and its first term is an upwind difference of the SPECIFIC internal
      ! energy e = E(x_H2,T)/mu, E being the energy per particle of the
      ! caloric EOS. What the flow carries into a cell is the specific energy
      ! OF THE UPSTREAM GAS, so the upstream end of that difference is
      ! e_up = E(x_H2,up, T_up)/mu_up and not this cell's caloric state
      ! evaluated at the upstream temperature. The rovibrational heat
      ! capacity of H2 belongs to the gas that holds the molecules: across a
      ! dissociation front two cells at the same temperature store different
      ! energy, and an atomic cell fed by a molecular neighbor receives
      ! molecular energy although it has no H2 of its own.
      !
      ! The rows below pin that on the residual itself.
      !
      !   U0  the atomic cell keeps the constant gamma_ad arithmetic exactly,
      !       so a run with no molecules is unaffected to the bit.
      !   U1  an isothermal column with an H2 gradient and a stationary mass
      !       flux is a steady state of the equation, so the residual
      !       vanishes on it.
      !   U2  the same with the CURRENT cell atomic and the upstream cell
      !       molecular, the case a residual that reads only this cell's
      !       composition cannot represent at all.
      !   U3  a cell whose two sides carry the same composition, the control:
      !       there the upstream composition is this cell's and the residual
      !       is unchanged.
      !   U4  a stationary flow across a dissociation layer, x_H2 = 0.8 down
      !       to 0: the temperature the solver returns against the analytic
      !       solution of the same discrete equation.
      !
      ! WHY U1 AND U2 CAN BE ISOTHERMAL. With a stationary mass flux
      ! (div(rho v) = 0), a constant density (no pressure work) and no net
      ! radiative term, the equation reduces to de/dr = 0: the steady state is
      ! the column of UNIFORM SPECIFIC ENERGY, which is isothermal exactly
      ! when E(x_H2,T)/mu is the same on both sides. That condition fixes the
      ! ratio mu_up/mu_j = E(x_H2,up,T)/E(x_H2,j,T) of the two mean molecular
      ! weights, and the fixtures below use it. Composition alone does not
      ! make a column isothermal, and that is the point: the equation must
      ! vanish on its own steady state whatever the two compositions are.
      !
      ! One line per assertion:
      !   PASS|FAIL <name> measured= reference= tol=
      use global_parameters, only: T0, q0, n0, gamma_ad, thereis_He
      use ion_cell_state, only: teq_cell
      use caloric_eos, only: internal_energy_of_mixture,                   &
                             temperature_of_mixture
      use equation_T, only: Tres, solve_T_brent
      implicit none

      real*8  :: params(40)
      real*8  :: T_iso, x_j, x_up, res, legacy, res_atomic
      real*8  :: worst, T_num, T_exact, d
      real*8  :: x_front(5), T_march(5), T_disc(5)
      real*8  :: res_wrong_e_up
      logical :: ok
      integer :: nfail, j

      nfail  = 0
      params = 0.0d0

      ! Scales of the residual. q0 divides the cooling sum, so it has to be
      ! nonzero even where every density is; T0 sets the temperature unit.
      T0 = 1.0d3
      q0 = 1.0d0
      n0 = 1.0d0
      thereis_He = .false.

      ! ---------------------------------------------------------------!
      ! Row U0: the atomic cell is the constant gamma_ad arithmetic.
      !
      ! With no H2 on either side of the interface the residual is the
      ! monatomic form written out in full, and the row asserts EXACT
      ! equality with that expression: an atomic run must reproduce it to
      ! the bit, so the branch is structural and not a limit of the mixture
      ! formula (the same gate caloric_eos states for its own maps).
      ! ---------------------------------------------------------------!
      call one_interface(0.0d0, 0.0d0, 1.2d0, 1.3d0, 1.0d0)
      teq_cell%coeff    = -0.4d0
      teq_cell%div_rhov = 0.37d0
      res_atomic = Tres(1.1d0, params)
      legacy = teq_cell%mum*teq_cell%rhov*1.1d0                            &
             - teq_cell%mup*teq_cell%rhov*teq_cell%Told                    &
             + gamma_ad*teq_cell%mum*teq_cell%dr*teq_cell%div_rhov*1.1d0   &
             - (gamma_ad - 1.0d0)*(teq_cell%coeff*1.1d0                    &
               + teq_cell%mup*teq_cell%mum*teq_cell%dr*teq_cell%heaold)
      call chk_bitwise('atomic_interface_is_the_monatomic_arithmetic',      &
                       res_atomic, legacy, nfail)
      ! And it does not read the upstream caloric energy at all, whatever is
      ! in that field: the guarantee that a run with no molecules keeps the
      ! arithmetic it had is a property of the branch, not of the value.
      teq_cell%e_up  = 1.0d30
      res_wrong_e_up = Tres(1.1d0, params)
      call chk_bitwise('atomic_interface_ignores_the_upstream_energy',      &
                       res_wrong_e_up, res_atomic, nfail)

      ! ---------------------------------------------------------------!
      ! Row U1: an isothermal column with an H2 gradient carries no energy
      ! change, so the residual vanishes on it.
      !
      ! The old residual evaluated the upstream energy at THIS cell's H2
      ! fraction, which on this fixture leaves
      !     rhov*E(x_j,T)*(mum - mup)
      ! behind: the whole molecular part of the energy the flow brings in is
      ! replaced by the molecular part this cell would have at the upstream
      ! temperature. That number is printed as the reference the row is
      ! measured against.
      ! ---------------------------------------------------------------!
      T_iso = 1.0d0                        ! 1000 K in units of T0
      x_j   = 0.2d0
      x_up  = 0.8d0
      call isothermal_interface(x_j, x_up, T_iso)
      res = Tres(T_iso, params)
      write(*,'(a,es13.6)') '  the residual of the same fixture with the'// &
         ' upstream energy at this cell''s H2 fraction: ',                  &
         teq_cell%rhov*internal_energy_of_mixture(x_j, T_iso)              &
         *(teq_cell%mum - teq_cell%mup)
      call chk('isothermal_h2_gradient_carries_no_energy_change',           &
               res, 0.0d0, 1.0d-13, nfail)

      ! ---------------------------------------------------------------!
      ! Row U2: the current cell atomic, the upstream cell molecular. This
      ! is the interface a residual built from this cell's composition alone
      ! cannot represent: it has no H2 to read, so it charges the upstream
      ! gas the monatomic energy and loses the rovibrational part entirely.
      ! ---------------------------------------------------------------!
      call isothermal_interface(0.0d0, 0.6d0, T_iso)
      res = Tres(T_iso, params)
      ! An H2 fraction of zero in the cell selected the monatomic form, whose
      ! own scaling carries the factor (gamma_ad - 1), so what a residual
      ! reading only this cell's composition leaves behind here is
      ! rhov*T*(mum - mup).
      write(*,'(a,es13.6)') '  the residual of the same fixture with the'// &
         ' upstream energy at this cell''s H2 fraction: ',                  &
         teq_cell%rhov*T_iso*(teq_cell%mum - teq_cell%mup)
      call chk('atomic_cell_below_a_molecular_upstream',                    &
               res, 0.0d0, 1.0d-13, nfail)

      ! ---------------------------------------------------------------!
      ! Row U3: the control. Where the two sides carry the same H2 fraction
      ! the upstream composition IS this cell's, the two mean molecular
      ! weights are equal, and the residual of the isothermal column
      ! vanishes with either reading of the upstream energy.
      ! ---------------------------------------------------------------!
      call isothermal_interface(0.4d0, 0.4d0, T_iso)
      res = Tres(T_iso, params)
      call chk('uniform_composition_interface_is_unchanged',                &
               res, 0.0d0, 1.0d-13, nfail)

      ! ---------------------------------------------------------------!
      ! Row U4: a stationary flow across a dissociation layer.
      !
      ! The column is pure hydrogen, so its mean molecular weight per
      ! particle is mu = 1 + x_H2 (a fraction x_H2 of the particles are two
      ! nuclei each and the rest are one). Density is constant and the mass
      ! flux stationary, so the discrete equation of each cell is
      ! e_j = e_up exactly, whose solution is
      !     T_j = temperature_of_mixture(x_j, mu_j*e_up) ,
      ! the inverse of the same caloric EOS. The row marches the layer with
      ! the solver the post-process uses and measures the departure from
      ! that solution.
      ! ---------------------------------------------------------------!
      x_front = (/ 0.8d0, 0.6d0, 0.4d0, 0.1d0, 0.0d0 /)
      T_march(1) = 1.0d0                   ! 1000 K at the foot of the layer
      T_disc(1)  = T_march(1)
      worst = 0.0d0
      do j = 2, 5
         call one_interface(x_front(j), x_front(j-1),                       &
                            1.0d0 + x_front(j), 1.0d0 + x_front(j-1),       &
                            T_march(j-1))
         T_exact = temperature_of_mixture(x_front(j),                       &
                      teq_cell%mup*teq_cell%e_up)
         call solve_T_brent(params, T_march(j-1), T_num, ok)
         if (.not. ok) then
            call chk_lt('dissociation_layer_root_is_bracketed',            &
                        1.0d0, 0.0d0, nfail)
            exit
         endif
         T_march(j) = T_num
         T_disc(j)  = T_exact
         d = abs(T_num - T_exact)/T_exact
         if (d .gt. worst) worst = d
      enddo
      write(*,'(a,5(f9.2,1x))') '  marched T [K] across the layer: ',       &
         T_march*T0
      write(*,'(a,5(f9.2,1x))') '  discrete solution     [K]:      ',      &
         T_disc*T0
      call chk_lt('dissociation_layer_matches_the_discrete_solution',       &
                  worst, 1.0d-8, nfail)

      if (nfail .gt. 0) stop 1

      contains

      ! One interface of a wind, in the units of the residual: the cell's H2
      ! fraction and mean molecular weight, the upstream cell's, and the
      ! upstream temperature. No heating, no cooling and no densities, so the
      ! residual is the advected balance alone; the density is constant, so
      ! the pressure-work coefficient is zero and the mass flux is
      ! stationary, so the enthalpy flux term is absent.
      subroutine one_interface(x_cell, x_upstream, mu_cell, mu_upstream,   &
                               T_upstream)
      real*8, intent(in) :: x_cell, x_upstream, mu_cell, mu_upstream
      real*8, intent(in) :: T_upstream
      teq_cell%nhi    = 0.0d0
      teq_cell%nhii   = 0.0d0
      teq_cell%nhei   = 0.0d0
      teq_cell%nheii  = 0.0d0
      teq_cell%nheiii = 0.0d0
      teq_cell%mup    = mu_cell
      teq_cell%mum    = mu_upstream
      teq_cell%rhov   = 2.5d0
      teq_cell%coeff  = 0.0d0
      teq_cell%dr     = 3.0d-3
      teq_cell%Told   = T_upstream
      teq_cell%heaold = 0.0d0
      teq_cell%x_h2   = x_cell
      teq_cell%x_h2_up = x_upstream
      teq_cell%e_up   = internal_energy_of_mixture(x_upstream, T_upstream) &
                        /mu_upstream
      teq_cell%div_rhov = 0.0d0
      end subroutine one_interface

      ! The same interface with the two mean molecular weights in the ratio
      ! that makes the specific energy uniform at the temperature T, so that
      ! the isothermal column IS the steady state of the equation (see the
      ! statement at the head of this file). mu of the cell is 1; only the
      ! ratio enters.
      subroutine isothermal_interface(x_cell, x_upstream, T)
      real*8, intent(in) :: x_cell, x_upstream, T
      real*8 :: mu_upstream
      mu_upstream = internal_energy_of_mixture(x_upstream, T)              &
                    /internal_energy_of_mixture(x_cell, T)
      call one_interface(x_cell, x_upstream, 1.0d0, mu_upstream, T)
      end subroutine isothermal_interface

      subroutine chk(name, measured, reference, tol, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (abs(measured - reference)                                 &
                 .le. tol*max(abs(reference),1.0d0))) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a,es9.2)') trim(verdict),         &
         trim(name), ' measured=', measured, ' reference=', reference,     &
         ' tol=', tol
      end subroutine chk

      subroutine chk_bitwise(name, measured, reference, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      ! The equality test on a real is the statement being made: the two
      ! expressions must be the SAME double, not merely close.
      if (measured .ne. reference) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es24.17,a,es24.17,a)') trim(verdict), trim(name), &
         ' measured=', measured, ' reference=', reference, ' tol=0'
      end subroutine chk_bitwise

      subroutine chk_lt(name, measured, bound, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, bound
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (measured .lt. bound)) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),   &
         ' measured=', measured, ' reference=<', bound, ' tol=0'
      end subroutine chk_lt

      end program upstream_caloric_energy

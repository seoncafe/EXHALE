      module energy_update_cooling_models
      ! Analytic cooling laws for the tests of the implicit energy update.
      ! Each is a volumetric cooling rate in the units `heat` and `cool`
      ! carry in solve_energy_semi_implicit, as a function of the code
      ! temperature, so the temperature equation the solver is handed has a
      ! closed form or an independently integrable reference.
      !
      ! The shapes are chosen for what they do to the iteration, not to
      ! imitate a particular coolant:
      !   linear      cool = c (T - T_eq); the implicit step has a closed
      !               form, so the solve can be checked against an exact
      !               number rather than against another iteration.
      !   falling     cool = Cp exp(-(T - Tp)/w), a cooling that DECREASES
      !               with temperature, as metal line cooling does above its
      !               peak near 2e4 K. This is the branch on which the
      !               previous update replaced the exact Newton derivative
      !               by 1 + c |dC/dT| and iterated exactly twice.
      !   constant    cool = C0, used to drive the balance below the floor by
      !               a known amount.
      use global_parameters
      implicit none

      ! linear
      real*8 :: lin_c   = 0.0d0
      real*8 :: lin_Teq = 0.0d0
      ! falling
      real*8 :: fall_Cp = 0.0d0
      real*8 :: fall_Tp = 0.0d0
      real*8 :: fall_w  = 1.0d0
      ! constant
      real*8 :: const_C0 = 0.0d0

      contains

      real*8 function cool_linear_scalar(T)
      real*8, intent(in) :: T
      cool_linear_scalar = lin_c*(T - lin_Teq)
      end function cool_linear_scalar

      real*8 function cool_falling_scalar(T)
      real*8, intent(in) :: T
      cool_falling_scalar = fall_Cp*exp(-(T - fall_Tp)/fall_w)
      end function cool_falling_scalar

      subroutine cool_linear(T, C)
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
      real*8, dimension(1-Ng:N+Ng), intent(out) :: C
      integer :: j
      do j = 1-Ng, N+Ng
         C(j) = cool_linear_scalar(T(j))
      enddo
      end subroutine cool_linear

      subroutine cool_falling(T, C)
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
      real*8, dimension(1-Ng:N+Ng), intent(out) :: C
      integer :: j
      do j = 1-Ng, N+Ng
         C(j) = cool_falling_scalar(T(j))
      enddo
      end subroutine cool_falling

      subroutine cool_constant(T, C)
      real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
      real*8, dimension(1-Ng:N+Ng), intent(out) :: C
      integer :: j
      do j = 1-Ng, N+Ng
         C(j) = const_C0
      enddo
      end subroutine cool_constant

      end module energy_update_cooling_models

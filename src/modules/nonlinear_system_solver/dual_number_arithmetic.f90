      ! FORWARD-MODE DIFFERENTIATION BY DUAL NUMBERS.
      !
      ! A dual number holds a value v and its first partial derivatives
      ! d(1:n_dual_partials) with respect to a fixed set of independent
      ! variables.  Every arithmetic operation defined here propagates the
      ! derivatives by the chain rule of that operation, so a routine
      ! written once and evaluated on dual arguments returns its value AND
      ! its exact first derivatives (to the rounding of the arithmetic, with
      ! no difference step and no truncation error).
      !
      ! WHAT IT IS USED FOR.  The carrier transport solve
      ! (diffusive_photochemistry::solve_carriers) needs the derivative of
      ! the chemical source of every carrier with respect to every carrier
      ! density.  The rows of that source are written once, in include
      ! files, and instantiated both on real*8 arguments (the residual and
      ! the local ionization sweep) and on dual arguments (the Jacobian), so
      ! the residual and its Jacobian are the same expressions.  The
      ! independent variables are the carrier densities, and
      ! n_dual_partials is the size of that set (n_carrier_max of
      ! diffusive_photochemistry, which checks the equality at compile
      ! time).
      !
      ! THE ONE NON-SMOOTH OPERATION, positive_part(x) = max(x, 0): of a dual
      ! it returns x itself (value and derivatives) when x%v >= 0, and zero
      ! with zero derivatives otherwise, i.e. the one-sided derivative of
      ! the branch the value is taken from.  At x%v = 0 no derivative
      ! exists; the branch of x is returned, a SELECTED one-sided
      ! derivative: the limit from the side x >= 0, which for the clips of
      ! the carrier source (nonnegative remainders and densities,
      ! diffusive_photochemistry) is the admissible side
      ! (src/tests/carrier_source_kinks).  It is a named function
      ! and not an extension of the intrinsic max, so that no module that
      ! uses this one acquires a generic max (an OpenMP reduction(max:...)
      ! in such a module would name the generic and not the intrinsic).
      module dual_number_arithmetic
      implicit none
      private

      integer, parameter, public :: n_dual_partials = 11

      type, public :: dual_number
         real*8 :: v
         real*8 :: d(n_dual_partials)
      end type dual_number

      public :: operator(+), operator(-), operator(*), operator(/)
      public :: assignment(=), positive_part, value_of

      interface operator(+)
         module procedure dual_add_dual, dual_add_real, real_add_dual,     &
                          dual_unary_plus
      end interface
      interface operator(-)
         module procedure dual_sub_dual, dual_sub_real, real_sub_dual,     &
                          dual_negate
      end interface
      interface operator(*)
         module procedure dual_mul_dual, dual_mul_real, real_mul_dual
      end interface
      interface operator(/)
         module procedure dual_div_real
      end interface
      interface assignment(=)
         module procedure dual_from_real
      end interface
      ! max(x, 0) of a number of either kind (one source text for both).
      interface positive_part
         module procedure positive_part_dual, positive_part_real
      end interface
      ! The value part of a number of either kind, so that one source text
      ! can test a value whatever kind it was instantiated at.
      interface value_of
         module procedure value_of_dual, value_of_real
      end interface

      contains

      elemental subroutine dual_from_real(a, r)
      type(dual_number), intent(out) :: a
      real*8,            intent(in)  :: r
      a%v = r
      a%d = 0.0d0
      end subroutine dual_from_real

      elemental function dual_add_dual(a, b) result(c)
      type(dual_number), intent(in) :: a, b
      type(dual_number) :: c
      c%v = a%v + b%v
      c%d = a%d + b%d
      end function dual_add_dual

      elemental function dual_add_real(a, r) result(c)
      type(dual_number), intent(in) :: a
      real*8,            intent(in) :: r
      type(dual_number) :: c
      c%v = a%v + r
      c%d = a%d
      end function dual_add_real

      elemental function real_add_dual(r, a) result(c)
      real*8,            intent(in) :: r
      type(dual_number), intent(in) :: a
      type(dual_number) :: c
      c%v = r + a%v
      c%d = a%d
      end function real_add_dual

      elemental function dual_unary_plus(a) result(c)
      type(dual_number), intent(in) :: a
      type(dual_number) :: c
      c = a
      end function dual_unary_plus

      elemental function dual_sub_dual(a, b) result(c)
      type(dual_number), intent(in) :: a, b
      type(dual_number) :: c
      c%v = a%v - b%v
      c%d = a%d - b%d
      end function dual_sub_dual

      elemental function dual_sub_real(a, r) result(c)
      type(dual_number), intent(in) :: a
      real*8,            intent(in) :: r
      type(dual_number) :: c
      c%v = a%v - r
      c%d = a%d
      end function dual_sub_real

      elemental function real_sub_dual(r, a) result(c)
      real*8,            intent(in) :: r
      type(dual_number), intent(in) :: a
      type(dual_number) :: c
      c%v = r - a%v
      c%d = -a%d
      end function real_sub_dual

      elemental function dual_negate(a) result(c)
      type(dual_number), intent(in) :: a
      type(dual_number) :: c
      c%v = -a%v
      c%d = -a%d
      end function dual_negate

      elemental function dual_mul_dual(a, b) result(c)
      type(dual_number), intent(in) :: a, b
      type(dual_number) :: c
      c%v = a%v*b%v
      c%d = a%d*b%v + a%v*b%d
      end function dual_mul_dual

      elemental function dual_mul_real(a, r) result(c)
      type(dual_number), intent(in) :: a
      real*8,            intent(in) :: r
      type(dual_number) :: c
      c%v = a%v*r
      c%d = a%d*r
      end function dual_mul_real

      elemental function real_mul_dual(r, a) result(c)
      real*8,            intent(in) :: r
      type(dual_number), intent(in) :: a
      type(dual_number) :: c
      c%v = r*a%v
      c%d = r*a%d
      end function real_mul_dual

      elemental function dual_div_real(a, r) result(c)
      type(dual_number), intent(in) :: a
      real*8,            intent(in) :: r
      type(dual_number) :: c
      c%v = a%v/r
      c%d = a%d/r
      end function dual_div_real

      elemental function positive_part_dual(a) result(c)
      type(dual_number), intent(in) :: a
      type(dual_number) :: c
      if (a%v .ge. 0.0d0) then
         c = a
      else
         c%v = 0.0d0
         c%d = 0.0d0
      endif
      end function positive_part_dual

      elemental real*8 function positive_part_real(r) result(c)
      real*8, intent(in) :: r
      c = max(r, 0.0d0)
      end function positive_part_real

      elemental real*8 function value_of_dual(a) result(v)
      type(dual_number), intent(in) :: a
      v = a%v
      end function value_of_dual

      elemental real*8 function value_of_real(r) result(v)
      real*8, intent(in) :: r
      v = r
      end function value_of_real

      end module dual_number_arithmetic
